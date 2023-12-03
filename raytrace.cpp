#include <fstream>
#include <sstream>
#include <tuple>
#include "raytrace.h"
#include "ray.h"
#include "matrix_ops.h"
#include "file_processing.h"

bool shadowRay(Intersection& intersection, Light& light, const Sphere& sphere)
{
    // ray from intersection point to light
    Ray shadow_ray = Ray(intersection.point, glm::normalize(glm::vec3(light.posX, light.posY, light.posZ) - intersection.point));

    glm::mat4 inverse_translate = inverseTranslate(sphere.posX, sphere.posY, sphere.posZ);
    glm::mat4 inverse_scale = inverseScale(sphere.scaleX, sphere.scaleY, sphere.scaleZ);

    glm::vec4 S_homo = glm::vec4(shadow_ray.get_origin(), 1);
    glm::vec4 c_homo = glm::vec4(shadow_ray.get_direction(), 0);
    glm::vec4 S_t_homo = inverse_scale * inverse_translate * S_homo;
    glm::vec4 c_t_homo = inverse_scale * inverse_translate * c_homo;
    glm::vec3 S_t = glm::vec3(S_t_homo.x, S_t_homo.y, S_t_homo.z);
    glm::vec3 c_t = glm::vec3(c_t_homo.x, c_t_homo.y, c_t_homo.z);
    Ray shadow_ray_t = Ray(S_t, c_t);

    float A = glm::dot(c_t, c_t);
    float B = glm::dot(c_t, S_t);
    float C = glm::dot(S_t, S_t) - 1;

    float discriminant = B * B - A * C;

    if (discriminant > 0)
    {
        float t1 = ((-B/A) + (sqrt(discriminant)/A));
        float t2 = ((-B/A) - (sqrt(discriminant)/A));

        if(intersection.concave_flag && (intersection.point.z > light.posZ))
        {
            return true;
        }

        if( (t1 > 0.0001 && (shadow_ray.at(t1).z) < -1 ) || (t2 > 0.0001 && (shadow_ray.at(t2).z < -1) ))
        {
            return true;
        }
    }
    return false;
}

bool isLightBehind(Intersection& intersection, Light& light)
{
    if(intersection.concave_flag == 1 && intersection.point.z > light.posZ)
    {
        return true;
    }
    return false;
}

void closestIntersection(Ray& r, Intersection& intersection, const ImageInfo &inputImage)
{
    glm::vec3 S = r.get_origin();
    glm::vec3 c = r.get_direction();

    intersection.no_intersect_flag = 1;

    for(int i = 0; i < inputImage.spheres.size(); i++)
    {
        Sphere sphere = inputImage.spheres[i];

        glm::mat4 inverse_translate = inverseTranslate(sphere.posX, sphere.posY, sphere.posZ);
        glm::mat4 inverse_scale = inverseScale(sphere.scaleX, sphere.scaleY, sphere.scaleZ);

        glm::vec4 S_homo = glm::vec4(r.get_origin(), 1);
        glm::vec4 c_homo = glm::vec4(r.get_direction(), 0);
        glm::vec4 S_t_homo = inverse_scale * inverse_translate * S_homo;
        glm::vec4 c_t_homo = inverse_scale * inverse_translate * c_homo;
        glm::vec3 S_t = glm::vec3(S_t_homo.x, S_t_homo.y, S_t_homo.z);
        glm::vec3 c_t = glm::vec3(c_t_homo.x, c_t_homo.y, c_t_homo.z);
        Ray r_t = Ray(S_t, c_t);

        float A = glm::dot(c_t, c_t);
        float B = glm::dot(c_t, S_t);
        float C = glm::dot(S_t, S_t) - 1;

        float discriminant = B * B - A * C;

        if (discriminant > 0)
        {
            /*
             * Case 1: no hits
             *
             * Case 2: t1 is negative, and t2 is positive
             * This happens when the image plane intersects the sphere, so it will be concave
             *
             * Case 3: t1 and t2 are both positive
             * t1 will be the closest one, and sphere will convex
             *
             * So if t_close is positive, we do not need to check whether t_far is also positive
             */

            float t_close = ((-B/A) - (sqrt(discriminant)/A));
            float t_far = ((-B/A) + (sqrt(discriminant)/A));

            /*
             * The rays stemming from the eye originate at their pixels, so all rays stem from their first valid intersection point
             */
            if(t_close > 0.0001)
            {
                // If t_close is positive, then both intersections are positive but t_close is closer
                if(intersection.no_intersect_flag || r.at(t_close).z > intersection.point.z )
                {
                    intersection.t = t_close;
                    intersection.sphere = inputImage.spheres[i];
                    intersection.point = r.at(intersection.t);

                    /* the sphere has been transformed to be size one at the origin, so the normal is the point of
                     * intersection - (0, 0, 0).
                     *
                     * To convert the normal back to world coordinates, apply the inverse transpose.
                    */
                    glm::vec3 normal = r_t.at(intersection.t);
                    glm::vec4 normal_homo = glm::transpose( inverse_scale * inverse_translate) * glm::vec4(normal, 0);
                    intersection.normal = glm::normalize(glm::vec3(normal_homo.x, normal_homo.y, normal_homo.z));
                    intersection.concave_flag = 0;
                    intersection.no_intersect_flag = 0;
                }
            }
            else if(t_far > 0.0001)
            {
                if(intersection.no_intersect_flag || r.at(t_far).z > intersection.point.z )
                {
                    // if t_close is negative and t_far is positive, then there is an intersection with the near plane
                    intersection.t = t_far;
                    intersection.sphere = inputImage.spheres[i];
                    intersection.point = r.at(intersection.t);

                    // the normal will be negative because the surface is concave
                    glm::vec3 normal = -r_t.at(intersection.t);
                    glm::vec4 normal_homo = glm::transpose(inverse_scale * inverse_translate) * glm::vec4(normal, 0);
                    intersection.normal = glm::normalize(glm::vec3(normal_homo.x, normal_homo.y, normal_homo.z));
                    intersection.concave_flag = 1;
                    intersection.no_intersect_flag = 0;
                }
            }
        }
    }
}

/*
 * function raytrace(r)
 *      if (ray.depth() > MAX_DEPTH) return black
 *      P = closest intersection of ray with all objects
 *      if( no intersection )
 *          return backgroundColor
 *      clocal = Sum(shadowRays(P,Lighti))
 *      cre = raytrace(rre)
 *      return (clocal+kre*cre+kra*cra)
 *   end
 */
glm::vec3 raytrace(Ray& r, const ImageInfo& inputImage, int depth)
{
    if( depth > MAX_DEPTH)
    {
        // return black
        return {0, 0, 0};
    }

    if(-0.001 < r.get_direction().x && r.get_direction().x < 0.001 && -0.001 < r.get_direction().y && r.get_direction().y < 0.001)
    {
        int x = 1;
    }

    Intersection intersection;
    closestIntersection(r, intersection, inputImage);

    if( intersection.no_intersect_flag)
    {
        if(depth == 1)
        {
            return {inputImage.back_r, inputImage.back_g, inputImage.back_b};
        }
        else
        {
            return {0, 0, 0};
        }
    }

    // ambient - Ka * Ia[c] * O[c]
    glm::vec3 c_local = glm::vec3(intersection.sphere.ka * inputImage.ambient_ir * intersection.sphere.r,
                                  intersection.sphere.ka * inputImage.ambient_ig * intersection.sphere.g,
                                  intersection.sphere.ka * inputImage.ambient_ib * intersection.sphere.b);

    glm::vec3 diffuse = glm::vec3(0, 0, 0);
    glm::vec3 specular = glm::vec3(0, 0, 0);
    glm::vec3 reflected = glm::vec3(0, 0, 0);

    for(int i = 0; i < inputImage.lights.size(); i++)
    {
        Light light = inputImage.lights[i];

        bool is_in_shadow = false;
        for(int s = 0; s < inputImage.spheres.size(); s++)
        {
            if(shadowRay(intersection, light, inputImage.spheres[s]))
            {
                is_in_shadow = true;
            }
        }

        bool is_light_behind = isLightBehind(intersection, light);

        if(!is_in_shadow && !is_light_behind)
        {
            glm::vec3 N = glm::normalize(intersection.normal);
            glm::vec3 L = glm::normalize(glm::vec3(light.posX, light.posY, light.posZ) - intersection.point);

            glm::vec3 V = glm::normalize(r.get_origin() - intersection.point);
            glm::vec3 R = glm::normalize(2 * glm::dot(N, L) * N - L);

            // diffuse : Kd * Ip[c] * (N dot L) *O[c]
            diffuse = glm::vec3(intersection.sphere.kd * light.ir * glm::dot(N, L) * intersection.sphere.r,
                                          intersection.sphere.kd * light.ig * glm::dot(N, L) * intersection.sphere.g,
                                          intersection.sphere.kd * light.ib * glm::dot(N, L) * intersection.sphere.b);

            // specular - Ks*Ip[c]*(R dot V)n
            specular = glm::vec3(
                    intersection.sphere.ks * light.ir * glm::pow(glm::dot(R, V), intersection.sphere.n),
                    intersection.sphere.ks * light.ig * glm::pow(glm::dot(R, V), intersection.sphere.n),
                    intersection.sphere.ks * light.ib * glm::pow(glm::dot(R, V), intersection.sphere.n));

            Ray reflection_ray = Ray(intersection.point, R);
            if(intersection.sphere.kr != 0)
            {
                reflected += intersection.sphere.kr * raytrace(reflection_ray, inputImage, depth + 1);
            }
            c_local += diffuse + specular + reflected;
        }
    }

    if(c_local.x > 1) c_local.x = 1;
    if(c_local.y > 1) c_local.y = 1;
    if(c_local.z > 1) c_local.z = 1;

    return c_local;
}

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input file>.txt" << std::endl;
        return 1;
    }

    // parse the input file and store the information in input_image
    ImageInfo input_image;
    getImageInfo(argv[1], input_image);

    // pixels[0] is the top left of the image and pixels[3*Width*Height-1] is the bottom right of the image.
    unsigned char *pixels;
    pixels = new unsigned char [3*input_image.width*input_image.height];

    // N distance from camer

    // TODO what if these are odd?
    float H = 1;
    float W = 1;
    float N = 1;

    int nCols = input_image.width;
    int nRows = input_image.height;
    glm::vec3 eye = glm::vec3(0, 0, 0);

    // Calculate the vectors across the horizontal and down the vertical viewport edges.
    glm::vec3 u = glm::vec3(1, 0, 0);
    glm::vec3 v = glm::vec3(0, 1, 0);
    glm::vec3 n = glm::vec3(0, 0, 1);

    glm::vec3 upper_left_vector = eye + glm::vec3(0, 0, -N) - float(W)*u + float(H)*v;

    int k = 0 ;
    for(int r = nRows - 1; r >= 0; r--)
    {
        for (int c = 0; c < input_image.width; c++)
        {
            float pixel_u_c = -W + W * (2 * float(c) / float(nCols));
            float pixel_v_r = -W + W * (2 * float(r) / float(nRows));
            // Point of pixel location in camera coordinates
            glm::vec3 P_pixel_camera = glm::vec3(pixel_u_c, pixel_v_r, -N);
            glm::vec3 P_pixel_world = eye - N*n + pixel_u_c*u + pixel_v_r*v;

            // TODO textbook also gives this in camera coordinates, what should I use
            glm::vec3 ray_direction = P_pixel_world - eye;

            Ray myRay = Ray(P_pixel_world, ray_direction);

            glm::vec3 color = raytrace(myRay, input_image, 1);

            pixels[k] = (unsigned char) (color.r * 255);
            pixels[k+1] = (unsigned char) (color.g * 255);
            pixels[k+2] = (unsigned char) (color.b * 255);
            k = k + 3 ;
        }
    }

    // TODO this is kinda a strange way to do this, trying to use c style strings
    std::string filename = "outputs/" + input_image.output;
    save_imageP3(input_image.width, input_image.height, &filename[0], pixels);
}
