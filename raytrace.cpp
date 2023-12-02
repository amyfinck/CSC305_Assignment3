#include <fstream>
#include <sstream>
#include <tuple>
#include "raytrace.h"
#include "ray.h"
#include "matrix_ops.h"
#include "file_processing.h"

// returns the t value of the closest intersection, or -1 if there is no intersection
void closestIntersection(Ray& r, Intersection& intersection, const ImageInfo &inputImage)
{
    glm::vec3 S = r.get_origin();
    glm::vec3 c = r.get_direction();

    for(int i = 0; i < inputImage.spheres.size(); i++)
    {
        Sphere sphere = inputImage.spheres[i];
        glm::vec4 S_homo = glm::vec4(r.get_origin(), 1);
        glm::vec4 c_homo = glm::vec4(r.get_direction(), 0);

        glm::mat4 inverse_translate = inverseTranslate(sphere.posX, sphere.posY, sphere.posZ);
        glm::mat4 inverse_scale = inverseScale(sphere.scaleX, sphere.scaleY, sphere.scaleZ);
        glm::vec4 S_t_homo = inverse_scale * inverse_translate * S_homo;
        glm::vec4 c_t_homo = inverse_scale * inverse_translate * c_homo;
        glm::vec3 S_t = glm::vec3(S_t_homo.x, S_t_homo.y, S_t_homo.z);
        glm::vec3 c_t = glm::vec3(c_t_homo.x, c_t_homo.y, c_t_homo.z);
        Ray r_t = Ray(S_t, c_t, 1);

        float A = glm::dot(c_t, c_t);
        float B = glm::dot(c_t, S_t);
        float C = glm::dot(S_t, S_t) - 1;

        float discriminant = B * B - A * C;

        if (discriminant > 0)
        {
            float t1 = ((-B/A) + (sqrt(discriminant)/A));
            float t2 = ((-B/A) - (sqrt(discriminant)/A));

            float t_candidate = std::min(t1, t2);
            if(t_candidate < 1) t_candidate = std::max(t1, t2);

            if(t_candidate >= 1)
            {
                intersection.sphere = inputImage.spheres[i];
                intersection.t = t_candidate;

                intersection.no_intersect_flag = 0;

                // TODO - this gives me the correct answer but I dont understand why I don't need to transform it first
                intersection.point = r.at(intersection.t);

                glm::vec3 normal = r_t.at(intersection.t);
                glm::vec4 normal_homo = glm::transpose( inverse_scale * inverse_translate) * glm::vec4(normal, 0);
                intersection.normal = glm::normalize(glm::vec3(normal_homo.x, normal_homo.y, normal_homo.z));

                return;
            }
        }
    }
    intersection.no_intersect_flag = 1;
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

    Intersection intersection;
    closestIntersection(r, intersection, inputImage);

    if( intersection.no_intersect_flag)
    {
        return {inputImage.back_r, inputImage.back_g, inputImage.back_b};
    }

    // ambient - Ka * Ia[c] * O[c]
    glm::vec3 c_local = glm::vec3(intersection.sphere.ka * inputImage.ambient_ir * intersection.sphere.r,
                                  intersection.sphere.ka * inputImage.ambient_ig * intersection.sphere.g,
                                  intersection.sphere.ka * inputImage.ambient_ib * intersection.sphere.b);

    for(int i = 0; i < inputImage.lights.size(); i++)
    {
        Light light = inputImage.lights[i];

        glm::vec3 N = glm::normalize(intersection.normal);
        glm::vec3 L = glm::normalize(glm::vec3(light.posX, light.posY, light.posZ) - intersection.point);
        glm::vec3 V = glm::normalize(intersection.point - r.get_origin());
        glm::vec3 R = glm::normalize(2 * glm::dot(N, L) * N - L);

        // diffuse : Kd * Ip[c] * (N dot L) *O[c]
        glm::vec3 diffuse = glm::vec3(intersection.sphere.kd * light.ir * glm::dot(N, L) * intersection.sphere.r,
                                      intersection.sphere.kd * light.ig * glm::dot(N, L) * intersection.sphere.g ,
                                      intersection.sphere.kd * light.ib * glm::dot(N, L) * intersection.sphere.b );

        // specular - Ks*Ip[c]*(R dot V)n
        glm::vec3 specular = glm::vec3(intersection.sphere.ks * light.ir * glm::pow(glm::dot(R, V), intersection.sphere.n),
                                       intersection.sphere.ks * light.ig * glm::pow(glm::dot(R, V), intersection.sphere.n),
                                       intersection.sphere.ks * light.ib * glm::pow(glm::dot(R, V), intersection.sphere.n));

        c_local += diffuse + specular;
    }

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

            Ray myRay = Ray(eye, ray_direction, 1);

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
