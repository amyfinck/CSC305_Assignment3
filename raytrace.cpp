#include <fstream>
#include <sstream>
#include "raytrace.h"
#include "ray.h"


#include "external/glm/glm/glm.hpp"
//#include "external/glm/glm/gtc/matrix_transform.hpp"

void save_imageP3(int Width, int Height, char* fname,unsigned char* pixels)
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname,Width,Height);
    fp = fopen(fname,"w");
    if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
    }
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    int k = 0 ;
    for(int j = 0; j < Height; j++) {

        for( int i = 0 ; i < Width; i++)
        {
            fprintf(fp," %d %d %d", pixels[k],pixels[k+1],pixels[k+2]) ;
            k = k + 3 ;
        }
        fprintf(fp,"\n") ;
    }
    fclose(fp);
}

/*
 * file_name is argv[1] from main, we want to pass by reference and promise not to modify it, hence the const and &
 * input_image is the struct we want to fill with the information from the input file
 */
void getImageInfo(const std::string& file_name, ImageInfo& input_image)
{
    std::ifstream inputFile(file_name);

    if (!inputFile.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        exit(1);
    }

    std::string keyword;
    int value;
    while(inputFile >> keyword)
    {
        if(keyword == "NEAR")
        {
            inputFile >> input_image.near;
        }
        if(keyword == "LEFT")
        {
            inputFile >> input_image.left;
        }
        if(keyword == "RIGHT")
        {
            inputFile >> input_image.right;
        }
        if(keyword == "BOTTOM")
        {
            inputFile >> input_image.bottom;
        }
        if(keyword == "TOP")
        {
            inputFile >> input_image.top;
        }
        if(keyword == "RES")
        {
            inputFile >> input_image.width;
            inputFile >> input_image.height;
        }
        if(keyword == "SPHERE")
        {
            Sphere sphere;
            inputFile >> sphere.name >> sphere.posX >> sphere.posY >> sphere.posZ
                      >> sphere.scaleX >> sphere.scaleY >> sphere.scaleZ
                      >> sphere.r >> sphere.g >> sphere.b
                      >> sphere.ka >> sphere.kd >> sphere.ks >> sphere.kr >> sphere.n;
            sphere.initDefaults();
            input_image.spheres.push_back(sphere);
        }
        if(keyword == "LIGHT")
        {
            Light light;
            inputFile >> light.name >> light.posX >> light.posY >> light.posZ
                      >> light.ir >> light.ig >> light.ib;
            input_image.lights.push_back(light);
        }
        if(keyword == "BACK")
        {
            inputFile >> input_image.back_r >> input_image.back_g >> input_image.back_b;
        }
        if(keyword == "AMBIENT")
        {
            inputFile >> input_image.ambient_ir >> input_image.ambient_ig >> input_image.ambient_ib;
        }
        if(keyword == "OUTPUT")
        {
            inputFile >> input_image.output;
        }
    }
    inputFile.close();
}

glm::mat4 inverseTranslate(float x, float y, float z)
{
    glm::mat4 inverseTranslateMatrix = glm::mat4(0);
    inverseTranslateMatrix[0][0] = 1;
    inverseTranslateMatrix[1][1] = 1;
    inverseTranslateMatrix[2][2] = 1;
    inverseTranslateMatrix[3][3] = 1;
    inverseTranslateMatrix[3][0] = -x;
    inverseTranslateMatrix[3][1] = -y;
    inverseTranslateMatrix[3][2] = -z;

    return inverseTranslateMatrix;
}

glm::mat4 inverseScale(float x, float y, float z)
{
    glm::mat4 inverseScaleMatrix = glm::mat4(0);
    inverseScaleMatrix[0][0] = 1 / x;
    inverseScaleMatrix[1][1] = 1 / y;
    inverseScaleMatrix[2][2] = 1 / z;
    inverseScaleMatrix[3][3] = 1;

    return inverseScaleMatrix;
}

// returns the t value of the closest intersection, or -1 if there is no intersection
float closestIntersection(Ray& r, const ImageInfo& inputImage)
{
    glm::vec3 S = r.get_origin();
    glm::vec3 c = r.get_direction();

    float t = -1;

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

        float A = glm::dot(c_t, c_t);
        float B = glm::dot(c_t, S_t);
        float C = glm::dot(S_t, S_t) - 1;

        float discriminant = B * B - A * C;

        if (discriminant > 0)
        {
            float t1 = (-B + sqrt(discriminant)) / (2 * A);
            float t2 = (-B - sqrt(discriminant)) / (2 * A);
            t = std::min(t1, t2);
        }
    }
    return t;
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
glm::vec3 raytrace(Ray& r, const ImageInfo& inputImage)
{
    if( r.get_depth() > MAX_DEPTH)
    {
        // return black
        return glm::vec3(0, 0, 0);
    }

    float closest_intersection_t = closestIntersection(r, inputImage);
    if( closest_intersection_t < 0)
    {
        // no intersection
        return glm::vec3(inputImage.back_r, inputImage.back_g, inputImage.back_b);
    }

    // Intersects at this point
    glm::vec3 P = r.at(closest_intersection_t);

    // TODO here we are assuming sphere 1 because we are hardcoding.
    // ambient - Ka * Ia[c] * O[c]
    glm::vec3 c_local = glm::vec3(inputImage.spheres[0].ka * inputImage.ambient_ir * inputImage.spheres[0].r,
                                  inputImage.spheres[0].ka * inputImage.ambient_ig * inputImage.spheres[0].g,
                                  inputImage.spheres[0].ka * inputImage.ambient_ib * inputImage.spheres[0].b);

    for(int i = 0; i < inputImage.lights.size(); i++)
    {
        Light light = inputImage.lights[i];

        // TODO think about these
        glm::vec3 N = glm::normalize(P);
        glm::vec3 L = glm::normalize(glm::vec3(light.posX, light.posY, light.posZ) - P);
        glm::vec3 V = glm::normalize(r.get_origin() - P);
        glm::vec3 R = glm::normalize(2 * glm::dot(N, L) * N - L);

        // diffuse - Kd * Ip[c] * (N dot L) *O[c]
        glm::vec3 diffuse = glm::vec3(inputImage.spheres[0].kd * light.ir * glm::dot(N, L) * inputImage.spheres[0].r,
                            inputImage.spheres[0].kd * light.ig * inputImage.spheres[0].g * glm::dot(N, L),
                            inputImage.spheres[0].kd * light.ib * inputImage.spheres[0].b * glm::dot(N, L));

        // specular - Ks*Ip[c]*(R dot V)n
        glm::vec3 specular = glm::vec3(inputImage.spheres[0].ks * light.ir * glm::pow(glm::dot(R, V), inputImage.spheres[0].n),
                             inputImage.spheres[0].ks * light.ig * inputImage.spheres[0].g * glm::pow(glm::dot(R, V), inputImage.spheres[0].n),
                             inputImage.spheres[0].ks * light.ib * inputImage.spheres[0].b * glm::pow(glm::dot(R, V), inputImage.spheres[0].n));

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

            glm::vec3 color = raytrace(myRay, input_image);

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
