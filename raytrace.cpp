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

glm::vec3 ray_color(const ray& r)
{
    glm::vec3 unit_direction = r.direction() / glm::length(r.direction());
    float a = 0.5f *(unit_direction.y + 1.0);
    glm::vec3 retval =  (1.0f-a)*glm::vec3(1.0, 1.0, 1.0) + a*glm::vec3(0.5, 0.7, 1.0);
    return retval;
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

    // N distance from camera
    float N = 1.0f;

    Sphere testSphere = input_image.spheres[0];

    // TODO what if these are odd?
    float H = float(input_image.width) / 2;
    float W = float(input_image.height) / 2;
    int nCols = input_image.width;
    int nRows = input_image.height;
    glm::vec3 eye = glm::vec3(0, 0, 0);

    // Calculate the vectors across the horizontal and down the vertical viewport edges.
    glm::vec3 u = glm::vec3(1, 0, 0);
    glm::vec3 v = glm::vec3(0, 1, 0);
    glm::vec3 n = glm::vec3(0, 0, -1);

    // Calculate the location of the upper left pixel.
    //int pixel00_u = - W + W * (2 * 0/ nCols);
    //int pixel00_v = - W + W * (2 * (nRows - 1) / nRows);

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

            // TODO textbook also gives this in camera coordinates, what should I use???
            glm::vec3 ray_direction = P_pixel_camera - eye;

            if(k % 10 == 0) {
                std::cout << "pixel_u_c: " << pixel_u_c << std::endl;
                std::cout << "pixel_v_r: " << pixel_v_r << std::endl;
                std::cout << "ray_direction: " << ray_direction.x << ", " << ray_direction.y << ", " << ray_direction.z << std::endl;
            }

            //ray myRay = ray(camera_center, ray_direction);

            //glm::vec3 color = ray_color(myRay);

            pixels[k] = (unsigned char) ((pixel_u_c + W) / (2 * W) * 255);
            pixels[k+1] = (unsigned char) ((pixel_v_r + H) / (2 * H) * 255);
            pixels[k+2] = (unsigned char) 0;
            k = k + 3 ;
        }
    }

    // TODO this is kinda a strange way to do this, trying to use c style strings
    std::string filename = "outputs/" + input_image.output;
    save_imageP3(input_image.width, input_image.height, &filename[0], pixels);
}
