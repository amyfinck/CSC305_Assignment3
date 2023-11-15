#include <fstream>
#include <sstream>
#include "raytrace.h"

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

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input file>.txt" << std::endl;
        return 1;
    }

    ImageInfo input_image;

    getImageInfo(argv[1], input_image);

    std::cout << "Back: " << input_image.back_r << " " << input_image.back_g << " " << input_image.back_b << std::endl;
    std::cout << "Ambient: " << input_image.ambient_ir << " " << input_image.ambient_ig << " " << input_image.ambient_ib << std::endl;
    std::cout << "Output: " << input_image.output << std::endl;

    //int width = input_image.width;	// Move these to your setup function. The actual resolution will be
    char fname3[20] = "outputs/sceneP3.ppm"; //This should be set based on the input file
    unsigned char *pixels;

    // pixels[0] is the top left of the image and
    // pixels[3*Width*Height-1] is the bottom right of the image.
    pixels = new unsigned char [3*input_image.width*input_image.height];

    // This loop just creates a gradient for illustration purposes only. You will not use it.
    float scale = 128.0 / (float) input_image.width ;
    int k = 0 ;
    for(int i = 0; i < input_image.height; i++) {
        std::clog << "\rScanlines remaining: " << (input_image.width - i) << ' ' << std::flush;
        for (int j = 0; j < input_image.width; j++) {
            int c = (i+j)*scale ;
            pixels[k] = c;
            pixels[k+1] = c;
            pixels[k+2] = c;
            k = k + 3 ;
        }
    }
    // TODO this is kinda a strange way to do this, trying to use c style strings
    std::string filename = "outputs/" + input_image.output;
    save_imageP3(input_image.width, input_image.height, &filename[0], pixels);
}
