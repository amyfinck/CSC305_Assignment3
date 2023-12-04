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
