//
// Created by Amy Finck on 2023-11-11.
//

#ifndef ASSIGNMENT3_RAYTRACE_H
#define ASSIGNMENT3_RAYTRACE_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "external/glm/glm/glm.hpp"

#define MAX_DEPTH 5

void save_imageP6(int Width, int Height, char* fname,unsigned char* pixels);

struct Sphere
{
    std::string name;
    float posX, posY, posZ;
    float scaleX, scaleY, scaleZ;
    float r, g, b;
    float ka, kd, ks, kr, n;
    glm::mat4 transform;
    glm::mat4 inverseTransform;

    void initDefaults() {
        // This function initializes the transform array
        transformArray();
        inverseTransformArray();
    }

    private:
        // Function to set the transform array based on member variables
        void transformArray()
        {
            transform[0][0] = (float)scaleX;
            transform[0][1] = 0;
            transform[0][2] = 0;
            transform[0][3] = 0;

            transform[1][1] = (float)scaleY;
            transform[1][0] = 0;
            transform[1][2] = 0;
            transform[1][3] = 0;

            transform[2][2] = (float)scaleZ;
            transform[2][0] = 0;
            transform[2][1] = 0;
            transform[2][3] = 0;

            transform[3][0] = (float)posX;
            transform[3][1] = (float)posY;
            transform[3][2] = (float)posZ;
            transform[3][3] = 1;
        }

        void inverseTransformArray()
        {
            inverseTransform[0][0] = 1 / (float)scaleX;
            inverseTransform[0][1] = 0;
            inverseTransform[0][2] = 0;
            inverseTransform[0][3] = (float) -posX;

            inverseTransform[1][1] = 1 / (float)scaleY;
            inverseTransform[1][0] = 0;
            inverseTransform[1][2] = 0;
            inverseTransform[1][3] = (float) -posY;

            inverseTransform[2][2] = 1 / (float)scaleZ;
            inverseTransform[2][0] = 0;
            inverseTransform[2][1] = 0;
            inverseTransform[2][3] = (float) -posZ;

            inverseTransform[3][0] = 0;
            inverseTransform[3][1] = 0;
            inverseTransform[3][2] = 0;
            inverseTransform[3][3] = 1;
        }
};

// LIGHT <name> <pos x=""> <pos y=""> <pos z=""> <ir> <ig> <ib>
struct Light
{
    std::string name;
    float posX, posY, posZ;
    float ir, ig, ib;
};

struct ImageInfo
{
    std::string fname;
    int near;
    int left;
    int right;
    int bottom;
    int top;
    int width;
    int height;
    std::vector<Sphere> spheres;
    std::vector<Light> lights;
    float back_r, back_g, back_b;
    float ambient_ir, ambient_ig, ambient_ib;
    std::string output;
};

#endif //ASSIGNMENT3_RAYTRACE_H
