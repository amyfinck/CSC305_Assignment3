//
// Created by Amy Finck on 2023-11-11.
//

#ifndef ASSIGNMENT3_RAYTRACE_H
#define ASSIGNMENT3_RAYTRACE_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstdlib>

#define MAX_DEPTH 5

void save_imageP6(int Width, int Height, char* fname,unsigned char* pixels);

struct Sphere
{
    std::string name;
    double posX, posY, posZ;
    double scaleX, scaleY, scaleZ;
    double r, g, b;
    double ka, kd, ks, kr, n;
};

// LIGHT <name> <pos x=""> <pos y=""> <pos z=""> <ir> <ig> <ib>
struct Light
{
    std::string name;
    double posX, posY, posZ;
    double ir, ig, ib;
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
    double back_r, back_g, back_b;
    double ambient_ir, ambient_ig, ambient_ib;
    std::string output;
};

#endif //ASSIGNMENT3_RAYTRACE_H
