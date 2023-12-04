#ifndef ASSIGNMENT3_RAYTRACE_H
#define ASSIGNMENT3_RAYTRACE_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <tuple>
#include "external/glm/glm/glm.hpp"
#include "external/glm/glm/gtc/matrix_transform.hpp"
#include "external/glm/glm/gtc/matrix_inverse.hpp"
#include "external/glm/glm/mat4x4.hpp"

#define MAX_DEPTH 3

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

struct Intersection
{
    float t;
    glm::vec3 normal;
    glm::vec3 point;
    Sphere sphere;
    int no_intersect_flag;
    int concave_flag;
};

#endif //ASSIGNMENT3_RAYTRACE_H
