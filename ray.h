#ifndef ASSIGNMENT3_RAY_H
#define ASSIGNMENT3_RAY_H

// TODO is this the best way to include this?
#include "external/glm/glm/glm.hpp"

class ray
{
public:
    // empty construcor
    ray() {}

    // origin is staring point, direction is vector of location
    ray(const glm::vec3& origin, const glm::vec3& direction) : orig(origin), dir(direction) {}

    glm::vec3 origin() const  { return orig; }
    glm::vec3 direction() const { return dir; }

    // inputs traversal length and outputs new ray position
    glm::vec3 at(float t) const
    {
        return orig + t * dir;
    }

private:
    glm::vec3 orig;
    glm::vec3 dir;
};


#endif //ASSIGNMENT3_RAY_H
