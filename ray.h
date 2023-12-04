#ifndef ASSIGNMENT3_RAY_H
#define ASSIGNMENT3_RAY_H

#include "external/glm/glm/glm.hpp"

class Ray
{
    public:
        Ray(const glm::vec3& o, const glm::vec3& dir)
        {
            origin = o;
            direction = dir;
        }

        glm::vec3 get_origin() const
        {
            return origin;
        }

        glm::vec3 get_direction() const
        {
            return direction;
        }

        glm::vec3 at(float t) const
        {
            return origin + t * direction;
        }

    private:
        glm::vec3 origin{};
        glm::vec3 direction{};
};

#endif //ASSIGNMENT3_RAY_H
