#ifndef ASSIGNMENT3_RAY_H
#define ASSIGNMENT3_RAY_H

// TODO is this the best way to include this?
#include "external/glm/glm/glm.hpp"

class Ray
{
    public:
        Ray(const glm::vec3& o, const glm::vec3& dir, int dep)
        {
            origin = o;
            direction = dir;
            depth = dep;
        }

        glm::vec3 get_origin() const
        {
            return origin;
        }

        glm::vec3 get_direction() const
        {
            return direction;
        }

        int get_depth()
        {
            return depth;
        }

        void set_depth(int dep)
        {
            depth = dep;
        }

        glm::vec3 at(float t) const
        {
            return origin + t * direction;
        }

    private:
        glm::vec3 origin;
        glm::vec3 direction;
        int depth;
};

#endif //ASSIGNMENT3_RAY_H
