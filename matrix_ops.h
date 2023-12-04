#include "raytrace.h"

glm::mat4 inverseTranslate(float x, float y, float z)
{
    glm::mat4 Identity = glm::mat4( 1. );
    glm::mat4 inverseTranslateMatrix = glm::translate( Identity, glm::vec3( -x,-y,-z ) );
    return inverseTranslateMatrix;
}

glm::mat4 inverseScale(float x, float y, float z)
{
    glm::mat4 Identity = glm::mat4( 1. );
    glm::mat4 inverseScaleMatrix = glm::scale( Identity, glm::vec3( 1/x,1/y,1/z ) );
    return inverseScaleMatrix;
}