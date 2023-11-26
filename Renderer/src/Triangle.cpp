#include "Triangle.h"
#include "UtilityFunctions.h"

Triangle::Triangle(size_t v0, size_t v1, size_t v2, glm::vec3 normal) {
    this->vertices.push_back(v0);
    this->vertices.push_back(v1);
    this->vertices.push_back(v2);
    this->normal = normal;
}

std::pair<glm::vec3, glm::vec3> Triangle::BoundingBox(std::vector<glm::vec3>& coords) {
    glm::vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
    glm::vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (size_t v : this->vertices) {
        for (size_t axis = 0; axis < 3; axis++) {
            vec3SetAt(min, axis, fmin(vec3At(min, axis), vec3At(coords[v], axis)));
            vec3SetAt(max, axis, fmax(vec3At(max, axis), vec3At(coords[v], axis)));
        }
    }

    this->bbox = std::make_pair(min, max);
    return this->bbox;
}

glm::vec3 Triangle::Centre(const std::vector<glm::vec3>& coords) {
    return (coords[vertices[0]] + coords[vertices[1]] + coords[vertices[2]]) / 3.0f;
}