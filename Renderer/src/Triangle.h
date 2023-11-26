#include <vector>
#include <string>
#include "../libs/glm-0.9.7.2/glm/glm.hpp"
#pragma once

struct Triangle {
	std::vector<size_t> vertices;
	glm::vec3 normal;
	std::pair<glm::vec3, glm::vec3> bbox;

    Triangle(size_t v0, size_t v1, size_t v2, glm::vec3 normal);
	std::pair<glm::vec3, glm::vec3> BoundingBox(std::vector<glm::vec3>& coords);
    glm::vec3 Centre(const std::vector<glm::vec3>& coords);
};