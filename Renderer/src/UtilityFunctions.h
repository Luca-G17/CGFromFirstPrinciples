#include <vector>
#include <string>
#include "../libs/glm-0.9.7.2/glm/glm.hpp"
#include "../libs/sdw/Colour.h"

#pragma once

std::vector<std::string> Split(std::string str, std::string delimiter);

Colour Uint32ToColour(uint32_t i);

glm::vec3 Uint32ToVec3(uint32_t i);

Colour ScaleColour(Colour c, float s);

// Component wise vec3 scaling
glm::vec3 ScaleVec3(glm::vec3 v, glm::vec3 scale);

// 1 / vector
glm::vec3 InvertVec3(glm::vec3 v);

Colour AddColours(Colour c0, Colour c1);

double vec3At(const glm::vec3& v, size_t axis);

void vec3SetAt(glm::vec3& v, size_t axis, double d);

float CyclicalValue(const double framesPerCycle, const double step);

glm::mat3 RotationMatrix(const float x, const float y, const float z);
