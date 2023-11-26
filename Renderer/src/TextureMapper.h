#include <vector>
#include <string>
#include "../libs/glm-0.9.7.2/glm/glm.hpp"
#include "../libs/sdw/Colour.h"
#include "UtilityFunctions.h"
#include "../libs/sdw/TextureMap.h"
#include "Triangle.h"
#pragma once

enum class TextureMappingType {
	NONE,
	FLAT,
	NORMAL
};

enum class TextureMappingPrimative {
	RECTANGLE
};

struct TextureVertex {
	size_t vertex;
	glm::vec2 textureCoord;
};

struct TextureMapper {
	TextureMappingPrimative primative = TextureMappingPrimative::RECTANGLE;
	TextureMappingType type = TextureMappingType::NONE;
	std::vector<std::vector<Colour>> textureMap = {};
	std::vector<std::vector<glm::vec3>> normalMap = {};
	std::vector<TextureVertex> vertexMapping = {};
	double scale = 1; // 1 => Stretch surface over the entire texture

	std::vector<std::vector<Colour>> TextureMapToColourArray(TextureMap& textureMap);
	
	// Transforms the normal map with the up direction as +ve in the z axis into a 2d vector of vec3s rotated to align the z axis with the normal
	// Rodrigues' rotation formula (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula):
	//
	// A = Normal Map Vector
	// V = Face Normal Vector
	// x = Angle between A and V
	// Z = Unit vector of world Z axis
	// N = Z X V/(||Z X V||)
	// A'= Acos(x) + (N X A)sin(x) + N(N.A)(1-cos(x)) = Transformed normal map vector
	std::vector<std::vector<glm::vec3>> NormalMapToVec3Array(const TextureMap& normalMap, const glm::vec3 normal);

	TextureMapper(TextureMap texture, TextureMap normal, glm::vec3 faceNormal, TextureMappingPrimative primative);
	TextureMapper(TextureMap texture, TextureMappingPrimative primative);
	TextureMapper();

	TextureVertex GetVertexMapping(size_t v) const;

	glm::vec2 BarycentricToMap(const Triangle& t, const glm::vec3 barycentric) const;

	Colour BarycentricTexture(const Triangle& t, const glm::vec3 barycentric) const;

	glm::vec3 BarycentricNormal(const Triangle& t, const glm::vec3 barycentric) const;

	void VertexMapping(const std::vector<glm::vec3>& vertices, std::vector<Triangle> triangles);
};