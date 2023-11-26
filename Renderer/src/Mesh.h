#include <vector>
#include <string>
#include "../libs/glm-0.9.7.2/glm/glm.hpp"
#include "Triangle.h"
#include "UtilityFunctions.h"
#include "TextureMapper.h"
#pragma once

enum class Shading {
	PROXIMITY,
	ANGLE_OF_INCIDENCE,
	SPECULAR,
	PROXY_AOI_AND_SPECULAR,
	GOURAUD,
	PHONG,
	NONE
};

struct MeshPropertyUpdate {
	double* value;
	double start;
	double finish;
	double speed;
	double step;
	bool cyclical;

	MeshPropertyUpdate(double* value, double start, double finish, double speed, bool cyclical);

	MeshPropertyUpdate(double* value, double finish, double speed, bool cyclical);

	bool Update();

	bool NotFinished();
};

struct MeshVec3PropertyUpdate {
	glm::vec3* value;
	glm::vec3 start;
	glm::vec3 finish;
	glm::vec3 bezierMid;
	float speed;
	float step;
	bool cyclical;
	bool bezier = false;

	MeshVec3PropertyUpdate(glm::vec3* value, glm::vec3 start, glm::vec3 finish, double speed, bool cyclical);

	MeshVec3PropertyUpdate(glm::vec3* value, glm::vec3 start, glm::vec3 bezierMid ,glm::vec3 finish, double speed, bool cyclical);

	bool Update();

	bool NotFinished(double v);
};

struct Mesh {
	std::vector<MeshPropertyUpdate> propertyUpdates = {};
	std::vector<MeshVec3PropertyUpdate> vec3PropertyUpdates= {};
	std::vector<Triangle> triangles;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> vertexNormals;
	std::pair<glm::vec3, glm::vec3> bbox;
	TextureMapper texture = {};

	Colour colour;
	double smoothness;
	double refractiveIndex;
	double transparancy;
	std::string name;
	glm::vec3 currentScale = {1, 1, 1};
	glm::vec3 scaleTransition = {1, 1, 1};
	glm::vec3 translation;
	glm::vec3 position;
	Shading shadingType = Shading::PROXY_AOI_AND_SPECULAR;

	void SetTextureMapper(TextureMap textureMap, TextureMap normalMap);

	void SetTextureMapper(TextureMap map);

	const glm::vec3& GetVertex(const size_t t, const size_t v) const;

	void SetVertex(const size_t t, const size_t v, const glm::vec3 vertex);

	// Vertex and face normals will require recomputation if I add rotations
	const glm::vec3& GetVertexNormal(size_t t, size_t v) const;

	const Colour GetColour() const;

	void UpdateMesh();

	std::pair<glm::vec3, glm::vec3> MeshBoundingBox();

	glm::vec3 Centre();

	void Translate(glm::vec3 translation);

	void Scale(glm::vec3 scale);

	void Rotate(glm::vec3 rotation);

	void ComputeBoundingBoxes();
};