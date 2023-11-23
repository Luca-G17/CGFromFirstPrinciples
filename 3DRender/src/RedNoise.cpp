#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include <cmath>
#include <Colour.h>
#include <TextureMap.h>
#include <iostream>
#include <sstream>
#include <map>
#include <ModelTriangle.h>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#include <filesystem>
#include <sys/stat.h>
#include <chrono>

#define WIDTH 440
#define HEIGHT 300
#define FPS_CAP 60

struct Shape2D {
	std::vector<CanvasPoint> points;
	Colour colour;
};

struct Shape3D {
	std::vector<glm::vec3> points;
	Colour colour;
};

struct DepthPoint {
	Colour colour;
	float depth;
};

Colour ScaleColour(Colour c, float s) {
	return Colour(
		fmax(fmin(c.red * s, 255), 0),
		fmax(fmin(c.green * s, 255), 0),
		fmax(fmin(c.blue * s, 255), 0)
		);
}

glm::vec3 ScaleVec3(glm::vec3 v, glm::vec3 scale) {
	return glm::vec3(
		v.x * scale.x,
		v.y * scale.y,
		v.z * scale.z
	);
}

glm::vec3 InvertVec3(glm::vec3 v) {
	return glm::vec3(
		1.0f / v.x,
		1.0f / v.y,
		1.0f / v.z
	);
}

Colour AddColours(Colour c0, Colour c1) {
	return Colour(
		fmax(fmin(c0.red + c1.red, 255), 0),
		fmax(fmin(c0.green + c1.green, 255), 0),
		fmax(fmin(c0.blue + c1.blue, 255), 0)
	);
}

std::vector<std::string> Split(std::string str, std::string delimiter) {
	size_t pos = 0;
	std::string token;
	std::vector<std::string> segments;
	while ((pos = str.find(delimiter)) != std::string::npos) { // while string contains delimiters
		token = str.substr(0, pos);
		segments.push_back(token);
		str.erase(0, pos + delimiter.length());
	}
	segments.push_back(str);
	return segments;
}

void RemoveChar(std::string str, char c) {
	str.erase(std::remove(str.begin(), str.end(), c), str.end());
	// std::remove moves all elements that are not c to the front of the range, it returns an iterator to the first of the unmoved elements
	// std::erase(iterator first, iterator last) erases the portion of the string inside the iterator range 
}

double vec3At(glm::vec3& v, size_t axis) {
	if (axis == 0) return v.x;
	if (axis == 1) return v.y;
	return v.z;
}

void vec3SetAt(glm::vec3& v, size_t axis, double d) {
	switch (axis)
	{
		case 0:
			v.x = d;
			break;
		case 1:
			v.y = d;
			break;
		case 2:
			v.z = d;
			break;
		default:
			break;
	}
}

struct Triangle {
	std::vector<size_t> vertices;
	glm::vec3 normal;
	std::pair<glm::vec3, glm::vec3> bbox;

	Triangle(size_t v0, size_t v1, size_t v2, glm::vec3 normal) {
		this->vertices.push_back(v0);
		this->vertices.push_back(v1);
		this->vertices.push_back(v2);
		this->normal = normal;
	}

	std::pair<glm::vec3, glm::vec3> BoundingBox(std::vector<glm::vec3>& coords) {
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

	glm::vec3 Centre(const std::vector<glm::vec3>& coords) {
		return (coords[vertices[0]] + coords[vertices[1]] + coords[vertices[2]]) / 3.0f;
	}
};

float CyclicalValue(const double framesPerCycle, const double step) {
	return -(1.0f / 2.0) * cos(2.0 * M_PI * step / framesPerCycle) + 1.0 / 2.0;
}

// Container for a stepwise update of a property on a mesh
// Assuming 24 Frame/s, a 10 second transition requires 240 Frames
struct MeshPropertyUpdate {
	double* value;
	double start;
	double finish;
	double speed;
	double step;
	bool cyclical;

	MeshPropertyUpdate(double* value, double start, double finish, double speed, bool cyclical) {
		this->value = value;
		this->start = start;
		this->finish = finish;
		this->speed = speed;
		this->step = 0;
		this->cyclical = cyclical;
	}

	MeshPropertyUpdate(double* value, double finish, double speed, bool cyclical) {
		this->value = value;
		this->start = *value;
		this->finish = finish;
		this->speed = speed;
		this->step = 0;
		this->cyclical = cyclical;
	}

	bool Update() {
		step++;
		double v = CyclicalValue(speed, step);
		*value = start + (finish - start) * v; 
		return NotFinished();
	}

	bool NotFinished() {
		if (cyclical)
			return true;
		
		return ((finish - *value) * (finish - start) > 0);
	}
};

struct MeshVec3PropertyUpdate {
	glm::vec3* value;
	glm::vec3 start;
	glm::vec3 finish;
	float speed;
	float step;
	bool cyclical;
	
	MeshVec3PropertyUpdate(glm::vec3* value, glm::vec3 start, glm::vec3 finish, double speed, bool cyclical) {
		this->value = value;
		this->start = start;
		this->finish = finish;
		this->speed = speed;
		this->step = 0;
		this->cyclical = cyclical;
	}

	bool Update() {
		step++;
		float v = CyclicalValue(speed, step);
		*value = start + (finish - start) * v;
		return NotFinished();
	}

	bool NotFinished() {
		if (cyclical)
			return true;
		
		return !(glm::dot(finish - start, finish - *value) <= 0); 
	}
};

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

// TODO: this is potentially stupid
std::vector<std::pair<size_t, glm::vec3>> GetRectangleCorners(const std::vector<glm::vec3>& vertices) {
	std::vector<std::pair<size_t, glm::vec3>> corners;
	for (size_t i = 0; i < vertices.size(); i++) {
		corners.push_back(std::make_pair(i, vertices[i]));
	}
	return corners;
}

int DetermineQuadrant(const glm::vec3& normal) {
	return (normal.z > 0) ? 1 : -1;
}

Colour Uint32ToColour(uint32_t i) {
	int r = (i & 0x00FF0000) >> 16;
	int g = (i & 0x0000FF00) >> 8;
	int b = i & 0x000000FF;
	return Colour(r, g, b);
}

glm::vec3 Uint32ToVec3(uint32_t i) {
	int x = (i & 0x00FF0000) >> 16;
	int y = (i & 0x0000FF00) >> 8;
	int z = i & 0x000000FF;
	return glm::vec3(x, y, z);
}

struct TextureMapper {
	TextureMappingPrimative primative = TextureMappingPrimative::RECTANGLE;
	TextureMappingType type = TextureMappingType::NONE;
	std::vector<std::vector<Colour>> textureMap = {};
	std::vector<std::vector<glm::vec3>> normalMap = {};
	std::vector<TextureVertex> vertexMapping = {};
	double scale = 1; // 1 => Stretch surface over the entire texture

	std::vector<std::vector<Colour>> TextureMapToColourArray(TextureMap& textureMap) {
		std::vector<std::vector<Colour>> colourArray;
		for (size_t y = 0; y < textureMap.height; y++) {
			colourArray.push_back({});
			for (size_t x = 0; x < textureMap.width; x++) {
				colourArray[y].push_back(Uint32ToColour(textureMap.pixels[y * textureMap.width + x]));
			}
		}
		return colourArray;
	}
	
	// Transforms the normal map with the up direction as +ve in the z axis into a 2d vector of vec3s rotated to align the z axis with the normal
	// Rodrigues' rotation formula (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula):
	//
	// A = Normal Map Vector
	// V = Face Normal Vector
	// x = Angle between A and V
	// Z = Unit vector of world Z axis
	// N = Z X V/(||Z X V||)
	// A'= Acos(x) + (N X A)sin(x) + N(N.A)(1-cos(x)) = Transformed normal map vector
	std::vector<std::vector<glm::vec3>> NormalMapToVec3Array(const TextureMap& normalMap, const glm::vec3 normal) {
		std::vector<std::vector<glm::vec3>> vec3Array;
		const glm::vec3 Z(0, 0, 1);
		for (size_t y = 0; y < normalMap.height; y++) {
			vec3Array.push_back({});
			for (size_t x = 0; x < normalMap.width; x++) {
				glm::vec3 A = (Uint32ToVec3(normalMap.pixels[y * normalMap.width + x]) * (1.0f / 127.5f) + (-1.0f * glm::vec3(1, 1, 1)));
				const float cosTheta = glm::dot(Z, normal);
				const float sinTheta = sqrt(1 - (cosTheta * cosTheta));
				const glm::vec3 N = glm::normalize(glm::cross(Z, normal));
				const glm::vec3 vTransformed = A * cosTheta + glm::cross(N, A) * sinTheta + N * glm::dot(N, A) * (1 - cosTheta);
				vec3Array[y].push_back(vTransformed);
			}
		}
		return vec3Array;
	}

	TextureMapper(TextureMap texture, TextureMap normal, glm::vec3 faceNormal, TextureMappingPrimative primative) {
		this->primative = primative;
		this->textureMap = TextureMapToColourArray(texture);
		this->normalMap = NormalMapToVec3Array(normal, faceNormal);
		this->type = TextureMappingType::NORMAL;
	}

	TextureMapper(TextureMap texture, TextureMappingPrimative primative) {
		this->primative = primative;
		this->textureMap = TextureMapToColourArray(texture);
		this->type = TextureMappingType::FLAT;
	}

	TextureMapper() {
	}

	TextureVertex GetVertexMapping(size_t v) const {
		for (const TextureVertex& t : vertexMapping) {
			if (t.vertex == v)
				return t;
		}
		return {};
	}

	glm::vec2 BarycentricToMap(const Triangle& t, const glm::vec3 barycentric) const {
		return GetVertexMapping(t.vertices[0]).textureCoord * barycentric[0] +  
			   GetVertexMapping(t.vertices[1]).textureCoord * barycentric[1] +
			   GetVertexMapping(t.vertices[2]).textureCoord * barycentric[2];
	}

	Colour BarycentricTexture(const Triangle& t, const glm::vec3 barycentric) const {
		glm::vec2 c = BarycentricToMap(t, barycentric);
		return textureMap[(int)c.y][(int)c.x];
	}

	glm::vec3 BarycentricNormal(const Triangle& t, const glm::vec3 barycentric) const {
		glm::vec2 c = BarycentricToMap(t, barycentric);
		return normalMap[(int)c.y][(int)c.x];
	}

	void VertexMapping(const std::vector<glm::vec3>& vertices, std::vector<Triangle> triangles) {
		switch (primative) {
			case TextureMappingPrimative::RECTANGLE: {
				Triangle t1 = triangles[0];
				Triangle t2 = triangles[1];
				glm::vec3 centre = (t1.Centre(vertices) + t2.Centre(vertices)) / 2.0f;
				std::vector<std::pair<size_t, glm::vec3>> corners = GetRectangleCorners(vertices);
				std::vector<glm::vec3> cornersToCentre;
				for (std::pair<size_t, glm::vec3> pair : corners) {
					cornersToCentre.push_back(centre - pair.second);
				}
				size_t bottomLeft; // (0, H * scale)
				size_t bottomRight; // (W * scale, H * scale)
				size_t topLeft; // (0, 0)
				size_t topRight; // (W * scale, 0)
				const glm::vec3 unitX(1, 0, 0);
				const glm::vec3 unitY(0, 1, 0);
				const glm::vec3 unitZ(0, 0, 1);
				for (size_t i = 0; i < cornersToCentre.size(); i++) {
					float projX = glm::dot(cornersToCentre[i], unitX);
					float projY = glm::dot(cornersToCentre[i], unitY);
					float projZ = glm::dot(cornersToCentre[i], unitZ);
					// TODO: This will have issues if projX or projY is zero, i.e. the rectangle is parallel to the XZ or YZ plane
					if (projY == 0) {
						std::swap(projZ, projY);
					} 
					if (projX > 0) {
						// Left side
						if (projY > 0) {
							bottomLeft = i;
						}
						else {
							topLeft = i;
						}
					}
					else {
						// Right side
						if (projY > 0) {
							bottomRight = i;
						} else {
							topRight = i;
						}
					}
				}
				int height = textureMap.size();
				int width = textureMap[0].size();
				int w = (int)(width * scale);
				int h = (int)(height * scale);
				vertexMapping.push_back({ bottomLeft, {0, h} });
				vertexMapping.push_back({ bottomRight, {w, h} });
				vertexMapping.push_back({ topLeft, {0, 0} });
				vertexMapping.push_back({ topRight, {w, 0} });
				break;
			}
			default: {
				break;
			}
		}
	}
};

struct Mesh {
	std::vector<MeshPropertyUpdate> propertyUpdates = {};
	std::vector<MeshVec3PropertyUpdate> vec3PropertyUpdates= {};
	std::vector<Triangle> triangles;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> vertexNormals;
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

	void SetTextureMapper(TextureMap textureMap, TextureMap normalMap) {
		this->texture = TextureMapper(textureMap, normalMap, this->triangles[0].normal, TextureMappingPrimative::RECTANGLE);
		this->texture.VertexMapping(vertices, triangles);
	}

	void SetTextureMapper(TextureMap map) {
		this->texture = TextureMapper(map, TextureMappingPrimative::RECTANGLE);
		this->texture.VertexMapping(vertices, triangles);
	}

	const glm::vec3& GetVertex(const size_t t, const size_t v) const {
		return vertices[triangles[t].vertices[v]];
	}

	void SetVertex(const size_t t, const size_t v, const glm::vec3 vertex) {
		vertices[triangles[t].vertices[v]] = vertex;
	}

	// Vertex and face normals will require recomputation if I add rotations
	const glm::vec3& GetVertexNormal(size_t t, size_t v) const {
		return vertexNormals[triangles[t].vertices[v]];
	}

	const Colour GetColour() const {
		return colour;
	}

	void UpdateMesh() {
		std::vector<MeshPropertyUpdate> updaters;
		for (MeshPropertyUpdate& propertyUpdate : propertyUpdates) {
			if (propertyUpdate.Update()) {
				updaters.push_back(propertyUpdate);
			}
		}
		this->propertyUpdates = updaters;

		std::vector<MeshVec3PropertyUpdate> vec3Updaters;
		for (MeshVec3PropertyUpdate& propertyUpdate : vec3PropertyUpdates) {
			if (propertyUpdate.Update()) {
				vec3Updaters.push_back(propertyUpdate);
			}
		}
		this->vec3PropertyUpdates = vec3Updaters;

		// --------------------------------------------------------------------------
		// VERTEX UPDATES:
		bool updated = false;

		// Position = Previous Position, Translation = New Position after translation
		glm::vec3 movement = translation - position;
		if (glm::dot(movement, movement) > 0.000001f) {
			Translate(movement);
			position += movement;
			updated = true;
		}

		if (fabs(glm::dot(scaleTransition, scaleTransition) - 3) > 0.000001f) {
			Scale(InvertVec3(currentScale));
			Scale(scaleTransition);
			updated = true;
			currentScale = scaleTransition;
		}

		if (updated)
			ComputeBoundingBoxes();
	}

	std::pair<glm::vec3, glm::vec3> MeshBoundingBox() {
		glm::vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
		glm::vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);

		for (Triangle& t : triangles) {
			for (size_t axis = 0; axis < 3; axis++) {
				vec3SetAt(min, axis, fmin(vec3At(min, axis), vec3At(t.bbox.first, axis)));
				vec3SetAt(max, axis, fmax(vec3At(max, axis), vec3At(t.bbox.second, axis)));
			}
		}
		return std::make_pair(min, max);
	}

	glm::vec3 Centre() {
		return (MeshBoundingBox().first + MeshBoundingBox().second) / 2.0f;
	}

	void Translate(glm::vec3 translation) {
		for (glm::vec3& vertex : vertices) {
			vertex += translation;
		}
	}

	void Scale(glm::vec3 scale) {
		glm::vec3 centre = this->Centre();
		this->Translate(-centre);
		for (glm::vec3& vertex : vertices) {
			vertex = ScaleVec3(vertex, scale);
		}
		this->Translate(centre);
	}

	void ComputeBoundingBoxes() {
		for (Triangle& t : triangles) {
			t.BoundingBox(vertices);
		}
	}
};

struct VertexNormal {
	glm::vec3 total;
	int n;
	glm::vec3 vertex;
};

// Mesh Container
// Material Properties:
// 		- Smoothness: 		0 = Diffuse, 1 = Reflective
//  	- Refractive Index: 1 = No Refraction
// 		- Transparancy: 	0 = Opaque, 1 = Transparent
// TODO: Just remove this struct its a pointless intermediary for Mesh
struct OBJ {
	Mesh mesh;
	std::string name;
	OBJ(std::vector<glm::vec3> coords, std::vector<std::vector<int>> faces, std::string name, Colour material, size_t coordsStart) : 
			name(name) {
		mesh.colour = material;
		mesh.vertices = {};
		mesh.smoothness = 0;
		mesh.refractiveIndex = 1;
		mesh.name = name;

		for (size_t c = coordsStart; c < coords.size(); c++) {
			mesh.vertices.push_back(coords[c]);
		}

		for (std::vector<int> face : faces) {
			glm::vec3 normal = glm::normalize(glm::cross(coords[face[0] - 1] - coords[face[1] - 1], coords[face[0] - 1] - coords[face[2] - 1]));
			Triangle t(face[0] - 1 - coordsStart, face[1] - 1 - coordsStart, face[2] - 1 - coordsStart, normal);
			t.BoundingBox(mesh.vertices);
			mesh.triangles.push_back(t);
		}

		std::vector<VertexNormal> vertexNormals;
		for (size_t t = 0; t < mesh.triangles.size(); t++) {
			glm::vec3& faceNormal = mesh.triangles[t].normal; 
			for (size_t v = 0; v < 3; v++) {
				const glm::vec3& vertex = mesh.GetVertex(t, v);
				// First check if v is in vertexNormals, if not then we add it with the current face normal
				// If v is in vertexNormals then take the dot product of the current vertex normal and v, if its negative, negate v and add to vertex normals
				bool found = false;
				for (VertexNormal& vn : vertexNormals) {
					if (vn.vertex == vertex) {
						if (glm::dot(vn.total, faceNormal) < 0)
							faceNormal *= -1;
						vn.total += faceNormal;
						vn.n++;
						found = true;
						break;
					}
				}
				if (!found) {
					vertexNormals.push_back({ faceNormal, 1, vertex });
				}
			}
		}
		for (VertexNormal& vn : vertexNormals) {
			vn.total = vn.total * (1.0f / vn.n);
		}
		for (size_t v = 0; v < mesh.vertices.size(); v++) {
			for (VertexNormal& vn : vertexNormals) {
				if (mesh.vertices[v] == vn.vertex) {
					mesh.vertexNormals.push_back(glm::normalize(vn.total));
					break;
				}
			}
		}
	}
};

struct OBJFile {
	std::vector<OBJ> objects;
	std::map<std::string, Colour> materials;

	OBJFile(std::string filename, std::string objfolder, float scaling) {
		std::ifstream File(objfolder + "/" + filename);
		std::string line;
		if (File.is_open()) {
			std::string matfilename;
			getline(File, line);
			matfilename = Split(line, " ")[1];

			// Read Materials
			std::string matpath = objfolder + "/" + matfilename;
			std::ifstream MaterialFile(matpath);
			if (MaterialFile.is_open()) {
				while(getline(MaterialFile, line)) {
					std::vector<std::string> words = Split(line, " ");
					if (words[0] == "newmtl") {
						std::string colourName = words[1];
						getline(MaterialFile, line);
						words = Split(line, " ");
						materials[colourName] = Colour(std::stof(words[1]) * 255, std::stof(words[2]) * 255, std::stof(words[3]) * 255);
					}
				}
				MaterialFile.close();
			}
			else {
				std::cout << "Filename: " << matfilename << ", not found in " << objfolder << "/" << std::endl;
				return;
			}

			// Read Objects
			std::vector<glm::vec3> coords;
			int coordsCount = 0;
			while (getline(File, line)) {
				std::vector<std::string> words = Split(line, " ");
				if (words[0] == "o") { // start parsing new object
					std::string name = words[1];

					getline(File, line);
					std::string material = Split(line, " ")[1];
					getline(File, line);
					words = Split(line, " ");
					int meshCoordCount = 0;
					while (words[0] == "v" ) {
						glm::vec3 coord = { std::stof(words[1]), std::stof(words[2]), std::stof(words[3]) };
						coords.push_back(coord * scaling);
						getline(File, line);
						words = Split(line, " ");
						meshCoordCount++;
					}
					RemoveChar(line, '/');
					words = Split(line, " ");
					std::vector<std::vector<int>> faces;
					while (words[0] == "f") {
						std::vector<int> face = { std::stoi(words[1]), std::stoi(words[2]), std::stoi(words[3]) };
						faces.push_back(face);
						getline(File, line);
						RemoveChar(line, '/');
						words = Split(line, " ");
					}
					objects.push_back(OBJ(coords, faces, name, materials[material], coordsCount));
					coordsCount += meshCoordCount;
				}
			}
			File.close();
		} 
		else {
			std::cout << "Filename: " << filename << ", not found in " << objfolder << "/" << std::endl;
		}

		// Filter triangles for debugging
		/*
		std::vector<OBJ> filteredObjects;
		for (size_t i = 0; i < objects.size(); i++) {
			OBJ object = objects[i];
			if (object.name == "back_wall") {
				filteredObjects.push_back(object);
			}
		}
		objects = filteredObjects;
		*/
	}

	std::vector<Mesh> GetMeshes() {
		std::vector<Mesh> meshes;
		for (OBJ obj : objects) {
				meshes.push_back(obj.mesh);
		}
		return meshes;
	}
};

double Lerp(const double a, const double b, const double s) {
	return a + (b - a) * s;
}

Colour LerpColour(Colour c0, Colour c1, double a) {
	// I wish Colour had operator overloading :(
	int red = c0.red + (c1.red - c0.red) * a;
	int green = c0.green + (c1.green - c0.green) * a;
	int blue = c0.blue + (c1.blue - c0.blue) * a;
	return Colour(red, green, blue);
}

std::vector<float> InterpolateSingleFloats(const float from, const float to, const int steps) {
	if (steps <= 1) {
		return { from };
	}
	const float step = (to - from) / (steps - 1); 
	std::vector<float> results;
	for (int i = 0; i < steps; i++) {
		results.push_back(from + (step * i));
	}
	return results;
}

std::vector<glm::vec3> InterpolateThreeElementValues(const glm::vec3 from, const glm::vec3 to, const int steps) {
	glm::vec3 dir = (to - from) * (1.0f / (steps - 1)); 
	std::vector<glm::vec3> results;
	for (int i = 0; i < steps; i++) {
		results.push_back(from + (dir * static_cast<float>(i)));
	}
	return results;
}

void InterpolateSingleFloatsTest() {
	std::vector<float> result;
	result = InterpolateSingleFloats(2.2, 8.5, 7);
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
}

void InterpolateThreeElementsTest() {
	std::vector<glm::vec3> result;
	result = InterpolateThreeElementValues(glm::vec3(1.0, 4.0, 9.2), glm::vec3(4.0, 1.0, 9.8), 4);
	for (size_t i = 0; i < result.size(); i++) std::cout << result[i].r << ", " << result[i].g << ", " << result[i].b << std::endl;
}

uint32_t PackColour(const uint8_t r, const uint8_t g, const uint8_t b) {
	return (255 << 24) + (r << 16) + (g << 8) + b;
}

glm::vec2 CanvasPointToVec2(CanvasPoint p) {
	return { p.x, p.y };
}

Shape2D CreateLine2D(CanvasPoint from, CanvasPoint to, Colour c) {
	float xdiff = to.x - from.x;
	float ydiff = to.y - from.y;
	float steps = fmax(fabs(xdiff), fabs(ydiff));
	float xstep = xdiff / steps;
	float ystep = ydiff / steps;
	std::vector<CanvasPoint> ps;
	for (float i = 0.0; i < steps; i++) {
		float x = from.x + xstep * i;
		float y = from.y + ystep * i;
		ps.push_back(CanvasPoint(x, y));
	}
	return { ps, c };
}

Shape2D CreateStrokedTriangle2D(CanvasTriangle verticies, Colour c) {
	CanvasPoint v0 = verticies.v0();
	CanvasPoint v1 = verticies.v1();
	CanvasPoint v2 = verticies.v2();

	Shape2D l1 = CreateLine2D(v0, v1, c);
	Shape2D l2 = CreateLine2D(v1, v2, c);
	Shape2D l3 = CreateLine2D(v2, v0, c);
	std::vector<CanvasPoint> ps;
	ps.insert(ps.end(), l1.points.begin(), l1.points.end());
	ps.insert(ps.end(), l2.points.begin(), l2.points.end());
	ps.insert(ps.end(), l3.points.begin(), l3.points.end());
	return { ps, c };
}

glm::vec3 CanvasPointToVec3(CanvasPoint p) {
	return { p.x, p.y, 0 };
}

Shape2D CreateFilledTriangle2D(CanvasTriangle verticies, const Colour c, const ModelTriangle* pModelTriangle, std::vector<std::vector<DepthPoint>>* const pDepthBuffer) {
	CanvasPoint v0 = verticies.v0();
	CanvasPoint v1 = verticies.v1();
	CanvasPoint v2 = verticies.v2();
	std::vector<CanvasPoint> ps;
	int maxX = fmax(v0.x, fmax(v1.x, v2.x));
	int minX = fmin(v0.x, fmin(v1.x, v2.x));
	int maxY = fmax(v0.y, fmax(v1.y, v2.y));
	int minY = fmin(v0.y, fmin(v1.y, v2.y));
	glm::vec3 A = CanvasPointToVec3(v0);
	glm::vec3 B = CanvasPointToVec3(v1);
	glm::vec3 C = CanvasPointToVec3(v2);

	float areaScalar = glm::length(glm::cross(B - A, C - A));
	for (int x = minX; x <= maxX; x++) {
  		for (int y = minY; y <= maxY; y++) {
			glm::vec3 P(x, y, 0);
			glm::vec3 PC = C - P;
			glm::vec3 PB = B - P;
			glm::vec3 PA = A - P;
			float alpha = glm::length(glm::cross(PB, PC)) / areaScalar;
			float beta = glm::length(glm::cross(PC, PA)) / areaScalar;
			float gamma = glm::length(glm::cross(PA, PB)) / areaScalar;
			if ((alpha >= 0 && alpha <= 1) && (beta >= 0 && beta <= 1) && (gamma >= 0 && gamma <= 1) && (fabs(alpha + beta + gamma - 1) <= 0.01)) {
				ps.push_back(CanvasPoint(x, y));
				if (pModelTriangle && pDepthBuffer) {
					// Convert barcyentric to cartesian
					std::vector<std::vector<DepthPoint>>& depthBuffer = *pDepthBuffer;
					glm::vec3 coords3d = (pModelTriangle->vertices[0] * alpha) + (pModelTriangle->vertices[1] * beta) + (pModelTriangle->vertices[2] * gamma);
					float depth = -1 / coords3d.z;
					if (y < HEIGHT && y >= 0 && x < WIDTH && x > 0) {
						if (depthBuffer[y][x].depth < depth) {
							depthBuffer[y][x] = { c, depth };
						}
					}
				}
			}
		}
	}
	return { ps, c };
}

std::pair<Shape2D, Shape2D> FilledTriangle2DWithOutline(CanvasTriangle verticies, Colour c) {
	Shape2D shaded = CreateFilledTriangle2D(verticies, c, nullptr, nullptr);
	Shape2D outline = CreateStrokedTriangle2D(verticies, Colour(255, 255, 255));
	return std::make_pair(shaded, outline);
}

void GrayscaleGradient(DrawingWindow &window) {
	const std::vector<float> row = InterpolateSingleFloats(255, 0, window.width);
	std::vector<uint32_t> rowColours;
	std::transform(row.cbegin(), row.cend(), std::back_inserter(rowColours), [](float f) -> uint32_t { return PackColour(f, f, f); });
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			window.setPixelColour(x, y, rowColours.at(x));
		}
	}
}

void TwoDGradient(DrawingWindow &window) {
	const glm::vec3 topLeft(255, 0, 0);        // red 
	const glm::vec3 topRight(0, 0, 255);       // blue 
	const glm::vec3 bottomRight(0, 255, 0);    // green 
	const glm::vec3 bottomLeft(255, 255, 0);   // yellow

	const std::vector<glm::vec3> left = InterpolateThreeElementValues(topLeft, bottomLeft, window.height);
	const std::vector<glm::vec3> right = InterpolateThreeElementValues(topRight, bottomRight, window.height);
	
	for (size_t y = 0; y < window.height; y++) {
		const std::vector<glm::vec3> row = InterpolateThreeElementValues(left[y], right[y], window.width);
		for (size_t x = 0; x < window.width; x++) {
			window.setPixelColour(x, y, PackColour(row[x].r, row[x].g, row[x].b));
		}
	}
}

CanvasPoint RandomCanvasPoint() {
	return CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
}

Colour RandomColour() {
	return Colour(rand() % 256, rand() % 256, rand() % 256);
}

void AddRandomTriangle(const bool filled, std::vector<Shape2D> &shapes) {
	if (filled) {
		std::pair<Shape2D, Shape2D> filled = FilledTriangle2DWithOutline(
			CanvasTriangle(
				RandomCanvasPoint(),
				RandomCanvasPoint(),
				RandomCanvasPoint()),
			RandomColour());
		shapes.push_back(filled.first);
		shapes.push_back(filled.second);
	} 
	else {
		shapes.push_back(CreateStrokedTriangle2D(
			CanvasTriangle(
				RandomCanvasPoint(),
				RandomCanvasPoint(),
				RandomCanvasPoint()),
			RandomColour()));
	}
}

void DrawPoints(DrawingWindow &window, const std::vector<glm::vec2> points, const int colour) {
	for (const glm::vec2 p : points) {
		window.setPixelColour(p.x, p.y, colour);
	}
}

void DrawTexturedShape2D(DrawingWindow &window, const Shape2D shape, const TextureMap texture) {
	for (CanvasPoint p : shape.points) {
		const int pixelIndex = p.texturePoint.y * texture.width + p.texturePoint.x;
		if (pixelIndex < texture.pixels.size()) {
			const uint32_t col = texture.pixels.at(pixelIndex);
			window.setPixelColour(p.x, p.y, col);
		}
	}
}

void DrawShape2D(DrawingWindow &window, Shape2D shape) {
	int col = PackColour(shape.colour.red, shape.colour.green, shape.colour.blue);
	for (CanvasPoint p : shape.points) {
		window.setPixelColour(p.x, p.y, col);
	}
}

void WitchSymbol(DrawingWindow &window) {
	float midX = WIDTH / 2;
	float midY = HEIGHT / 2;
	Colour white(255, 255, 255);
	DrawShape2D(window, CreateLine2D(CanvasPoint(0, 0), CanvasPoint(midX, midY), white));
	DrawShape2D(window, CreateLine2D(CanvasPoint(WIDTH - 1, 0), CanvasPoint(midX, midY), white));
	DrawShape2D(window, CreateLine2D(CanvasPoint(midX, 0), CanvasPoint(midX, HEIGHT), white));
	DrawShape2D(window, CreateLine2D(CanvasPoint(midX - WIDTH / 6, midY), CanvasPoint(midX + WIDTH / 6, midY), white));
}

glm::mat3 RotationMatrix(const float x, const float y, const float z) {
	//  Roll
	//  | cos(x) | -sin(x) | 0 |
	//  | sin(x) |  cos(x) | 0 |    = Rx(x)
	//  |   0    |    0    | 1 |
	//
	//  Yaw
	//  |  cos(y) | 0 | sin(y) |
	//  |    0    | 1 |   0    |    = Ry(y)
	//  | -sin(y) | 0 | cos(y) |
	//
	//  Pitch
	//  | 1 |   0    |    0    |
	//  | 0 | cos(z) | -sin(z) |    = Rz(z)
	//  | 0 | sin(z) |  cos(z) |
	//
	// R = Rx(x) * Ry(y) * Rz(z)
	const glm::mat3 roll(glm::vec3(cos(x), sin(x), 0), glm::vec3(-sin(x), cos(x), 0), glm::vec3(0, 0, 1));
	const glm::mat3 yaw(glm::vec3(cos(y), 0, -sin(y)), glm::vec3(0, 1, 0), glm::vec3(sin(y), 0, cos(y)));
	const glm::mat3 pitch(glm::vec3(1, 0, 0), glm::vec3(0, cos(z), sin(z)), glm::vec3(0, -sin(z), cos(z)));
	return yaw * pitch * roll;
}

struct Camera {
	bool orbiting = false;
	glm::vec3 orbitCentre;
	float orbitRadius;
	float orbitAngle;

	glm::vec3 position;
	glm::vec3 rotation;
	glm::mat3 rotationMatrix;
	float focalLength;

	Camera(const glm::vec3 p, const glm::vec3 r, const float f) {
		position = p;
		rotation = r;
		focalLength = f;
		rotationMatrix = RotationMatrix(rotation.x, rotation.y, rotation.z);
	}

	// X = Roll, Y = Yaw, Z = Pitch
	void AddRotation(float x, float y, float z) {
		rotation.x += x;
		rotation.y += y;
		rotation.z += z;
		rotationMatrix = RotationMatrix(rotation.x, rotation.y, rotation.z);
	}

	void AddTranslation(float x, float y, float z) {
		glm::vec3 delta = RotationMatrix(-rotation.x, -rotation.y, -rotation.z) * glm::vec3(x, y, z);
		position += delta;
	}

	void lookAt(glm::vec3 p) {
		glm::vec3 cameraToP = p - position;
		float pitch = atan2(-cameraToP.y, fabs(cameraToP.z));
		float yaw = atan2(cameraToP.x, -cameraToP.z);

		this->rotation = { 0, yaw, pitch };
		rotationMatrix = RotationMatrix(rotation.x, rotation.y, rotation.z);
	}

	void startOrbit(glm::vec3 centre, float radius) {
		position = centre;
		position.z += radius;
		orbitRadius = radius;
		rotation = glm::vec3(0, 0, 0);
		rotationMatrix = RotationMatrix(0, 0, 0);
		orbitAngle = 0;
		orbiting = true;
	}

	void stopOrbit() {
		orbiting = false;
	}

	void orbitStep() {
		orbitAngle += M_PI / 180; // Add 1 degree
		position.x = orbitRadius * sin(orbitAngle) + orbitCentre.x;
		position.z = orbitRadius * cos(orbitAngle) + orbitCentre.z;
		lookAt(orbitCentre);
	}

	void output() {
		std::cout << "Rotation: " << glm::to_string(rotation) << " Position: " << glm::to_string(position) << std::endl;
	}
};

CanvasPoint getCanvasIntersectionPoint(const Camera &camera, glm::vec3 vertexPosition) {
	// transform the vertex such that the camera is the origin
	glm::vec3 vPos = camera.rotationMatrix * (vertexPosition - camera.position);
	const float u = (150 * camera.focalLength * vPos.x / -vPos.z) + WIDTH / 2.0;
	const float v = (150 * camera.focalLength * vPos.y / vPos.z) + HEIGHT / 2.0;
	return CanvasPoint(u, v);
}

glm::vec3 CanvasPointToWorld(const Camera &camera, const float u, const float v) {
	const float x = (u - (WIDTH / 2.0)) / 150.0;
	const float y = -(v - (HEIGHT / 2.0)) / 150.0;
	return glm::transpose(camera.rotationMatrix) * (glm::vec3(x, y, -camera.focalLength));
}


Shape2D Pointcloud(std::vector<Mesh>& meshes, Camera camera) {
	std::vector<CanvasPoint> ps;
	for (Mesh& m : meshes) {
		for (size_t t = 0; t < m.triangles.size(); t++) {
			for (size_t v = 0; v < 3; v++) {
				ps.push_back(getCanvasIntersectionPoint(camera, m.GetVertex(t, v)));
			}
		}
	}

	return { ps, Colour(255, 255, 255) };
}

std::vector<Shape2D> Wireframe(std::vector<Mesh>& meshes, Camera& camera) {
	std::vector<Shape2D> wireframe;
	for (Mesh& m : meshes) {
		for (size_t t = 0; t < m.triangles.size(); t++) {
			CanvasPoint p1 = getCanvasIntersectionPoint(camera, m.GetVertex(t, 0));
			CanvasPoint p2 = getCanvasIntersectionPoint(camera, m.GetVertex(t, 1));
			CanvasPoint p3 = getCanvasIntersectionPoint(camera, m.GetVertex(t, 2));
			Shape2D triangle2D = CreateStrokedTriangle2D(CanvasTriangle(p1, p2, p3), Colour(255, 255, 255));
			wireframe.push_back(triangle2D);
		}
	}
	return wireframe;
}

ModelTriangle TranslateTriangle(const ModelTriangle tri, const glm::vec3 t) {
	return ModelTriangle(tri.vertices[0] + t, tri.vertices[1] + t, tri.vertices[2] + t, tri.colour);
}

ModelTriangle RotateTriangle(const ModelTriangle tri, const glm::mat3 r) {
	return ModelTriangle(r * tri.vertices[0], r * tri.vertices[1], r * tri.vertices[2], tri.colour);
}

bool PointInsideCanvas(CanvasPoint p) {
	return (p.x >= 0 && p.x < WIDTH) && (p.y >= 0 && p.y < HEIGHT);
}

std::vector<std::vector<DepthPoint>> RasterisedRender(DrawingWindow& window, const std::vector<Mesh>& meshes, const Camera camera) {
	std::vector<std::vector<DepthPoint>> depthBuffer;
	for (int y = 0; y < HEIGHT; y++) {
		depthBuffer.push_back(std::vector<DepthPoint>());
		for (int x = 0; x < WIDTH; x++)
			depthBuffer[y].push_back({ Colour(), 0 });
	}

	for (const Mesh& m : meshes) {
		for (size_t t = 0; t < m.triangles.size(); t++) {
			// We must maintain a depth buffer as we build the triangles
			const CanvasPoint p1 = getCanvasIntersectionPoint(camera, m.GetVertex(t, 0));
			const CanvasPoint p2 = getCanvasIntersectionPoint(camera, m.GetVertex(t, 1));
			const CanvasPoint p3 = getCanvasIntersectionPoint(camera, m.GetVertex(t, 2));
			ModelTriangle tri(m.GetVertex(t, 0), m.GetVertex(t, 1), m.GetVertex(t, 2), m.GetColour());
			const ModelTriangle cameraSpaceTriangle = RotateTriangle(TranslateTriangle(tri, -camera.position), camera.rotationMatrix);
			if (PointInsideCanvas(p1) || PointInsideCanvas(p2) || PointInsideCanvas(p3)) {
				CreateFilledTriangle2D(CanvasTriangle(p1, p2, p3), tri.colour, &cameraSpaceTriangle, &depthBuffer);
			}
		}
	}

	return depthBuffer;
}

struct RayCollision {
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec3 barycentric;
	size_t triangle;
	size_t mesh;
	float distanceToCamera;
};

struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;
	bool isInsideMesh;
};

bool TriangleHitLoc(const Mesh& mesh, const int& triangle, const Ray& r, glm::vec3& loc, glm::vec3& barycentric) {
	// Get ray intersection with the plane described by two edges of the triangle

	// TODO: Figure out how to remove this because normal is precomputed
	glm::vec3 normal = glm::cross(mesh.GetVertex(triangle, 1) - mesh.GetVertex(triangle, 0), mesh.GetVertex(triangle, 2) - mesh.GetVertex(triangle, 0));
	float d = glm::dot(mesh.GetVertex(triangle, 0), normal);
	float n = glm::dot(normal, r.origin); // n.p_0
	float m = glm::dot(normal, r.direction); // n.u
	if (m == 0) 
		return false;
	float t = (d - n) / m; // tri = (d - n.p0) / n.u
	loc = r.origin + (r.direction * t);

	// Check that loc is inside triangle
	float area2 = glm::length(normal);
	glm::vec3 PC = mesh.GetVertex(triangle, 2) - loc;
	glm::vec3 PB = mesh.GetVertex(triangle, 1) - loc;
	glm::vec3 PA = mesh.GetVertex(triangle, 0) - loc;
	float alpha = glm::length(glm::cross(PB, PC)) / area2;
	float beta = glm::length(glm::cross(PC, PA)) / area2;
	float gamma = glm::length(glm::cross(PA, PB)) / area2;	
	barycentric = glm::vec3(alpha, beta, gamma);
	return (alpha >= 0 && alpha <= 1) &&
		   (beta  >= 0 && beta  <= 1) &&
		   (gamma >= 0 && gamma <= 1) &&
		   (fabs(alpha + beta + gamma - 1) <= 0.0001);
}

bool SlabTest(const Ray& r, const Mesh& mesh, const Triangle& t) {
	double tx1 = (t.bbox.first.x - r.origin.x) / r.direction.x;
	double tx2 = (t.bbox.second.x - r.origin.x) / r.direction.x;

	double tmin = fmin(tx1, tx2);
	double tmax = fmax(tx1, tx2);

	double ty1 = (t.bbox.first.y - r.origin.y) / r.direction.y;
	double ty2 = (t.bbox.second.y - r.origin.y) / r.direction.y;

	tmin = fmax(tmin, fmin(ty1, ty2));
	tmax = fmin(tmax, fmax(ty1, ty2));
	return tmax >= tmin;
}

void Collision(const Ray& r, const std::vector<Mesh>& meshes, const size_t mesh_i, const size_t triangle_i, RayCollision& nearest) {
	const Mesh& mesh = meshes[mesh_i];
	const Triangle& t = mesh.triangles[triangle_i];
	if (!SlabTest(r, mesh, t)) {
		return;
	}

	glm::vec3 loc;
	glm::vec3 barycentric;
	if (TriangleHitLoc(mesh, triangle_i, r, loc, barycentric)) {
		float sqrDist = glm::length2(loc - r.origin);
		// Check that the collision is after the start of the ray & collision is closer
		if (glm::dot(r.direction, loc - r.origin) > 0.0 && sqrDist < nearest.distanceToCamera) {
			glm::vec3 triNorm = t.normal;
			if (glm::dot(triNorm, r.direction) > 0.0)
				triNorm *= -1;
			nearest = { loc, triNorm, barycentric, triangle_i, mesh_i, sqrDist };
		}
	}
}

// Returning a bool here cause someone decided we are using C++11 and therefore I can't use std::optional :(
bool NearestRayCollision(const Ray& r, const std::vector<Mesh>& meshes, RayCollision& collision) {
	RayCollision nearest = { glm::vec3(), glm::vec3(), glm::vec3(), 0, 0, FLT_MAX };
	for (size_t m = 0; m < meshes.size(); m++) {
		const Mesh& mesh = meshes[m];
		for (size_t i = 0; i < mesh.triangles.size(); i++) {
			Collision(r, meshes, m, i, nearest);
		}
	}
	if (nearest.distanceToCamera == FLT_MAX)
		return false;
	collision = nearest;
	return true; 
}

enum Shading {
	PROXIMITY,
	ANGLE_OF_INCIDENCE,
	SPECULAR,
	PROXY_AOI_AND_SPECULAR,
	GOURAUD,
	PHONG,
	NONE
};

struct Light {
	glm::vec3 position;
	double radius;
};

double ProximityLighting(const std::vector<Light>& lights, const RayCollision& rc) {
	// 1 / 4pir^2
	double dist2 = -FLT_MAX;
	for (const Light& light : lights) {
		dist2 = fmax(dist2, glm::length2(rc.position - light.position));
	}
	return fmin(1 / (4 * M_PI * dist2) * 8, 1);
}

double AngleOfIncidenceLighting(const std::vector<Light>& lights, const RayCollision& rc) {
	double lighting = -FLT_MAX;
	for (const Light& light : lights) {
		glm::vec3 dir = glm::normalize(light.position - rc.position);
		lighting = fmax(lighting, glm::dot(dir, rc.normal));
	}
	return lighting > 0 ? lighting/2 : 0;
}

// Returns the reflected vector given normalised incident and normal vectors
glm::vec3 ReflectionVector(const glm::vec3& incident, const glm::vec3& normal) {
	return incident - 2.0f * normal * (glm::dot(normal, incident));
}

double SpecularLighting(const std::vector<Light>& lights, const RayCollision& rc, const glm::vec3& cameraP) {
	double lighting = -FLT_MAX;
	for (const Light& light : lights) {
		const glm::vec3 colToCamera = glm::normalize(cameraP - rc.position);
		const glm::vec3 incident = glm::normalize(rc.position - light.position);
		const glm::vec3 reflected = ReflectionVector(incident, rc.normal);
		lighting = fmax(lighting, glm::dot(reflected, colToCamera));
	}

	return lighting > 0 ? pow(lighting, 256) : 0;
}

double GouraudShading(const std::vector<Light>& lights, const RayCollision& rc, const glm::vec3& cameraP, const Mesh& mesh, const int& triangle) {
	// First apply AOI lighting + Specular Shading at each vertex, may have to invert normals
	std::vector<double> vertexLighting;
	for (size_t i = 0; i < 3; i++) {
		const glm::vec3& vertex = mesh.GetVertex(triangle, i);
		RayCollision c;
		const glm::vec3& vertexNormal = mesh.GetVertexNormal(triangle, i);
		c.normal = glm::dot(rc.normal, vertexNormal) >= 0 ? vertexNormal : -vertexNormal;
		c.position = vertex;
		vertexLighting.push_back(fmin(1, fmax(SpecularLighting(lights, c, cameraP), AngleOfIncidenceLighting(lights, c))));
	}

	// Assume vertex order remains the same, (A, B, C)
	return vertexLighting[0] * rc.barycentric[0] + vertexLighting[1] * rc.barycentric[1] + vertexLighting[2] * rc.barycentric[2];
}
double PhongShading(const std::vector<Light>& lights, const RayCollision& rc, const glm::vec3& cameraP, const Mesh& mesh, const int& triangle) {
	glm::vec3 interpolatedNormal(0, 0, 0);
	for (size_t i = 0; i < 3; i++) {
		glm::vec3 vertexNormal = mesh.GetVertexNormal(triangle, i);
		vertexNormal = glm::dot(rc.normal, vertexNormal) >= 0 ? vertexNormal : -vertexNormal;
		interpolatedNormal += vertexNormal * rc.barycentric[i];
	}
	RayCollision c;
	c.normal = interpolatedNormal;
	c.position = rc.position;
	const double dist2 = glm::length2(rc.position - lights[0].position);
	return fmin(1, SpecularLighting(lights, c, cameraP) + AngleOfIncidenceLighting(lights, c) / (4 * M_PI * dist2) * 100);/// (4 * M_PI * dist2));
}

enum class ShadowType {
	HARD_SHADOWS,
	SOFT_SHADOWS,
	NONE
};

std::pair<Colour, double> RayCast(const Ray& ray, const Camera& camera, const std::vector<Light>& lights, const std::vector<Mesh>& meshes, const int depth, const ShadowType shadows=ShadowType::NONE, const Shading shader=Shading::NONE);

// Cast a new ray from the collision point of the previous ray
// Note: This is technically recursive if it keeps hitting mirrors
std::pair<Colour, double> MirrorReflect(
										const Ray& incoming, 
										const RayCollision& rc, 
										const Camera& camera, 
										const std::vector<Light>& lights, 
										const std::vector<Mesh>& meshes, 
										const int depth, 
										const ShadowType shadows=ShadowType::HARD_SHADOWS, 
										const Shading shader=Shading::NONE) {
	const glm::vec3 reflected = ReflectionVector(glm::normalize(incoming.direction), rc.normal);
	const Ray reflectedRay = { rc.position + (rc.normal * 0.0001f), reflected, incoming.isInsideMesh };
	return RayCast(reflectedRay, camera, lights, meshes, depth + 1, shadows, shader);
}

double SchlickReflectance(const double cosTheta, const double refractiveIndex) {
	double r0 = (refractiveIndex - 1) / (refractiveIndex + 1);
	r0 = r0 * r0;
	return r0 + (1 - r0) * pow(1 - cosTheta, 5);
}

double ProcessShader(const Shading shader, const RayCollision& rc, const glm::vec3 camera, const std::vector<Light>& lights, const std::vector<Mesh>& meshes) {
	switch (shader) {
		case PROXIMITY: {
			return ProximityLighting(lights, rc);
		}
		case ANGLE_OF_INCIDENCE: {
			return AngleOfIncidenceLighting(lights, rc);
		}
		case SPECULAR: {
			return SpecularLighting(lights, rc, camera);
		}
		case PROXY_AOI_AND_SPECULAR: {
			double dist2 = glm::length2(lights[0].position - rc.position);
			return fmin(1, (SpecularLighting(lights, rc, camera) + AngleOfIncidenceLighting(lights, rc)) / (4 * M_PI * dist2) * 80);/// (4 * M_PI * dist2));
		}
		case GOURAUD: {
			return GouraudShading(lights, rc, camera, meshes[rc.mesh], rc.triangle);
		}
		case PHONG: {
			return PhongShading(lights, rc, camera, meshes[rc.mesh], rc.triangle);
		}
		default: {
			return 1;
		}
	}
}

const Ray ComputeRefractedRay(const double refractiveIndex, const glm::vec3& incoming, const glm::vec3& normal, const glm::vec3& origin, const bool insideMesh) {
	const float cosTheta = fmin(glm::dot(-glm::normalize(incoming), normal), 1.0);
	float refractionRatio = refractiveIndex;
	if (!insideMesh)
		refractionRatio = 1.0 / refractiveIndex;

	const float sinR2 = refractionRatio * refractionRatio * (1.0f - (cosTheta * cosTheta));
	const float cosR = sqrt(1.0f - sinR2);
	const glm::vec3 refractedVector = refractionRatio * incoming + (refractionRatio * cosTheta - cosR) * normal;
	return { origin - (normal * 0.0001f), refractedVector, !insideMesh};
}

std::pair<Colour, double> DielectricTransmission(
												 const Ray& incoming,
												 const RayCollision& rc, 
												 const Camera& camera, 
												 const std::vector<Light>& lights, 
												 const std::vector<Mesh>& meshes, 
												 const int depth, 
												 const ShadowType shadows=ShadowType::HARD_SHADOWS, 
												 const Shading shader=Shading::NONE) {
	glm::vec3 normal = rc.normal;
	const double transparancy = meshes[rc.mesh].transparancy;
	const double refractiveIndex = meshes[rc.mesh].refractiveIndex;
	const double smoothness = meshes[rc.mesh].smoothness;

	const float cosTheta = fmin(glm::dot(-glm::normalize(incoming.direction), normal), 1.0);
	const double fresnelReflectance = SchlickReflectance(cosTheta, refractiveIndex);
	const double reflectionScalar = (1 - transparancy) * (smoothness) + transparancy * fresnelReflectance;
	const double transmissionScalar = transparancy * (1 - fresnelReflectance);
	const double localScalar = (1 - transparancy) * (1 - smoothness);

	const Ray refractedRay = ComputeRefractedRay(refractiveIndex, incoming.direction, normal, rc.position, incoming.isInsideMesh);

	const std::pair<Colour, double> reflected = MirrorReflect(incoming, rc, camera, lights, meshes, depth, shadows, shader);
	const std::pair<Colour, double> transmitted = RayCast(refractedRay, camera, lights, meshes, depth + 1, shadows, shader);
	const Colour localColour = meshes[rc.mesh].colour;
	const double localLighting = ProcessShader(shader, rc, camera.position, lights, meshes);

	const Colour combinedColour = AddColours(ScaleColour(transmitted.first, transmissionScalar), AddColours(ScaleColour(localColour, localScalar), ScaleColour(reflected.first, reflectionScalar)));
	const double combinedLighting = transmitted.second * transmissionScalar + localLighting * localScalar + reflected.second * reflectionScalar;
	return std::make_pair(combinedColour, combinedLighting);
}

std::pair<glm::vec3, glm::vec3> OrthornormalPlaneVectors(const glm::vec3 normal) {
	glm::vec3 offset = { normal.x + 0.1, normal.y, normal.z };
	glm::vec3 v1 = glm::cross(normal, offset);
	glm::vec3 v2 = glm::cross(normal, v1);
	return std::make_pair(glm::normalize(v1), glm::normalize(v2));
}

double InShadow(const glm::vec3& p, const glm::vec3& norm, const std::vector<Mesh>& meshes, const std::vector<Light>& lights, ShadowType type) {
	double min = 1;
	for (const Light& light : lights) {
		glm::vec3 startPosition = p + norm * 0.0001f;
		glm::vec3 v = light.position - p;
		if (type == ShadowType::SOFT_SHADOWS) {
			int occludedCount = 0;
			int lightCount = 0;
			const int steps = 10;
			float stepSize = light.radius / steps;
			std::pair<glm::vec3, glm::vec3> plane = OrthornormalPlaneVectors(v);

			for (int i = -steps / 2; i < steps / 2; i++) {
				for (int f = -steps / 2; f < steps / 2; f++) {
					glm::vec3 l = light.position + (i * stepSize * plane.first) + (f * stepSize * plane.second);
					if (glm::length2(l - light.position) < light.radius * light.radius) {
						Ray r = { startPosition, l - p};
						RayCollision c;
						lightCount++;
						if (NearestRayCollision(r, meshes, c))
							occludedCount++;
					}
				}
			}
			min = fmin(min, (float)occludedCount / lightCount);
		}
		else {
			Ray r = { startPosition, v, false };
			RayCollision c;
			int bounces = 0;
			double cummulativeOcclusion = 1;
			size_t prevMesh = -1;
			while(NearestRayCollision(r, meshes, c) && bounces < 10) {
				// If its a dielectric bounce up to 10 times
				bounces++;
				if (meshes[c.mesh].transparancy == 0)
					return 1;
				
				if (prevMesh != c.mesh)
					cummulativeOcclusion *= (1 - meshes[c.mesh].transparancy);
				prevMesh = c.mesh;
				const double refractiveIndex = meshes[c.mesh].refractiveIndex;
				r = ComputeRefractedRay(refractiveIndex, r.direction, c.normal, c.position, !r.isInsideMesh);
			}
			return cummulativeOcclusion == 1 ? 0 : cummulativeOcclusion;
		}
	}
	return min;
}

std::pair<Colour, double> RayCast(const Ray& ray, const Camera& camera, const std::vector<Light>& lights, const std::vector<Mesh>& meshes, const int depth, const ShadowType shadows, const Shading shader) {
	RayCollision c;
	Colour colour;
	double lighting = 0;
	constexpr int maxdepth = 10;
	constexpr float ambient = 0.2f;
	if (NearestRayCollision(ray, meshes, c)) {
		const Mesh& mesh = meshes[c.mesh];
		
		// Get colour of collision, either the colour of the mesh or of the texture map 
		if (mesh.texture.type != TextureMappingType::NONE) {
			// First map the world collision point to the texture map.
			// If the texture map primative is a rectangle
			colour = mesh.texture.BarycentricTexture(mesh.triangles[c.triangle], c.barycentric);
			if (mesh.texture.type == TextureMappingType::NORMAL)
				c.normal = mesh.texture.BarycentricNormal(mesh.triangles[c.triangle], c.barycentric);
		}
		else {
			colour = mesh.colour;
		}

		// Compute the direct point illumation shadows
		if (shadows != ShadowType::NONE){
			double occluded = InShadow(c.position, c.normal, meshes, lights, shadows);
			lighting = (ambient - 1) * occluded + 1;
		}

		// Compute the shader lighting, multiply by direct point illumination shadows
		lighting = ProcessShader(shader, c, camera.position, lights, meshes) * lighting;

		// If the material is refractive then compute the colour and lighting generated from the dielectric transmission 
		if (mesh.refractiveIndex != 1 && depth < maxdepth) {
			std::pair<Colour, double> refraction = DielectricTransmission(ray, c, camera, lights, meshes, depth, shadows, shader);
			colour = refraction.first;
			lighting = refraction.second;
		}
		else if (mesh.smoothness > 0 && depth < maxdepth) {
			// Reflective Material
			std::pair<Colour, double> reflection = MirrorReflect(ray, c, camera, lights, meshes, depth, shadows, shader);
			colour = LerpColour(colour, reflection.first, mesh.smoothness);
			lighting = reflection.second;
		}
	}
	else
		colour = Colour(0, 0, 0);
	lighting = fmax(lighting, ambient);
	return std::make_pair(colour, lighting);
}

void Raytrace(DrawingWindow& window, const Camera& camera, const std::vector<Light>& lights, const std::vector<Mesh>& meshes, const ShadowType shadows=ShadowType::HARD_SHADOWS, const Shading shader=Shading::NONE) {
	for (int v = 0; v < HEIGHT; v++) {
		for (int u = 0; u < WIDTH; u++) {
			const glm::vec3 rayV = CanvasPointToWorld(camera, u, v);
			const Ray ray = { camera.position, rayV, false };
			const std::pair<Colour, double> pixel = RayCast(ray, camera, lights, meshes, 0, shadows, shader);
			const Colour colour = pixel.first;
			const double lighting = pixel.second;
			window.setPixelColour(u, v, PackColour(colour.red * lighting, colour.green * lighting, colour.blue * lighting));
		}
	}
}

enum Drawing {
	GRAYSCALE_GRADIENT,
	TWO_D_GRADIENT,
	WITCH_SYMBOL,
	RANDOM_TRIANGLES,
	POINT_CLOUD,
	WIRE_FRAME,
	RASTERISED_3D,
	RAYTRACE
};

bool handleEvent(const SDL_Event event, DrawingWindow &window, std::vector<Shape2D> &shapes, Drawing& d, Camera* const pCamera, std::vector<Light>& lights) {
	Drawing Drawings3D[] = { POINT_CLOUD, WIRE_FRAME, RASTERISED_3D, RAYTRACE };
	if (event.type == SDL_KEYDOWN) {
	 	if (event.key.keysym.sym == SDLK_u && d == RANDOM_TRIANGLES) AddRandomTriangle(false, shapes);
		else if (event.key.keysym.sym == SDLK_f && d == RANDOM_TRIANGLES) AddRandomTriangle(true, shapes);
		else if (pCamera) {
			Camera& c = *pCamera;
			constexpr float step = 0.1;
			constexpr float rotStep = M_PI / 180;
			// Translations
			if (event.key.keysym.sym == SDLK_w) c.AddTranslation(0, 0, -step);
			else if (event.key.keysym.sym == SDLK_s) c.AddTranslation(0, 0, step);
			else if (event.key.keysym.sym == SDLK_a) c.AddTranslation(-step, 0, 0);
			else if (event.key.keysym.sym == SDLK_d) c.AddTranslation(step, 0, 0);
			else if (event.key.keysym.sym == SDLK_SPACE) c.AddTranslation(0, step, 0);
			else if (event.key.keysym.sym == SDLK_c) c.AddTranslation(0, -step, 0);

			// Rotations
			else if (event.key.keysym.sym == SDLK_LEFT) c.AddRotation(0, -rotStep, 0);
			else if (event.key.keysym.sym == SDLK_RIGHT) c.AddRotation(0, rotStep, 0);
			else if (event.key.keysym.sym == SDLK_UP) c.AddRotation(0, 0, rotStep);
			else if (event.key.keysym.sym == SDLK_DOWN) c.AddRotation(0, 0, -rotStep);

			// Light Translations
			else if (event.key.keysym.sym == SDLK_j) lights[0].position += glm::vec3(-step, 0, 0);
			else if (event.key.keysym.sym == SDLK_l) lights[0].position += glm::vec3(step, 0, 0);
			else if (event.key.keysym.sym == SDLK_k) lights[0].position += glm::vec3(0, 0, step);
			else if (event.key.keysym.sym == SDLK_i) lights[0].position += glm::vec3(0, 0, -step);
			else if (event.key.keysym.sym == SDLK_m) lights[0].position += glm::vec3(0, step, 0);
			else if (event.key.keysym.sym == SDLK_n) lights[0].position += glm::vec3(0, -step, 0);
		
			else if (event.key.keysym.sym == SDLK_l) c.lookAt(glm::vec3(0, 0, 0));
			else if (event.key.keysym.sym == SDLK_o) {
				if (c.orbiting) {
					c.stopOrbit();
				} else {
					c.startOrbit(glm::vec3(0, 0, 0), 4);
				}
			}
			else if (event.key.keysym.sym == SDLK_q) {
				for (int i = 0; i < 4; i++) {
					if (Drawings3D[i] == d) {
						d = Drawings3D[(i + 1) % 4];
						break;
					}
				}
			}
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		std::cout << "(" << event.button.x << ", " << event.button.y << ")" << std::endl;
	} else if (event.type == SDL_QUIT) {
		return true;
	}
	return false;
}

void DrawRasterized3D(DrawingWindow& window, Camera& camera, std::vector<Mesh>& meshes) {
	std::vector<std::vector<DepthPoint>> depthBuffer = RasterisedRender(window, meshes, camera);
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			Colour c  = depthBuffer[y][x].colour;
			window.setPixelColour(x, y, PackColour(c.red, c.green, c.blue));
		}
	}
}

void draw2D(DrawingWindow &window, std::vector<Shape2D> shapes, Drawing d) {
	switch (d) {
		case GRAYSCALE_GRADIENT:
			GrayscaleGradient(window);
			break;
		case TWO_D_GRADIENT:
			TwoDGradient(window);
			break;
		case WITCH_SYMBOL:
			WitchSymbol(window);
			break;
		default: {
			for (Shape2D shape : shapes) {
				DrawShape2D(window, shape);
			}
		}
	}
}

void AddMeshPropertyUpdaters(std::vector<Mesh>& meshes, int frame) {
	for (Mesh& mesh : meshes) {
		if (mesh.name == "tall_box" && frame == 0) {
			/*
			MeshPropertyUpdate smoothnessUpdater(&(mesh.smoothness), 1, 120, true);
			mesh.propertyUpdates.push_back(smoothnessUpdater);
			*/

			/*
			mesh.refractiveIndex = 1.2f;
			MeshPropertyUpdate transparencyUpdater(&(mesh.transparancy), 0, 1, 60, true);
			mesh.propertyUpdates.push_back(transparencyUpdater);
			*/

			MeshVec3PropertyUpdate scalingUpdater(&(mesh.scaleTransition), glm::vec3(1, 1, 1), glm::vec3(0.2, 1, 0.2), 120, false);
			mesh.vec3PropertyUpdates.push_back(scalingUpdater);
			
			MeshVec3PropertyUpdate positionUpdater(&(mesh.translation), glm::vec3(0, 0, 0), glm::vec3(0.8, 0, 0), 120, false);
			mesh.vec3PropertyUpdates.push_back(positionUpdater);
		} else if (mesh.name == "floor" && frame == 0) {
			TextureMap brickWall = TextureMap("Textures/brickwall.ppm");
			TextureMap brickWallNormal = TextureMap("Textures/brickwall_normal.ppm");
			mesh.SetTextureMapper(brickWall, brickWallNormal);
		}
		else if (mesh.name == "short_box") {
			if (frame == 0) {
				MeshVec3PropertyUpdate positionUpdater(&(mesh.translation), glm::vec3(0, 0, 0), glm::vec3(0, 0.4, 0), 120, false);
				mesh.vec3PropertyUpdates.push_back(positionUpdater);
			}
			else if (frame == 60) {
				mesh.refractiveIndex = 2.0f;
				MeshPropertyUpdate transparencyUpdater(&(mesh.transparancy), 0, 1, 60, false);
				mesh.propertyUpdates.push_back(transparencyUpdater);
			}
			else if (frame == 100) {
				MeshPropertyUpdate transparencyUpdater(&(mesh.transparancy), 1, 0, 60, false);
				mesh.propertyUpdates.push_back(transparencyUpdater);
			}
			else if (frame == 130) {
				mesh.refractiveIndex = 1;
			}
		}

	}
}

void run(Drawing draw) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<Shape2D> shapes;
	Camera *camera = nullptr;
	std::vector<Mesh> meshes;
	std::vector<Light> lights = {{{ 0, 0.4f, 2 }, 0.3f }};
	bool FPS_OUTPUT = false;
	bool RECORDING = false;
	bool DYNAMIC_MESH = true;
	int frame = 0;
	std::string directory;

	if (RECORDING) {
		unsigned int ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		char format[] = "Renders/%u/";
		char out[100];
		sprintf(out, format, ms);
		directory = out;
		if (mkdir(out, 0777) == -1) {
			std::cerr << "Error: " << strerror(errno) << std::endl;
		}
	}

	if (draw == POINT_CLOUD || draw == WIRE_FRAME || draw == RASTERISED_3D || draw == RAYTRACE) {
		OBJFile objs("cornell-box.obj", "objs/", 0.35); // 0.35
		meshes = objs.GetMeshes();
		camera = new Camera(glm::vec3(0.0, 0.0, 3.0), glm::vec3(0.0, 0.0, 0.0), 2.0);
	}

	if (DYNAMIC_MESH) {
		AddMeshPropertyUpdaters(meshes, frame);
	}

	while (true) {
		const uint64_t start = SDL_GetPerformanceCounter();
		if (window.pollForInputEvents(event)) handleEvent(event, window, shapes, draw, camera, lights);
		window.clearPixels();
		switch (draw) {
			case RASTERISED_3D:
				DrawRasterized3D(window, *camera, meshes);
				break;
			case RAYTRACE:
				Raytrace(window, *camera, lights, meshes, ShadowType::HARD_SHADOWS, PROXY_AOI_AND_SPECULAR);
				break;
			case POINT_CLOUD:
				shapes = std::vector<Shape2D> { Pointcloud(meshes, *camera) };
				draw2D(window, shapes, draw);
				break;
			case WIRE_FRAME:
				shapes = Wireframe(meshes, *camera);
				draw2D(window, shapes, draw);
				break;
			default:
				draw2D(window, shapes, draw);
		}
		window.renderFrame();
		if (camera && camera->orbiting) {
			camera->orbitStep();
		}
		
		frame++;
		if (DYNAMIC_MESH) {
			for (Mesh& m : meshes) m.UpdateMesh();
			AddMeshPropertyUpdaters(meshes, frame);
		}
		if (RECORDING) {
			std::string filename = std::to_string(frame) + ".bmp";
			window.saveBMP(directory + filename);
		}

		uint64_t end = SDL_GetPerformanceCounter();
		float elapsed = (end - start) / (float)SDL_GetPerformanceFrequency() * 1000.0f;
		SDL_Delay(int(fmax((1.0f / FPS_CAP) * 1000.0f - elapsed, 0)));
		if (FPS_OUTPUT) {
			end = SDL_GetPerformanceCounter();
			elapsed = (end - start) / (float)SDL_GetPerformanceFrequency() * 1000.0f;
			std::cout << "FPS=" <<  1.0f / (elapsed / 1000.0f) << std::endl;
		}

	}
	delete camera;
}

int main(int argc, char *argv[]) {
	run(RAYTRACE);
}
