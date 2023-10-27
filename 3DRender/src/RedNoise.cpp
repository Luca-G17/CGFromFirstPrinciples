#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include <cmath>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <ModelTriangle.h>
#include <sys/resource.h>
#include <math.h>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#include <limits.h>

#define WIDTH 440
#define HEIGHT 300

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

struct OBJ {
	std::vector<ModelTriangle> triangles;
	std::string name;
	OBJ(std::vector<glm::vec3> coords, std::vector<std::vector<int>> faces, std::string name, Colour material) : 
			name(name) {
		for (std::vector<int> face : faces) {
			triangles.push_back(ModelTriangle(coords[face[0] - 1], coords[face[1] - 1], coords[face[2] - 1], material));
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
			while (getline(File, line)) {
				std::vector<std::string> words = Split(line, " ");
				if (words[0] == "o") { // start parsing new object
					std::string name = words[1];

					getline(File, line);
					std::string material = Split(line, " ")[1];
					getline(File, line);
					words = Split(line, " ");
					while (words[0] == "v" ) {
						glm::vec3 coord = { std::stof(words[1]), std::stof(words[2]), std::stof(words[3]) };
						coords.push_back(coord * scaling);
						getline(File, line);
						words = Split(line, " ");
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
					objects.push_back(OBJ(coords, faces, name, materials[material]));
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

	std::vector<ModelTriangle> GetTriangles() {
		std::vector<ModelTriangle> triangles;
		for (OBJ obj : objects) {
			for (ModelTriangle triangle : obj.triangles)
				triangles.push_back(triangle);
		}
		return triangles;
	}
};

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
		std::cout << glm::to_string(rotation) << std::endl;
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
	const float y = (v - (HEIGHT / 2.0)) / 150.0;
	return camera.rotationMatrix * glm::vec3(x, y, camera.focalLength);
}


Shape2D Pointcloud(std::vector<ModelTriangle> triangles, Camera camera) {
	std::vector<CanvasPoint> ps;
	for (ModelTriangle triangle : triangles) {
		for (glm::vec3 vertex : triangle.vertices) {
			ps.push_back(getCanvasIntersectionPoint(camera, vertex));
		}
	}
	return { ps, Colour(255, 255, 255) };
}

std::vector<Shape2D> Wireframe(std::vector<ModelTriangle> triangles, Camera camera) {
	std::vector<Shape2D> wireframe;
	for (ModelTriangle triangle : triangles) {
		CanvasPoint p1 = getCanvasIntersectionPoint(camera, triangle.vertices[0]);
		CanvasPoint p2 = getCanvasIntersectionPoint(camera, triangle.vertices[1]);
		CanvasPoint p3 = getCanvasIntersectionPoint(camera, triangle.vertices[2]);
		Shape2D triangle2D = CreateStrokedTriangle2D(CanvasTriangle(p1, p2, p3), Colour(255, 255, 255));
		wireframe.push_back(triangle2D);
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

std::vector<std::vector<DepthPoint>> RasterisedRender(DrawingWindow& window, const std::vector<ModelTriangle> triangles, const Camera camera) {
	std::vector<std::vector<DepthPoint>> depthBuffer;
	for (int y = 0; y < HEIGHT; y++) {
		depthBuffer.push_back(std::vector<DepthPoint>());
		for (int x = 0; x < WIDTH; x++)
			depthBuffer[y].push_back({ Colour(), 0 });
	}
	for (ModelTriangle triangle : triangles) {
		// We must maintain a depth buffer as we build the triangles
		const CanvasPoint p1 = getCanvasIntersectionPoint(camera, triangle.vertices[0]);
		const CanvasPoint p2 = getCanvasIntersectionPoint(camera, triangle.vertices[1]);
		const CanvasPoint p3 = getCanvasIntersectionPoint(camera, triangle.vertices[2]);
		const ModelTriangle cameraSpaceTriangle = RotateTriangle(TranslateTriangle(triangle, -camera.position), camera.rotationMatrix);
		if (PointInsideCanvas(p1) || PointInsideCanvas(p2) || PointInsideCanvas(p3)) {
			CreateFilledTriangle2D(CanvasTriangle(p1, p2, p3), triangle.colour, &cameraSpaceTriangle, &depthBuffer);
		}
	}
	return depthBuffer;
}

struct RayCollision {
	Colour color;
	glm::vec3 position;
	float distanceToCamera;
};

struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;
};

bool TriangleHitLoc(const ModelTriangle& tri, const Ray& r, glm::vec3& loc) {
	// Get ray intersection with the plane described by two verticies of the triangle
	glm::vec3 normal = glm::cross(tri.vertices[1] - tri.vertices[0], tri.vertices[2] - tri.vertices[0]);
	float d = glm::dot(tri.vertices[0], normal);
	float n = glm::dot(normal, r.origin); // n.p_0
	float m = glm::dot(normal, r.direction); // n.u
	if (m == 0) 
		return false;
	float t = (d - n) / m; // tri = (d - n.p0) / n.u
	loc = r.origin + (r.direction * t);

	// Check that loc is inside triangle
	float area2 = glm::length(normal);
	glm::vec3 PC = tri.vertices[2] - loc;
	glm::vec3 PB = tri.vertices[1] - loc;
	glm::vec3 PA = tri.vertices[0] - loc;
	float alpha = glm::length(glm::cross(PB, PC)) / area2;
	float beta = glm::length(glm::cross(PC, PA)) / area2;
	float gamma = glm::length(glm::cross(PA, PB)) / area2;	
	return (alpha >= 0 && alpha <= 1) &&
		   (beta  >= 0 && beta  <= 1) &&
		   (gamma >= 0 && gamma <= 1) &&
		   (fabs(alpha + beta + gamma - 1) <= 0.01);
}

// Returning a bool here cause someone decided we are using C++11 and therefore I can't use std::optional :(
bool NearestRayCollision(Ray& r, const std::vector<ModelTriangle>& triangles, RayCollision& collision) {
	RayCollision nearest = { Colour(), glm::vec3(), FLT_MAX };
	for (ModelTriangle t : triangles) {
		glm::vec3 loc;
		if (TriangleHitLoc(t, r, loc)) {
			float sqrDist = glm::length2(loc - r.origin);
			// Check that the collision is after the start of the ray & collision is closer
			if (glm::dot(r.direction, r.origin) > 0.0 && sqrDist < nearest.distanceToCamera) {
				nearest = { t.colour, loc, sqrDist };
			}
		}
	}
	if (nearest.distanceToCamera == FLT_MAX)
		return false;
	collision = nearest;
	return true; 
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

bool handleEvent(const SDL_Event event, DrawingWindow &window, std::vector<Shape2D> &shapes, const Drawing d, Camera* const pCamera) {
	if (event.type == SDL_KEYDOWN) {
	 	if (event.key.keysym.sym == SDLK_u && d == RANDOM_TRIANGLES) AddRandomTriangle(false, shapes);
		else if (event.key.keysym.sym == SDLK_f && d == RANDOM_TRIANGLES) AddRandomTriangle(true, shapes);
		else if (pCamera) {
			Camera& c = *pCamera;
			constexpr float step = 0.1;
			constexpr float rotStep = M_PI / 180;
			// Translations
			if (event.key.keysym.sym == SDLK_w && d == RASTERISED_3D) c.AddTranslation(0, 0, -step);
			else if (event.key.keysym.sym == SDLK_s && d == RASTERISED_3D) c.AddTranslation(0, 0, step);
			else if (event.key.keysym.sym == SDLK_a && d == RASTERISED_3D) c.AddTranslation(-step, 0, 0);
			else if (event.key.keysym.sym == SDLK_d && d == RASTERISED_3D) c.AddTranslation(step, 0, 0);
			else if (event.key.keysym.sym == SDLK_SPACE && d == RASTERISED_3D) c.AddTranslation(0, step, 0);
			else if (event.key.keysym.sym == SDLK_c && d == RASTERISED_3D) c.AddTranslation(0, -step, 0);

			// Rotations
			else if (event.key.keysym.sym == SDLK_LEFT) c.AddRotation(0, -rotStep, 0);
			else if (event.key.keysym.sym == SDLK_RIGHT) c.AddRotation(0, rotStep, 0);
			else if (event.key.keysym.sym == SDLK_UP) c.AddRotation(0, 0, rotStep);
			else if (event.key.keysym.sym == SDLK_DOWN) c.AddRotation(0, 0, -rotStep);

			else if (event.key.keysym.sym == SDLK_l) c.lookAt(glm::vec3(0, 0, 0));
			else if (event.key.keysym.sym == SDLK_o) {
				if (c.orbiting) {
					c.stopOrbit();
				} else {
					c.startOrbit(glm::vec3(0, 0, 0), 4);
				}
			}
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	} else if (event.type == SDL_QUIT) {
		return true;
	}
	return false;
}

void Raytrace(DrawingWindow& window, const Camera& camera, const std::vector<ModelTriangle>& triangles) {
	for (int v = 0; v < HEIGHT; v++) {
		for (int u = 0; u < WIDTH; u++) {
			// first find coords in 3D
			glm::vec3 rayV = CanvasPointToWorld(camera, u, v);
			glm::vec3 rayP = camera.position + rayV;
			Ray ray = { rayP, rayV };
			RayCollision c;
			uint32_t colour;
			if (NearestRayCollision(ray, triangles, c))
				colour = PackColour(c.color.red, c.color.green, c.color.green);
			else
				colour = PackColour(0, 0, 0);
			window.setPixelColour(u, v, colour);
		}
	}
}

void DrawRasterized3D(DrawingWindow& window, Camera& camera, std::vector<ModelTriangle>& triangles) {
	bool quit = false;
	SDL_Event event;
	std::vector<Shape2D> shapes;
	while (!quit) {
		if (window.pollForInputEvents(event)) quit = handleEvent(event, window, shapes, RASTERISED_3D, &camera);
		window.clearPixels();
		std::vector<std::vector<DepthPoint>> depthBuffer = RasterisedRender(window, triangles, camera);
		for (int y = 0; y < HEIGHT; y++) {
			for (int x = 0; x < WIDTH; x++) {
				Colour c  = depthBuffer[y][x].colour;
				window.setPixelColour(x, y, PackColour(c.red, c.green, c.blue));
			}
		}
		if (camera.orbiting) {
			camera.orbitStep();
		}
		// camera.output();
		window.renderFrame();
	}
	return;
}

void draw2D(DrawingWindow &window, std::vector<Shape2D> shapes, Drawing d) {
	window.clearPixels();
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

void run(Drawing draw) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<Shape2D> shapes;
	bool quit = false;

	if (draw == POINT_CLOUD || draw == WIRE_FRAME || draw == RASTERISED_3D || draw == RAYTRACE) {
		OBJFile objs("cornell-box.obj", "objs/", 0.35);
		std::vector<ModelTriangle> triangles = objs.GetTriangles();
		Camera camera(glm::vec3(0.0, 0.0, 4.0), glm::vec3(0.0, 0.0, 0.0), 2.0);
		switch (draw) {
			case RASTERISED_3D:
				DrawRasterized3D(window, camera, triangles);
				break;
			case RAYTRACE:
				Raytrace(window, camera, triangles);
				break;
			case POINT_CLOUD:
				shapes.push_back(Pointcloud(triangles, camera));
				break;
			case WIRE_FRAME:
				shapes = Wireframe(triangles, camera);
				break;
		}
	}

	while (!quit) {
		if (window.pollForInputEvents(event)) quit = handleEvent(event, window, shapes, draw, nullptr);
		draw2D(window, shapes, draw);
		window.renderFrame();
	}
}

int main(int argc, char *argv[]) {
	constexpr rlim_t stackSize = 16 * 1024 * 1024;
	struct rlimit r1;
	r1.rlim_cur = stackSize;
	if (setrlimit(RLIMIT_STACK, &r1) != 0)
		std::cerr << "setrlimit returned error\n";
	
	run(RAYTRACE);
}
