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

#define WIDTH 680
#define HEIGHT 520

struct Shape2D {
	std::vector<CanvasPoint> points;
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

Shape2D CreateFilledTriangle2D(CanvasTriangle verticies, Colour c, ModelTriangle* pModelTriangle, std::vector<std::vector<DepthPoint>> *pDepthBuffer) {
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

std::vector<float> InterpolateSingleFloats(float from, float to, int steps) {
	if (steps <= 1) {
		return { from };
	}
	float step = (to - from) / (steps - 1); 
	std::vector<float> results;
	for (int i = 0; i < steps; i++) {
		results.push_back(from + (step * i));
	}
	return results;
}

std::vector<glm::vec3> InterpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int steps) {
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

uint32_t PackColour(uint8_t r, uint8_t g, uint8_t b) {
	return (255 << 24) + (r << 16) + (g << 8) + b;
}

void GrayscaleGradient(DrawingWindow &window) {
	std::vector<float> row = InterpolateSingleFloats(255, 0, window.width);
	std::vector<uint32_t> rowColours;
	std::transform(row.cbegin(), row.cend(), std::back_inserter(rowColours), [](float f) -> uint32_t { return PackColour(f, f, f); });
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			window.setPixelColour(x, y, rowColours.at(x));
		}
	}
}

void TwoDGradient(DrawingWindow &window) {
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	std::vector<glm::vec3> left = InterpolateThreeElementValues(topLeft, bottomLeft, window.height);
	std::vector<glm::vec3> right = InterpolateThreeElementValues(topRight, bottomRight, window.height);
	
	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> row = InterpolateThreeElementValues(left[y], right[y], window.width);
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

void AddRandomTriangle(bool filled, std::vector<Shape2D> &shapes) {
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

void DrawPoints(DrawingWindow &window, std::vector<glm::vec2> points, int colour) {
	for (glm::vec2 p : points) {
		window.setPixelColour(p.x, p.y, colour);
	}
}

void DrawTexturedShape2D(DrawingWindow &window, Shape2D shape, TextureMap texture) {
	for (CanvasPoint p : shape.points) {
		int pixelIndex = p.texturePoint.y * texture.width + p.texturePoint.x;
		if (pixelIndex < texture.pixels.size()) {
			uint32_t col = texture.pixels.at(pixelIndex);
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

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
	// transform the vertex such that the camera is the origin
	glm::vec3 vPos = vertexPosition - cameraPosition;
	float u = (300 * focalLength * vPos.x / -vPos.z) + WIDTH / 2.0;
	float v = (300 * focalLength * vPos.y / vPos.z) + HEIGHT / 2.0;
	return CanvasPoint(u, v);
}

struct Camera {
	glm::vec3 position;
	float focalLength;
};

Shape2D Pointcloud(std::vector<ModelTriangle> triangles, Camera camera) {
	std::vector<CanvasPoint> ps;
	for (ModelTriangle triangle : triangles) {
		for (glm::vec3 vertex : triangle.vertices) {
			ps.push_back(getCanvasIntersectionPoint(camera.position, vertex, camera.focalLength));
		}
	}
	return { ps, Colour(255, 255, 255) };
}

std::vector<Shape2D> Wireframe(std::vector<ModelTriangle> triangles, Camera camera) {
	std::vector<Shape2D> wireframe;
	for (ModelTriangle triangle : triangles) {
		CanvasPoint p1 = getCanvasIntersectionPoint(camera.position, triangle.vertices[0], camera.focalLength);
		CanvasPoint p2 = getCanvasIntersectionPoint(camera.position, triangle.vertices[1], camera.focalLength);
		CanvasPoint p3 = getCanvasIntersectionPoint(camera.position, triangle.vertices[2], camera.focalLength);

		Shape2D triangle2D = CreateStrokedTriangle2D(CanvasTriangle(p1, p2, p3), Colour(255, 255, 255));
		wireframe.push_back(triangle2D);
	}
	return wireframe;
}

ModelTriangle TranslateTriangle(ModelTriangle tri, glm::vec3 t) {
	return ModelTriangle(tri.vertices[0] + t, tri.vertices[1] + t, tri.vertices[2] + t, tri.colour);
}

std::vector<std::vector<DepthPoint>> RasterisedRender(DrawingWindow& window, std::vector<ModelTriangle> triangles, Camera camera) {
	std::vector<std::vector<DepthPoint>> depthBuffer;
	for (int y = 0; y < HEIGHT; y++) {
		depthBuffer.push_back(std::vector<DepthPoint>());
		for (int x = 0; x < WIDTH; x++)
			depthBuffer[y].push_back({ Colour(), 0 });
	}
	for (ModelTriangle triangle : triangles) {
		// We must maintain a depth buffer as we build the triangles
		CanvasPoint p1 = getCanvasIntersectionPoint(camera.position, triangle.vertices[0], camera.focalLength);
		CanvasPoint p2 = getCanvasIntersectionPoint(camera.position, triangle.vertices[1], camera.focalLength);
		CanvasPoint p3 = getCanvasIntersectionPoint(camera.position, triangle.vertices[2], camera.focalLength);
		ModelTriangle cameraSpaceTriangle = TranslateTriangle(triangle, -camera.position);
		CreateFilledTriangle2D(CanvasTriangle(p1, p2, p3), triangle.colour, &cameraSpaceTriangle, &depthBuffer);
	}
	return depthBuffer;
}

enum Drawing {
	GRAYSCALE_GRADIENT,
	TWO_D_GRADIENT,
	WITCH_SYMBOL,
	RANDOM_TRIANGLES,
	POINT_CLOUD,
	WIRE_FRAME,
	RASTERISED_3D,
};

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

bool handleEvent(SDL_Event event, DrawingWindow &window, std::vector<Shape2D> &shapes, Drawing d) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u && d == RANDOM_TRIANGLES) AddRandomTriangle(false, shapes);
		else if (event.key.keysym.sym == SDLK_f && d == RANDOM_TRIANGLES) AddRandomTriangle(true, shapes);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	} else if (event.type == SDL_QUIT) {
		return true;
	}
	return false;
}

void run(Drawing draw) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<Shape2D> shapes;
	bool quit = false;

	if (draw == POINT_CLOUD || draw == WIRE_FRAME || draw == RASTERISED_3D) {
		OBJFile objs("cornell-box.obj", "objs/", 0.35);
		std::vector<ModelTriangle> triangles = objs.GetTriangles();
		Camera camera = { glm::vec3(0.0, 0.0, 4.0), 2.0 };
		if (draw == RASTERISED_3D) {
			std::vector<std::vector<DepthPoint>> depthBuffer = RasterisedRender(window, triangles, camera);
			while (!quit) {
				if (window.pollForInputEvents(event)) quit = handleEvent(event, window, shapes, draw);
				window.clearPixels();
				for (int y = 0; y < HEIGHT; y++) {
					for (int x = 0; x < WIDTH; x++) {
						Colour c  = depthBuffer[y][x].colour;
						window.setPixelColour(x, y, PackColour(c.red, c.green, c.blue));
					}
				}				
				window.renderFrame();
			}
			return;
		}
		else {
			if (draw == POINT_CLOUD) shapes.push_back(Pointcloud(triangles, camera));
			else shapes = Wireframe(triangles, camera);
		}
	}

	while (!quit) {
		if (window.pollForInputEvents(event)) quit = handleEvent(event, window, shapes, draw);
		draw2D(window, shapes, draw);
		window.renderFrame();
	}
}

int main(int argc, char *argv[]) {
	const rlim_t stackSize = 16 * 1024 * 1024;
	struct rlimit r1;
	r1.rlim_cur = stackSize;
	if (setrlimit(RLIMIT_STACK, &r1) != 0)
		std::cerr << "setrlimit returned error\n";
	
	run(WIRE_FRAME);
}
