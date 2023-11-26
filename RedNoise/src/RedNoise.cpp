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

#define WIDTH 320
#define HEIGHT 240

struct Shape2D {
	std::vector<CanvasPoint> points;
	Colour colour;
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

Shape2D FillFlatToppedTriangle(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour c) {
	float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
	float slope2 = (v2.x - v0.x) / (v2.y - v0.y);
	float x1 = v0.x;
	float x2 = v0.x;
	std::vector<CanvasPoint> ps;
	for (int y = v0.y; y >= v1.y; y--) {
		Shape2D l = CreateLine2D(CanvasPoint(x1, y), CanvasPoint(x2, y), Colour());
		x1 -= slope1;
		x2 -= slope2;
		ps.insert(ps.begin(), l.points.begin(), l.points.end());
	}
	return { ps, c };
}

Shape2D FillFlatBottomTriangle(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour c) {
	// v2 is the highest point
	float slope1 = (v2.x - v0.x) / (v2.y - v0.y);
	float slope2 = (v2.x - v1.x) / (v2.y - v1.y);
	float x1 = v2.x;
	float x2 = v2.x;
	std::vector<CanvasPoint> ps;
	for (int y = v2.y; y <= v0.y; y++) {
		Shape2D l = CreateLine2D(CanvasPoint(x1, y), CanvasPoint(x2, y), Colour());
		x1 += slope1;
		x2 += slope2;
		ps.insert(ps.begin(), l.points.begin(), l.points.end());
	}
	return { ps, c };
}

void SortVerticies(CanvasTriangle &triangle) {
	std::sort(triangle.vertices.begin(), triangle.vertices.end(), [](const CanvasPoint a, const CanvasPoint b) { return a.y > b.y; } );
}

std::pair<Shape2D, Shape2D> CreateFilledTriangle2D(CanvasTriangle verticies, Colour c) {
	SortVerticies(verticies);
	CanvasPoint v0 = verticies.v0();
	CanvasPoint v1 = verticies.v1();
	CanvasPoint v2 = verticies.v2();
	std::vector<CanvasPoint> ps;
	if (verticies.v1().y == verticies.v2().y) {
		ps = FillFlatBottomTriangle(v0, v1, v2, c).points;
	} 
	else if (v0.y == v1.y) {
		ps = FillFlatToppedTriangle(v0, v1, v2, c).points;
	} 
	else {
		float v3x = v0.x + ((v1.y - v0.y) / (v2.y - v0.y)) * (v2.x - v0.x);
		CanvasPoint v3 = CanvasPoint(v3x, v1.y);
		Shape2D bottom = FillFlatBottomTriangle(v1, v3, v2, c);
		Shape2D top = FillFlatToppedTriangle(v0, v1, v3, c);
		ps.insert(ps.begin(), bottom.points.begin(), bottom.points.end());
		ps.insert(ps.begin(), top.points.begin(), top.points.end());
	}
	Shape2D outline = CreateStrokedTriangle2D(verticies, Colour(255, 255, 255));
	Shape2D shaded = { ps, c };
	return std::make_pair(shaded, outline);
}

TexturePoint CanvasPointToTexturePoint(CanvasPoint p) {
	return TexturePoint(p.x, p.y);
}

glm::vec3 CanvasPointToVec3(CanvasPoint p) {
	return glm::vec3(p.x, p.y, 1.0);
}

std::vector<glm::vec2> InterpolateTwoElementValues(glm::vec2 from, glm::vec2 to, int steps) {
	glm::vec2 dir = (to - from) * (1.0f / (steps - 1));
	std::vector<glm::vec2> results;
	for (int i = 0; i < steps; i++) {
		results.push_back(from + (dir * static_cast<float>(i)));
	}
	return results;
}

Shape2D TexturedFlatToppedTriangle2D(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, CanvasPoint tv0, CanvasPoint tv1, CanvasPoint tv2) {
	// Assume v0 at the bottom
	float slope1 = (v1.x - v0.x) / (v1.y - v0.y);
	float slope2 = (v2.x - v0.x) / (v2.y - v0.y);
	float x1 = v0.x;
	float x2 = v0.x;
	std::vector<CanvasPoint> ps;

	glm::vec2 diff1 = CanvasPointToVec2(v1) - CanvasPointToVec2(v0);
	float steps1 = fmax(fabs(diff1.x), fabs(diff1.y));

	glm::vec2 diff2 = CanvasPointToVec2(v2) - CanvasPointToVec2(v0);
	float steps2 = fmax(fabs(diff2.x), fabs(diff2.y));
	float steps = fmin(steps1, steps2);
	std::vector<glm::vec2> tl1 = InterpolateTwoElementValues(CanvasPointToVec2(tv0), CanvasPointToVec2(tv1), steps);
	std::vector<glm::vec2> tl2 = InterpolateTwoElementValues(CanvasPointToVec2(tv0), CanvasPointToVec2(tv2), steps);
	for (int s = 0; s < steps; s++) {
		int y = v0.y - s;
		Shape2D l = CreateLine2D(CanvasPoint(x1, y), CanvasPoint(x2, y), Colour());
		x1 -= slope1;
		x2 -= slope2;
		std::vector<glm::vec2> textureScanline = InterpolateTwoElementValues(tl1[s], tl2[s], l.points.size());
		for (int i = 0; i < l.points.size(); i++) {
			l.points[i].texturePoint = TexturePoint(textureScanline[i].x, textureScanline[i].y);
		}
		ps.insert(ps.begin(), l.points.begin(), l.points.end());
	}
	return { ps, Colour() };
}

Shape2D TexturedFlatBottomTriangle2D(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, CanvasPoint tv0, CanvasPoint tv1, CanvasPoint tv2) {
	float slope1 = (v2.x - v0.x) / (v2.y - v0.y);
	float slope2 = (v2.x - v1.x) / (v2.y - v1.y);
	float x1 = v2.x;
	float x2 = v2.x;
	std::vector<CanvasPoint> ps;

	glm::vec2 diff1 = CanvasPointToVec2(v0) - CanvasPointToVec2(v2);
	float steps1 = fmax(fabs(diff1.x), fabs(diff1.y));

	glm::vec2 diff2 = CanvasPointToVec2(v1) - CanvasPointToVec2(v2);
	float steps2 = fmax(fabs(diff2.x), fabs(diff2.y));
	float steps = fmin(steps1, steps2);
	std::vector<glm::vec2> tl1 = InterpolateTwoElementValues(CanvasPointToVec2(tv2), CanvasPointToVec2(tv0), steps);
	std::vector<glm::vec2> tl2 = InterpolateTwoElementValues(CanvasPointToVec2(tv2), CanvasPointToVec2(tv1), steps);
	for (int s = 0; s < steps; s++) {
		int y = v2.y + s;
		Shape2D l = CreateLine2D(CanvasPoint(x1, y), CanvasPoint(x2, y), Colour());
		x1 += slope1;
		x2 += slope2;
		std::vector<glm::vec2> textureScanline = InterpolateTwoElementValues(tl1[s], tl2[s], l.points.size());
		for (int i = 0; i < l.points.size(); i++) {
			l.points[i].texturePoint = TexturePoint(textureScanline[i].x, textureScanline[i].y);
		}
		ps.insert(ps.begin(), l.points.begin(), l.points.end());
	}
	return { ps, Colour() };
}

CanvasPoint TexturePointToCanvasPoint(TexturePoint t) {
	return CanvasPoint(t.x, t.y);
} 

std::pair<Shape2D, Shape2D> CreateTexturedTriangle2D(CanvasTriangle verticies, CanvasTriangle texturedTriangle) {
	verticies.v0().texturePoint = CanvasPointToTexturePoint(texturedTriangle.v0());
	verticies.v1().texturePoint = CanvasPointToTexturePoint(texturedTriangle.v1());
	verticies.v2().texturePoint = CanvasPointToTexturePoint(texturedTriangle.v2());

	SortVerticies(verticies);
	CanvasPoint v0 = verticies.v0();
	CanvasPoint v1 = verticies.v1();
	CanvasPoint v2 = verticies.v2();
	CanvasPoint vt0 = TexturePointToCanvasPoint(v0.texturePoint);
	CanvasPoint vt1 = TexturePointToCanvasPoint(v1.texturePoint);
	CanvasPoint vt2 = TexturePointToCanvasPoint(v2.texturePoint);
	std::vector<CanvasPoint> ps;
	if (verticies.v1().y == verticies.v2().y) {
		ps = TexturedFlatBottomTriangle2D(v0, v1, v2, vt0, vt1, vt2).points;
	} 
	else if (v0.y == v1.y) {
		ps = TexturedFlatToppedTriangle2D(v0, v1, v2, vt0, vt1, vt2).points;
	} 
	else {
		float v3x = v0.x + ((v1.y - v0.y) / (v2.y - v0.y)) * (v2.x - v0.x);
		CanvasPoint v3 = CanvasPoint(v3x, v1.y);
		// Get v3 in texture coordinates.
		glm::vec2 v2v0 = CanvasPointToVec2(v0) - CanvasPointToVec2(v2);
		glm::vec2 v2v3 = CanvasPointToVec2(v3) - CanvasPointToVec2(v2);
		glm::vec2 vt2vt0 = CanvasPointToVec2(vt0) - CanvasPointToVec2(vt2);  
		float s = glm::length(v2v0) / glm::length(v2v3);
		glm::vec2 vt3vec2 = vt2vt0 * s;
		CanvasPoint vt3 = CanvasPoint(vt3vec2.x, vt3vec2.y);

		Shape2D bottom = TexturedFlatBottomTriangle2D(v1, v3, v2, vt1, vt3, vt2);
		Shape2D top = TexturedFlatToppedTriangle2D(v0, v1, v3, vt0, vt1, vt3);
		ps.insert(ps.begin(), bottom.points.begin(), bottom.points.end());
		ps.insert(ps.begin(), top.points.begin(), top.points.end());
	}
	Shape2D outline = CreateStrokedTriangle2D(verticies, Colour(255, 255, 255));
	Shape2D shaded = { ps, Colour() };
	return std::make_pair(shaded, outline);
}

Shape2D CreateTexturedTriangle2DBarycentric(CanvasTriangle verticies, CanvasTriangle texturedTriangle) {
	verticies.v0().texturePoint = CanvasPointToTexturePoint(texturedTriangle.v0());
	verticies.v1().texturePoint = CanvasPointToTexturePoint(texturedTriangle.v1());
	verticies.v2().texturePoint = CanvasPointToTexturePoint(texturedTriangle.v2());
	CanvasPoint v0 = verticies.v0();
	CanvasPoint v1 = verticies.v1();
	CanvasPoint v2 = verticies.v2();
	int maxX = fmax(v0.x, fmax(v1.x, v2.x));
	int minX = fmin(v0.x, fmin(v1.x, v2.x));
	int maxY = fmax(v0.y, fmax(v1.y, v2.y));
	int minY = fmin(v0.y, fmin(v1.y, v2.y));
	glm::vec3 A = CanvasPointToVec3(v0);
	glm::vec3 B = CanvasPointToVec3(v1);
	glm::vec3 C = CanvasPointToVec3(v2);
	float area = glm::length(glm::cross((B - A), (C - A)));

	std::vector<CanvasPoint> ps;
	for (int y = minY; y <= maxY; y++) {
		for (int x = minX; x <= maxX; x++) {
			glm::vec3 P = CanvasPointToVec3(CanvasPoint(x, y));
			glm::vec3 PC = C - P;
			glm::vec3 PB = B - P;
			glm::vec3 PA = A - P;
			float alpha = glm::length(glm::cross(PB, PC)) / area;
			float beta = glm::length(glm::cross(PC, PA)) / area;
			float gamma = glm::length(glm::cross(PA, PB)) / area;
			if ((alpha >= 0.0 && alpha <= 1.0) && (beta >= 0.0 && beta <= 1.0) && (gamma >= 0.0 && gamma <= 1.0) && fabs(alpha + beta + gamma - 1) <= 0.01) {
				glm::vec2 textureCoord = 
					CanvasPointToVec2(texturedTriangle.v0()) * alpha
					+ CanvasPointToVec2(texturedTriangle.v1()) * beta
					+ CanvasPointToVec2(texturedTriangle.v2()) * gamma;
				CanvasPoint p = CanvasPoint(x, y);
				if (x == 300 && y == 230) {
					std::cout << textureCoord.x << ", " << textureCoord.y << std::endl;
				}
				p.texturePoint = TexturePoint(textureCoord.x, textureCoord.y);
				ps.push_back(p);
			}
		}
	}
	return { ps, Colour() };
}

std::vector<float> InterpolateSingleFloats(float from, float to, int steps) {
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

void RedNoise(DrawingWindow &window) {
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
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
		std::pair<Shape2D, Shape2D> filled = CreateFilledTriangle2D(
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
	int col = PackColour(shape.colour.red, shape.colour.blue, shape.colour.green);
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
I’ve found that using the API gives you a lot more control and therefore less nonsense, but I’m also not trying to be an edgelord with it, so I’m not pushing it to go to weird places so YMMV.
void TestTextureMapping(DrawingWindow &window, TextureMap texture) {
	CanvasTriangle textureTriangle = CanvasTriangle(CanvasPoint(195, 10), CanvasPoint(395, 380), CanvasPoint(65, 330));
	CanvasTriangle triangle = CanvasTriangle(CanvasPoint(160, 10), CanvasPoint(300, 230), CanvasPoint(10, 150));
	std::pair<Shape2D, Shape2D> tri = CreateTexturedTriangle2D(triangle, textureTriangle);
	DrawTexturedShape2D(window, tri.first, texture);
	DrawShape2D(window, tri.second);
}

void draw(DrawingWindow &window, std::vector<Shape2D> shapes, TextureMap texture) {
	window.clearPixels();
	// RedNoise(window);
	// GrayscaleGradient(window);
	// TwoDGradient(window);
	WitchSymbol(window);

	// TestTextureMapping(window, texture);
	for (Shape2D shape : shapes) {
		DrawShape2D(window, shape);
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window, std::vector<Shape2D> &shapes) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) AddRandomTriangle(false, shapes);
		else if (event.key.keysym.sym == SDLK_f) AddRandomTriangle(true, shapes);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<Shape2D> shapes;
	std::string filename = "texture.ppm";
	TextureMap texture = TextureMap(filename);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, shapes);
		draw(window, shapes, texture);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
