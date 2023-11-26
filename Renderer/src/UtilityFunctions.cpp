#include "UtilityFunctions.h"

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

double vec3At(const glm::vec3& v, size_t axis) {
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

float CyclicalValue(const double framesPerCycle, const double step) {
	return -(1.0f / 2.0) * cos(2.0 * M_PI * step / framesPerCycle) + 1.0 / 2.0;
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