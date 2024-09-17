# Computer Graphics from First Principles
This repository documents my c++ computer graphics from first principles adventure

## Week:
- [Week 1 - Red Noise](#week-1)
- [Week 2 - Interpolation](#week-2)
- [Week 3 - Basic Geometry](#week-3)
- [Week 4 - 3D Rendering](#week-4)
- [Week 5 - Camera Position/Orientation](#week-5)
- [Week 6 - Raytracing](#week-6)
- [Final Project](#final-project)

## Week 1 ##
Red Noise:  
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/bd99a178-14cb-4f3f-9bc5-5454b463dcca)

## Week 2 ##
Linear Interpolation:  
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/e751f148-f84a-42df-b1b4-aa32e1202ae7)
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/3545f8f5-0b8c-4494-a44f-0ff5d8a0e2c1)

```c++
std::vector<float> interpolateSingleFloats(float from, float to, int steps) {
	float step = (to - from) / (steps - 1); 
	std::vector<float> results;
	for (int i = 0; i < steps; i++) {
		results.push_back(from + (step * i));
	}
	return results;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int steps) {
	glm::vec3 dir = (to - from) * (1.0f / (steps - 1)); 
	std::vector<glm::vec3> results;
	for (int i = 0; i < steps; i++) {
		results.push_back(from + (dir * static_cast<float>(i)));
	}
	return results;
}
```
## Week 3 ##

Basic Geometry:  
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/e0c8f5f0-e567-4f8d-a764-36c9ccff3007)
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/2deaf150-a8b2-4a63-b618-ab05a14fe608)
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/a1c35230-959f-4d59-8966-3c1164cf3efa)


```c++
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
```
## Week 4 ##

3D Rendering:  
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/f955a64e-8ec1-48d8-acb8-472387d67138)
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/362f6bc7-82a8-45a3-9144-fb6759cc8827)
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/2349d127-ff66-4968-9d9c-6b888f4bf221)

```c++
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
	// transform the vertex such that the camera is the origin
	glm::vec3 vPos = vertexPosition - cameraPosition;
	float u = (300 * focalLength * vPos.x / -vPos.z) + WIDTH / 2.0;
	float v = (300 * focalLength * vPos.y / vPos.z) + HEIGHT / 2.0;
	return CanvasPoint(u, v);
}
```
## Week 5 ##

Camera Translation + Rotation:  
![Peek 2023-10-21 18-55](https://github.com/LucaUoB/CGRepo/assets/63655147/4447388c-efba-430d-a13d-6526e7751626)

```c++
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
```

## Week 6 ##

Raytracing:  
![image](https://github.com/LucaUoB/CGRepo/assets/63655147/42318483-8b81-4d19-ba56-9c2c910d3639)  

```c++
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
	RayCollision nearest = { Colour(), glm::vec3(), glm::vec3(), FLT_MAX };
	for (const ModelTriangle& t : triangles) {
		glm::vec3 loc;
		if (TriangleHitLoc(t, r, loc)) {
			float sqrDist = glm::length2(loc - r.origin);
			// Check that the collision is after the start of the ray & collision is closer
			if (glm::dot(r.direction, loc - r.origin) > 0.0 && sqrDist < nearest.distanceToCamera) {
				glm::vec3 triNorm = glm::normalize(glm::cross(t.vertices[0] - t.vertices[1], t.vertices[0] - t.vertices[2]));
				if (glm::dot(triNorm, r.direction) > 0.0)
					triNorm *= -1;
				nearest = { t.colour, loc, triNorm, sqrDist };
			}
		}
	}
	if (nearest.distanceToCamera == FLT_MAX)
		return false;
	collision = nearest;
	return true; 
}

bool InShadow(const glm::vec3& p, const glm::vec3& norm, const std::vector<ModelTriangle>& triangles, const std::vector<glm::vec3> lights) {
	for (const glm::vec3& light : lights) {
		glm::vec3 v = light - p;
		Ray r = { p + norm * 0.001f, v };
		RayCollision c;
		if (!NearestRayCollision(r, triangles, c))
			return false;
	}
	return true;
}
```
## Final Project ##
### Features ###
- [Specular](https://en.wikipedia.org/wiki/Specular_highlight) + [AOI](https://en.wikipedia.org/wiki/Angle_of_incidence_(optics)) + Proximity Lighting
- [Gouraud](https://en.wikipedia.org/wiki/Gouraud_shading) + [Phong Shading](https://en.wikipedia.org/wiki/Phong_shading)
- [Shlick Reflectance](https://en.wikipedia.org/wiki/Schlick%27s_approximation)
- [Dielectric Transmission](https://en.wikipedia.org/wiki/Dielectric_Shader)
- Soft Shadows
- [Normal Mapping Textures](https://en.wikipedia.org/wiki/Normal_mapping)
- Parallelized Ray Tracer

https://github.com/user-attachments/assets/52e87e26-ce96-4f09-9a7e-e80a63dd49e2


