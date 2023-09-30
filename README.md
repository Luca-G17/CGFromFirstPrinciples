# CGRepo
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
