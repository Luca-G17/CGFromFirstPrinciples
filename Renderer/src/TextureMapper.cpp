#include "TextureMapper.h"

// TODO: this is potentially stupid
std::vector<std::pair<size_t, glm::vec3>> GetRectangleCorners(const std::vector<glm::vec3>& vertices) {
	std::vector<std::pair<size_t, glm::vec3>> corners;
	for (size_t i = 0; i < vertices.size(); i++) {
		corners.push_back(std::make_pair(i, vertices[i]));
	}
	return corners;
}

std::vector<std::vector<Colour>> TextureMapper::TextureMapToColourArray(TextureMap& textureMap) {
    std::vector<std::vector<Colour>> colourArray;
    for (size_t y = 0; y < textureMap.height; y++) {
        colourArray.push_back({});
        for (size_t x = 0; x < textureMap.width; x++) {
            colourArray[y].push_back(Uint32ToColour(textureMap.pixels[y * textureMap.width + x]));
        }
    }
    return colourArray;
}

std::vector<std::vector<glm::vec3>> TextureMapper::NormalMapToVec3Array(const TextureMap& normalMap, const glm::vec3 normal) {
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

TextureMapper::TextureMapper(TextureMap texture, TextureMap normal, glm::vec3 faceNormal, TextureMappingPrimative primative) {
    this->primative = primative;
    this->textureMap = TextureMapToColourArray(texture);
    this->normalMap = NormalMapToVec3Array(normal, faceNormal);
    this->type = TextureMappingType::NORMAL;
}


TextureMapper::TextureMapper(TextureMap texture, TextureMappingPrimative primative) {
    this->primative = primative;
    this->textureMap = TextureMapToColourArray(texture);
    this->type = TextureMappingType::FLAT;
}

TextureMapper::TextureMapper(){
}

TextureVertex TextureMapper::GetVertexMapping(size_t v) const {
    for (const TextureVertex& t : vertexMapping) {
        if (t.vertex == v)
            return t;
    }
    return {};
}

glm::vec2 TextureMapper::BarycentricToMap(const Triangle& t, const glm::vec3 barycentric) const {
    return GetVertexMapping(t.vertices[0]).textureCoord * barycentric[0] +  
            GetVertexMapping(t.vertices[1]).textureCoord * barycentric[1] +
            GetVertexMapping(t.vertices[2]).textureCoord * barycentric[2];
}

Colour TextureMapper::BarycentricTexture(const Triangle& t, const glm::vec3 barycentric) const {
    glm::vec2 c = BarycentricToMap(t, barycentric);
    return textureMap[(int)c.y][(int)c.x];
}

glm::vec3 TextureMapper::BarycentricNormal(const Triangle& t, const glm::vec3 barycentric) const {
    glm::vec2 c = BarycentricToMap(t, barycentric);
    return normalMap[(int)c.y][(int)c.x];
}

void TextureMapper::VertexMapping(const std::vector<glm::vec3>& vertices, std::vector<Triangle> triangles) {
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
            vertexMapping.push_back({ bottomLeft, {0, h - 1} });
            vertexMapping.push_back({ bottomRight, {w - 1, h - 1} });
            vertexMapping.push_back({ topLeft, {0, 0} });
            vertexMapping.push_back({ topRight, {w - 1, 0} });
            break;
        }
        default: {
            break;
        }
    }
}