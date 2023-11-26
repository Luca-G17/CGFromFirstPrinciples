#include "Mesh.h"

// --------------------------------------------------------------------------------------------
// MESH PROPERTY UPDATE
// --------------------------------------------------------------------------------------------

MeshPropertyUpdate::MeshPropertyUpdate(double* value, double start, double finish, double speed, bool cyclical) {
    this->value = value;
    this->start = start;
    this->finish = finish;
    this->speed = speed;
    this->step = 0;
    this->cyclical = cyclical;
}

MeshPropertyUpdate::MeshPropertyUpdate(double* value, double finish, double speed, bool cyclical) {
    this->value = value;
    this->start = *value;
    this->finish = finish;
    this->speed = speed;
    this->step = 0;
    this->cyclical = cyclical;
}

bool MeshPropertyUpdate::Update() {
    step++;
    double v = CyclicalValue(speed, step);
    *value = start + (finish - start) * v; 
    return NotFinished();
}

bool MeshPropertyUpdate::NotFinished() {
    if (cyclical)
        return true;
    
    return ((finish - *value) * (finish - start) > 0);
}

// --------------------------------------------------------------------------------------------
// MESH VEC3 PROPERTY UPDATE
// --------------------------------------------------------------------------------------------

MeshVec3PropertyUpdate::MeshVec3PropertyUpdate(glm::vec3* value, glm::vec3 start, glm::vec3 finish, double speed, bool cyclical) {
    this->value = value;
    this->start = start;
    this->finish = finish;
    this->speed = speed;
    this->step = 0;
    this->cyclical = cyclical;
}

MeshVec3PropertyUpdate::MeshVec3PropertyUpdate(glm::vec3* value, glm::vec3 start, glm::vec3 bezierMid ,glm::vec3 finish, double speed, bool cyclical) {
    this->value = value;
    this->start = start;
    this->finish = finish;
    this->speed = speed;
    this->step = 0;
    this->cyclical = cyclical;
    this->bezierMid = bezierMid;
    this->bezier = true;
}

bool MeshVec3PropertyUpdate::Update() {
    step++;
    float v = CyclicalValue(speed, step);
    if (bezier) {
        *value = ((1 - v) * (((1 - v) * start) + (v * bezierMid))) + (v * (((1 - v) * bezierMid) + (v * finish)));
    }
    else {
        *value = start + (finish - start) * v;
    }
    return NotFinished(v);
}

bool MeshVec3PropertyUpdate::NotFinished(double v) {
    if (cyclical)
        return true;
    
    return ((1 - v) > 0.0001f); 
}

// --------------------------------------------------------------------------------------------
// MESH
// --------------------------------------------------------------------------------------------

void Mesh::SetTextureMapper(TextureMap textureMap, TextureMap normalMap) {
    this->texture = TextureMapper(textureMap, normalMap, this->triangles[0].normal, TextureMappingPrimative::RECTANGLE);
    this->texture.VertexMapping(vertices, triangles);
}

void Mesh::SetTextureMapper(TextureMap map) {
    this->texture = TextureMapper(map, TextureMappingPrimative::RECTANGLE);
    this->texture.VertexMapping(vertices, triangles);
}

const glm::vec3& Mesh::GetVertex(const size_t t, const size_t v) const {
    return vertices[triangles[t].vertices[v]];
}

void Mesh::SetVertex(const size_t t, const size_t v, const glm::vec3 vertex) {
    vertices[triangles[t].vertices[v]] = vertex;
}

// Vertex and face normals will require recomputation if I add rotations
const glm::vec3& Mesh::GetVertexNormal(size_t t, size_t v) const {
    return vertexNormals[triangles[t].vertices[v]];
}

const Colour Mesh::GetColour() const {
    return colour;
}

void Mesh::UpdateMesh() {
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

    // Position = Previous Position, Translation = New Position after translation
    glm::vec3 movement = translation - position;
    if (glm::dot(movement, movement) > 0.000001f) {
        Translate(movement);
        position += movement;
    }

    if (fabs(glm::dot(scaleTransition, scaleTransition) - 3) > 0.000001f) {
        Scale(InvertVec3(currentScale));
        Scale(scaleTransition);
        currentScale = scaleTransition;
    }
    MeshBoundingBox();
}

std::pair<glm::vec3, glm::vec3> Mesh::MeshBoundingBox() {
    glm::vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
    glm::vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (Triangle& t : triangles) {
        for (size_t axis = 0; axis < 3; axis++) {
            vec3SetAt(min, axis, fmin(vec3At(min, axis), vec3At(t.bbox.first, axis)));
            vec3SetAt(max, axis, fmax(vec3At(max, axis), vec3At(t.bbox.second, axis)));
        }
    }
    this->bbox = std::make_pair(min, max);
    return this->bbox;
}

glm::vec3 Mesh::Centre() {
    return (MeshBoundingBox().first + MeshBoundingBox().second) / 2.0f;
}

void Mesh::Translate(glm::vec3 translation) {
    for (glm::vec3& vertex : vertices) {
        vertex += translation;
    }
    ComputeBoundingBoxes();
}

void Mesh::Scale(glm::vec3 scale) {
    glm::vec3 centre = this->Centre();
    this->Translate(-centre);
    for (glm::vec3& vertex : vertices) {
        vertex = ScaleVec3(vertex, scale);
    }
    this->Translate(centre);
    ComputeBoundingBoxes();
}

void Mesh::Rotate(glm::vec3 rotation) {
    glm::mat3 rotMat = RotationMatrix(rotation.x, rotation.y, rotation.z);
    glm::vec3 centre = this->Centre();
    this->Translate(-centre);
    for (glm::vec3& vertex : vertices) {
        vertex = rotMat * vertex;
    }
    for (glm::vec3& vertexNormal : vertexNormals) {
        vertexNormal = rotMat * vertexNormal;
    }
    for (Triangle& t : triangles) {
        t.normal = rotMat * t.normal;
    }
    this->Translate(centre);
    ComputeBoundingBoxes();
}

void Mesh::ComputeBoundingBoxes() {
    for (Triangle& t : triangles) {
        t.BoundingBox(vertices);
    }
}