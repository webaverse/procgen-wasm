#include "mesh.h"
#include <iostream>

uint8_t *TerrainVertexBuffer::getBuffer() const {
  // calculate size
  size_t neededSize = 
    // positions
    sizeof(uint32_t) +
    positions.size() * sizeof(positions[0]) +
    // normals
    sizeof(uint32_t) +
    normals.size() * sizeof(normals[0]) +
    // biomes
    sizeof(uint32_t) +
    biomes.size() * sizeof(biomes[0]) +
    // biomesWeights
    sizeof(uint32_t) +
    biomesWeights.size() * sizeof(biomesWeights[0]) +
    // biomesUvs
    sizeof(uint32_t) +
    biomesUvs.size() * sizeof(biomesUvs[0]) +
    // indices
    sizeof(uint32_t) +
    indices.size() * sizeof(indices[0]);

  // allocate buffer
  uint8_t *buffer = (uint8_t *)malloc(neededSize);
  int index = 0;

  // positions
  *((uint32_t *)(buffer + index)) = positions.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &positions[0], positions.size() * sizeof(positions[0]));
  index += positions.size() * sizeof(positions[0]);

  // normals
  *((uint32_t *)(buffer + index)) = normals.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &normals[0], normals.size() * sizeof(normals[0]));
  index += normals.size() * sizeof(normals[0]);

  // biomes
  *((uint32_t *)(buffer + index)) = biomes.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &biomes[0], biomes.size() * sizeof(biomes[0]));
  index += biomes.size() * sizeof(biomes[0]);

  // biomesWeights
  *((uint32_t *)(buffer + index)) = biomesWeights.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &biomesWeights[0], biomesWeights.size() * sizeof(biomesWeights[0]));
  index += biomesWeights.size() * sizeof(biomesWeights[0]);

  // biomesUvs
  *((uint32_t *)(buffer + index)) = biomesUvs.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &biomesUvs[0], biomesUvs.size() * sizeof(biomesUvs[0]));
  index += biomesUvs.size() * sizeof(biomesUvs[0]);

  // indices
  *((uint32_t *)(buffer + index)) = indices.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &indices[0], indices.size() * sizeof(indices[0]));
  index += indices.size() * sizeof(indices[0]);

  return buffer;
}
void TerrainVertexBuffer::pushVertexData(const VertexData &vertexData) {
  positions.push_back(vertexData.position);
  normals.push_back(vertexData.normal);
  biomes.push_back(vertexData.biomes);
  biomesWeights.push_back(vertexData.biomesWeights);
  biomesUvs.push_back(vertexData.biomeUvs);
  skylights.push_back(vertexData.skylight);
  aos.push_back(vertexData.ao);
}

//

uint8_t *LiquidVertexBuffer::getBuffer() const {
  // calculate size
  size_t neededSize =
    // positions
    sizeof(uint32_t) +
    positions.size() * sizeof(positions[0]) +
    // normals
    sizeof(uint32_t) +
    normals.size() * sizeof(normals[0]) +
    // biomes
    sizeof(uint32_t) +
    biomes.size() * sizeof(biomes[0]) +
    // indices
    sizeof(uint32_t) +
    indices.size() * sizeof(indices[0]);

  // allocate buffer
  uint8_t *buffer = (uint8_t *)malloc(neededSize);
  int index = 0;

  // positions
  *((uint32_t *)(buffer + index)) = positions.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &positions[0], positions.size() * sizeof(positions[0]));
  index += positions.size() * sizeof(positions[0]);

  // normals
  *((uint32_t *)(buffer + index)) = normals.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &normals[0], normals.size() * sizeof(normals[0]));
  index += normals.size() * sizeof(normals[0]);

  // biomes
  *((uint32_t *)(buffer + index)) = biomes.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &biomes[0], biomes.size() * sizeof(biomes[0]));
  index += biomes.size() * sizeof(biomes[0]);

  // indices
  *((uint32_t *)(buffer + index)) = indices.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &indices[0], indices.size() * sizeof(indices[0]));
  index += indices.size() * sizeof(indices[0]);

  return buffer;
}
void LiquidVertexBuffer::pushVertexData(const VertexData &vertexData) {
  positions.push_back(vertexData.position);
  normals.push_back(vertexData.normal);
  biomes.push_back(vertexData.biomes.x);
}