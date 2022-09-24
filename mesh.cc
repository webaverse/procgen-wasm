#include "mesh.h"
#include "util.h"
#include <iostream>

void TerrainGeometry::pushPointMetadata(const Heightfield &fieldValue) {
  // nothing
}
uint8_t *TerrainGeometry::getBuffer() const {
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
    // biomesUvs1
    sizeof(uint32_t) +
    biomesUvs1.size() * sizeof(biomesUvs1[0]) +
    // biomesUvs2
    sizeof(uint32_t) +
    biomesUvs2.size() * sizeof(biomesUvs2[0]) +
    sizeof(uint32_t) +
    materialsWeights.size() * sizeof(materialsWeights[0]) +
    // seeds
    // sizeof(uint32_t) +
    // seeds.size() * sizeof(seeds[0]) +
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

  // biomesUvs1
  *((uint32_t *)(buffer + index)) = biomesUvs1.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &biomesUvs1[0], biomesUvs1.size() * sizeof(biomesUvs1[0]));
  index += biomesUvs1.size() * sizeof(biomesUvs1[0]);

  // biomesUvs2
  *((uint32_t *)(buffer + index)) = biomesUvs2.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &biomesUvs2[0], biomesUvs2.size() * sizeof(biomesUvs2[0]));
  index += biomesUvs2.size() * sizeof(biomesUvs2[0]);

  // seeds
  // *((uint32_t *)(buffer + index)) = seeds.size();
  // index += sizeof(uint32_t);
  // std::memcpy(buffer + index, &seeds[0], seeds.size() * sizeof(seeds[0]));
  // index += seeds.size() * sizeof(seeds[0]);

  // materials weights
  *((uint32_t *)(buffer + index)) = materialsWeights.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &materialsWeights[0], materialsWeights.size() * sizeof(materialsWeights[0]));
  index += materialsWeights.size() * sizeof(materialsWeights[0]);

  // indices
  *((uint32_t *)(buffer + index)) = indices.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &indices[0], indices.size() * sizeof(indices[0]));
  index += indices.size() * sizeof(indices[0]);

  return buffer;
}

//

void WaterGeometry::pushPointMetadata(const Waterfield &fieldValue) {
  factors.push_back(fieldValue.waterFactor);
}
uint8_t *WaterGeometry::getBuffer() const {
  // calculate size
  size_t neededSize =
    // positions
    sizeof(uint32_t) +
    positions.size() * sizeof(positions[0]) +
    // normals
    sizeof(uint32_t) +
    normals.size() * sizeof(normals[0]) +
    // factors
    sizeof(uint32_t) +
    factors.size() * sizeof(factors[0]) +
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

  // factors
  *((uint32_t *)(buffer + index)) = factors.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &factors[0], factors.size() * sizeof(factors[0]));
  index += factors.size() * sizeof(factors[0]);

  // indices
  *((uint32_t *)(buffer + index)) = indices.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &indices[0], indices.size() * sizeof(indices[0]));
  index += indices.size() * sizeof(indices[0]);

  return buffer;
}

//

uint8_t *BarrierGeometry::getBuffer() const {
  // calculate size
  size_t neededSize =
    // positions
    sizeof(uint32_t) +
    positions.size() * sizeof(positions[0]) +
    // normals
    sizeof(uint32_t) +
    normals.size() * sizeof(normals[0]) +
    // uvs
    sizeof(uvs) +
    uvs.size() * sizeof(uvs[0]) +
    // positions2D
    sizeof(positions2D) +
    positions2D.size() * sizeof(positions2D[0]) +
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

  // uvs
  *((uint32_t *)(buffer + index)) = uvs.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &uvs[0], uvs.size() * sizeof(uvs[0]));
  index += uvs.size() * sizeof(uvs[0]);

  // positions2D
  *((uint32_t *)(buffer + index)) = positions2D.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &positions2D[0], positions2D.size() * sizeof(positions2D[0]));
  index += positions2D.size() * sizeof(positions2D[0]);

  // indices
  *((uint32_t *)(buffer + index)) = indices.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &indices[0], indices.size() * sizeof(indices[0]));
  index += indices.size() * sizeof(indices[0]);

  return buffer;
}