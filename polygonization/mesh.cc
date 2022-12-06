#include "mesh.h"
#include "../task/octree.h"
#include "../utils/util.h"
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
    // materials
    sizeof(uint32_t) +
    materials.size() * sizeof(materials[0]) +
    // materials weights
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

  // materials
  *((uint32_t *)(buffer + index)) = materials.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &materials[0], materials.size() * sizeof(materials[0]));
  index += materials.size() * sizeof(materials[0]);

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
  factors.push_back(fieldValue.field.liquidHeight);
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
    // liquids
    sizeof(uint32_t) +
    liquids.size() * sizeof(liquids[0]) +
    // liquids weights
    sizeof(uint32_t) +
    liquidsWeights.size() * sizeof(liquidsWeights[0]) +
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

  // liquids
  *((uint32_t *)(buffer + index)) = liquids.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &liquids[0], liquids.size() * sizeof(liquids[0]));
  index += liquids.size() * sizeof(liquids[0]);

  // liquids weights
  *((uint32_t *)(buffer + index)) = liquidsWeights.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &liquidsWeights[0], liquidsWeights.size() * sizeof(liquidsWeights[0]));
  index += liquidsWeights.size() * sizeof(liquidsWeights[0]);

  // indices
  *((uint32_t *)(buffer + index)) = indices.size();
  index += sizeof(uint32_t);
  std::memcpy(buffer + index, &indices[0], indices.size() * sizeof(indices[0]));
  index += indices.size() * sizeof(indices[0]);

  return buffer;
}

//

uint8_t *SplatInstanceGeometry::getBuffer() const {
  // serialize
  size_t size = sizeof(uint32_t); // numInstances
  for (auto &iter : instances) {
      const SplatInstance &instance = iter.second;

      size += sizeof(int); // instanceId
      
      size += sizeof(uint32_t); // numPs
      size += sizeof(float) * instance.ps.size(); // ps
      
      size += sizeof(uint32_t); // numQs
      size += sizeof(float) * instance.qs.size(); // qs
  }

  uint8_t *buffer = (uint8_t *)malloc(size);
  int index = 0;

  *((uint32_t *)(buffer + index)) = instances.size();
  index += sizeof(uint32_t);
  
  for (auto &iter : instances) {
      const SplatInstance &instance = iter.second;

      *((int *)(buffer + index)) = instance.instanceId;
      index += sizeof(int);

      *((uint32_t *)(buffer + index)) = instance.ps.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.ps.data(), sizeof(float) * instance.ps.size());
      index += sizeof(float) * instance.ps.size();

      *((uint32_t *)(buffer + index)) = instance.qs.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.qs.data(), sizeof(float) * instance.qs.size());
      index += sizeof(float) * instance.qs.size();
  }

  return buffer;
}

uint8_t *MaterialAwareSplatInstanceGeometry::getBuffer() const {
  // serialize
  size_t size = sizeof(uint32_t); // numInstances
  for (auto &iter : instances) {
      const MaterialAwareSplatInstance &instance = iter.second;

      size += sizeof(int); // instanceId
      
      size += sizeof(uint32_t); // numPs
      size += sizeof(float) * instance.ps.size(); // ps
      
      size += sizeof(uint32_t); // numQs
      size += sizeof(float) * instance.qs.size(); // qs

      size += sizeof(uint32_t); // num materials
      size += sizeof(instance.materials[0]) * instance.materials.size(); // materials

      size += sizeof(uint32_t); // num materials weights
      size += sizeof(instance.materialsWeights[0]) * instance.materialsWeights.size(); // materials weights
  }

  uint8_t *buffer = (uint8_t *)malloc(size);
  int index = 0;

  *((uint32_t *)(buffer + index)) = instances.size();
  index += sizeof(uint32_t);
  
  for (auto &iter : instances) {
      const MaterialAwareSplatInstance &instance = iter.second;

      *((int *)(buffer + index)) = instance.instanceId;
      index += sizeof(int);

      *((uint32_t *)(buffer + index)) = instance.ps.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.ps.data(), sizeof(float) * instance.ps.size());
      index += sizeof(float) * instance.ps.size();

      *((uint32_t *)(buffer + index)) = instance.qs.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.qs.data(), sizeof(float) * instance.qs.size());
      index += sizeof(float) * instance.qs.size();

      *((uint32_t *)(buffer + index)) = instance.materials.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.materials.data(), sizeof(instance.materials[0]) * instance.materials.size());
      index += sizeof(instance.materials[0]) * instance.materials.size();
      
      *((uint32_t *)(buffer + index)) = instance.materialsWeights.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.materialsWeights.data(), sizeof(instance.materialsWeights[0]) * instance.materialsWeights.size());
      index += sizeof(instance.materialsWeights[0]) * instance.materialsWeights.size();
  }

  return buffer;
}

uint8_t *GrassGeometry::getBuffer() const {
  // serialize
  size_t size = sizeof(uint32_t); // numInstances
  for (auto &iter : instances) {
      const GrassSplatInstance &instance = iter.second;

      size += sizeof(int); // instanceId
      
      size += sizeof(uint32_t); // numPs
      size += sizeof(float) * instance.ps.size(); // ps
      
      size += sizeof(uint32_t); // numQs
      size += sizeof(float) * instance.qs.size(); // qs

      size += sizeof(uint32_t); // num materials
      size += sizeof(instance.materials[0]) * instance.materials.size(); // materials

      size += sizeof(uint32_t); // num materials weights
      size += sizeof(instance.materialsWeights[0]) * instance.materialsWeights.size(); // materials weights

      size += sizeof(uint32_t); // num grass props
      size += sizeof(instance.grassProps[0]) * instance.grassProps.size(); // grass props
  }

  uint8_t *buffer = (uint8_t *)malloc(size);
  int index = 0;

  *((uint32_t *)(buffer + index)) = instances.size();
  index += sizeof(uint32_t);
  
  for (auto &iter : instances) {
      const GrassSplatInstance &instance = iter.second;

      *((int *)(buffer + index)) = instance.instanceId;
      index += sizeof(int);

      *((uint32_t *)(buffer + index)) = instance.ps.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.ps.data(), sizeof(float) * instance.ps.size());
      index += sizeof(float) * instance.ps.size();

      *((uint32_t *)(buffer + index)) = instance.qs.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.qs.data(), sizeof(float) * instance.qs.size());
      index += sizeof(float) * instance.qs.size();

      *((uint32_t *)(buffer + index)) = instance.materials.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.materials.data(), sizeof(instance.materials[0]) * instance.materials.size());
      index += sizeof(instance.materials[0]) * instance.materials.size();
      
      *((uint32_t *)(buffer + index)) = instance.materialsWeights.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.materialsWeights.data(), sizeof(instance.materialsWeights[0]) * instance.materialsWeights.size());
      index += sizeof(instance.materialsWeights[0]) * instance.materialsWeights.size();

      *((uint32_t *)(buffer + index)) = instance.grassProps.size();
      index += sizeof(uint32_t);
      memcpy(buffer + index, instance.grassProps.data(), sizeof(instance.grassProps[0]) * instance.grassProps.size());
      index += sizeof(instance.grassProps[0]) * instance.grassProps.size();
  }

  return buffer;
}

//

uint8_t *PoiGeometry::getBuffer() const {
  // serialize
  size_t size = sizeof(uint32_t) + // numPs
    sizeof(float) * ps.size() + // ps
    sizeof(uint32_t) + // numInstances
    sizeof(int32_t) * instances.size();

  uint8_t *buffer = (uint8_t *)malloc(size);
  int index = 0;

  *((uint32_t *)(buffer + index)) = ps.size();
  index += sizeof(uint32_t);

  memcpy(buffer + index, ps.data(), sizeof(float) * ps.size());
  index += sizeof(float) * ps.size();

  *((uint32_t *)(buffer + index)) = instances.size();
  index += sizeof(uint32_t);

  memcpy(buffer + index, instances.data(), sizeof(int32_t) * instances.size());
  index += sizeof(int32_t) * instances.size();

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
    indices.size() * sizeof(indices[0]) + 
    // numLeafNodes
    sizeof(uint32_t) +
    leafNodes.size() * (sizeof(vm::ivec2) + sizeof(int)) +
    // leafNodesMin
    sizeof(leafNodesMin) +
    // leafNodesMax
    sizeof(leafNodesMax) +
    // leafNodesIndex
    leafNodesIndex.size() * sizeof(leafNodesIndex[0]);

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

  // leaf nodes
  *((uint32_t *)(buffer + index)) = leafNodes.size();
  index += sizeof(uint32_t);
  for (size_t i = 0; i < leafNodes.size(); i++) {
    OctreeNodePtr nodePtr = leafNodes[i];
    OctreeNode &node = *nodePtr;

    std::memcpy(buffer + index, &node.min, sizeof(vm::ivec2));
    index += sizeof(vm::ivec2);
    std::memcpy(buffer + index, &node.lod, sizeof(int));
    index += sizeof(int);
  }

  // leafNodesMin
  std::memcpy(buffer + index, &leafNodesMin, sizeof(leafNodesMin));
  index += sizeof(leafNodesMin);
  
  // leafNodesMax
  std::memcpy(buffer + index, &leafNodesMax, sizeof(leafNodesMax));
  index += sizeof(leafNodesMax);
  
  // leafNodesIndex
  memcpy(buffer + index, leafNodesIndex.data(), leafNodesIndex.size() * sizeof(leafNodesIndex[0]));
  index += leafNodesIndex.size() * sizeof(leafNodesIndex[0]);

  return buffer;
}

uint8_t *HeightfieldGeometry::getBuffer() const {
  // calculate size
  size_t neededSize = 
    sizeof(uint32_t) + // numPixels
    heightfieldImage.size() * sizeof(heightfieldImage[0]); // pixels

  // allocate buffer
  uint8_t *buffer = (uint8_t *)malloc(neededSize);
  int index = 0;

  // numPixels
  *((uint32_t *)(buffer + index)) = heightfieldImage.size();
  index += sizeof(uint32_t);

  // pixels
  std::memcpy(buffer + index, heightfieldImage.data(), heightfieldImage.size() * sizeof(heightfieldImage[0]));
  index += heightfieldImage.size() * sizeof(heightfieldImage[0]);

  return buffer;
}

uint8_t *SplatInstanceGeometryGroup::getBuffer() const {
  // serialize
  size_t size = sizeof(uint32_t); // num geometries

  for (auto &iter : geometries) {
      size += sizeof(uint32_t); // geometry buffer address
  }

  uint8_t *buffer = (uint8_t *)malloc(size);
  int index = 0;

  *((uint32_t *)(buffer + index)) = geometries.size(); // num geometries
  index += sizeof(uint32_t);
  
  for (auto &iter : geometries) {
      const SplatInstanceGeometry &geometry = iter;

      *((uint32_t *)(buffer + index)) = (uint32_t)geometry.getBuffer(); // geometry buffer address
      index += sizeof(uint32_t);
  }

  return buffer;
}
