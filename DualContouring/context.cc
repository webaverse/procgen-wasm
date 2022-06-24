#include "context.h"
#include <iostream>

//

float TerrainDCContext::densityFn(const vm::vec3 &position, DCInstance *inst, Chunk3D &chunk) {
  return terrainDensityFn(position, inst, chunk);
}

//

float LiquidDCContext::densityFn(const vm::vec3 &position, DCInstance *inst, Chunk3D &chunk) {
  return liquidDensityFn(position, inst, chunk);
}