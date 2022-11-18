#include "biomes.h"

bool isLiquidBiome(unsigned char b) {
  return b == (unsigned char)BIOME::biOcean ||
    b == (unsigned char)BIOME::biRiver ||
    b == (unsigned char)BIOME::biFlowingRiver ||
    b == (unsigned char)BIOME::biSwampland ||
    b == (unsigned char)BIOME::biFrozenRiver ||
    b == (unsigned char)BIOME::biFrozenOcean || 
    b == (unsigned char)BIOME::biLava;
}