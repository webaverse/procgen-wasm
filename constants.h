#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

constexpr int numThreads = NUM_THREADS;

constexpr int WORLD_BASE_HEIGHT = 128;

constexpr int MIN_WORLD_HEIGHT = 0;
constexpr int MAX_WORLD_HEIGHT = 2048;

// * Terrain parameters
constexpr float OCEAN_THRESHOLD = 0.8f;
constexpr int WATER_BASE_HEIGHT = WORLD_BASE_HEIGHT;

//

constexpr float BIOME_DEBUG_MESH_BASE_HEIGHT = 0.1;

constexpr int CAVE_BASE_HEIGHT = 75;

// constexpr int chunkSize = 16;
// constexpr int CHUNK_CACHE_RANGE = 8;
// constexpr int cacheWidth = chunkSize * CHUNK_CACHE_RANGE;

constexpr double frustumCullDistancePenalty = 100000; // 100km
constexpr double priorityDistancePenalty = 100000; // 100km

/* enum class PEEK_FACES : int {
  FRONT = 0,
  BACK,
  LEFT,
  RIGHT,
  TOP,
  BOTTOM,
  NONE
}; */

#endif // _CONSTANTS_H_