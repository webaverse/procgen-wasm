#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

constexpr int numThreads = NUM_THREADS;

constexpr int WORLD_BASE_HEIGHT = 128;

constexpr int MIN_WORLD_HEIGHT = 0;
constexpr int MAX_WORLD_HEIGHT = 1024;

// * Terrain parameters
constexpr float COLD_WARM_BORDER = 0.5f;
constexpr float WARM_HOT_BORDER = 0.1f;
constexpr int MAX_NUM_GRASSES_PER_CHUNK = 2048;
constexpr float GRASS_THRESHOLD = 0.3f;
constexpr int MAX_NUM_VEGGIES_PER_CHUNK = 8;
constexpr float VEGGIE_THRESHOLD = 0.65;
constexpr float OCEAN_THRESHOLD = 0.85f;
constexpr float RIVER_THRESHOLD = 0.85f;
constexpr float RIVER_BASE = 0.3f;
constexpr float WATER_THRESHOLD = OCEAN_THRESHOLD;
constexpr float ROCK_THRESHOLD = 0.92f;
constexpr int NUM_STONES_AROUND_ROCK = 2;
constexpr int WATER_BASE_HEIGHT = WORLD_BASE_HEIGHT;
constexpr int OCEAN_DEPTH = WATER_BASE_HEIGHT * 8;
constexpr int RIVER_DEPTH = WATER_BASE_HEIGHT / 6;
constexpr int WATER_BASE_DEPTH = 8;
constexpr float WATER_OFFSET = 1.5f;
constexpr float WATER_HEIGHT_DIFFERENCE = 3.f;
constexpr int NUM_BUSHES_AROUND_TREE = 2;
constexpr float BUSH_AROUND_TREE_BASE_OFFSET = 2.f;
constexpr float BUSH_AROUND_TREE_OFFSET_RANGE = 3.f;
constexpr float GRASS_HEIGHT_VARIATION_RANGE = 0.3f;
constexpr float GRASS_COLOR_VARIATION_BASE = 1.f;
constexpr float GRASS_COLOR_VARIATION_RANGE = 0.2f;
constexpr float CRUSHED_GRASS_THRESHOLD = 1.f - ROCK_THRESHOLD;
constexpr float GRASS_MODEL_BASE_HEIGHT = 1.2f;

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