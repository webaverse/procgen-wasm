#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

constexpr int numThreads = NUM_THREADS;

constexpr int WORLD_BASE = 64;

constexpr int MIN_WORLD_HEIGHT = 0;
constexpr int MAX_WORLD_HEIGHT = 512;

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