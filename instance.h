#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <ctime>
#include <string.h>
#include <memory>
#include "chunk.h"
#include "context.h"
#include "task.h"
#include "vector.h"
#include "noises.h"
#include "tracker.h"

//

class Tracker;
class OctreeNodeSpec;
class OctreeNode;

//

class ChunkResult {
public:
    uint8_t *terrainMeshBuffer;
    uint8_t *waterMeshBuffer;
    uint8_t *barrierMeshBuffer;

    void free() {
        std::free(terrainMeshBuffer);
        std::free(waterMeshBuffer);
        std::free(barrierMeshBuffer);
        std::free(this);
    }
};

class SeedNoise {
public:
    float seed;
    float seedRadius;
};

class Line2 {
public:
    vm::vec2 start;
    vm::vec2 end;

    float distanceToPoint(const vm::vec2 &p) const {
        vm::vec2 v = end - start;
        vm::vec2 w = p - start;

        float c1 = vm::dot(w, v);
        if (c1 <= 0)
            return vm::length(p - start);

        float c2 = vm::dot(v, v);
        if (c2 <= c1)
            return vm::length(p - end);

        float b = c1 / c2;
        vm::vec2 pb = start + v * b;
        return vm::length(p - pb);
    }
};

class PGInstance {
public:
    int seed;
    int chunkSize;
    Noises noises;
    std::unique_ptr<vm::box3> clipRange;

    vm::vec3 worldPosition;
    vm::vec3 cameraPosition;
    Quat cameraQuaternion;
    std::array<float, 16> projectionMatrix;

    // 2d caches

    NoiseField getNoise(float bx, float by);
    uint8_t getBiome(float bx, float bz);
    // SeedNoise getSeedNoise(int bx, int bz);

    void getHeightFieldCenter(int bx, int bz, int lod, std::vector<Heightfield> &heightfield);
    void getHeightFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, std::vector<Heightfield> &heightfieldSeams);
    Heightfield getHeightField(float bx, float bz);
    float getHeight(float bx, float bz);
    
    void getWaterFieldCenter(int bx, int bz, int lod, std::vector<Waterfield> &waterfield);
    void getWaterFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, std::vector<Waterfield> &waterfieldSeams);
    Waterfield getWaterField(int bx, int bz, int lod);

    // void getCaveFieldCenter(int bx, int bz, int lod, std::vector<Cavefield> &cavefields);
    // void getCaveFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, std::vector<Cavefield> &cavefields);
    // Cavefield getCavefield(int bx, int bz);
    // void getCaveFieldChunk(int bx, int bz, int lod, const std::array<int, 2> &lodArray, std::vector<Cavefield> &cavefields);
    // Cavefield getCavefield(int bx, int bz, const std::vector<Line2> &caveLines);

    // float getSeed(int bx, int bz);

    float getComputedBiomeHeight(unsigned char b, const vm::vec2 &worldPosition);

    //

    PGInstance(int seed, int chunkSize);
    ~PGInstance();

    //
    
    // peek buffer

    /* Chunk3D &getChunk(const vm::ivec3 &min, const int lod, GenerateFlags flags);
    Chunk3D &getChunkInternal(const vm::ivec3 &min, int lod);
    Chunk3D &getChunkAt(const float x, const float y, const float z, const int lod, GenerateFlags flags);

    Chunk2D &getChunk(const vm::ivec2 &min, const int lod, GenerateFlags flags);
    Chunk2D &getChunkInternal(const vm::ivec2 &min, int lod);
    Chunk2D &getChunkAt(const float x, const float z, const int lod, GenerateFlags flags);
    
    //

    Mutex *getChunkLock(const vm::ivec2 &worldPos, const int lod, const int flags);
    Mutex *getChunkLock(const vm::ivec3 &worldPos, const int lod); */

    //

    // void getHeightfieldRange(const vm::ivec2 &worldPositionXZ, const vm::ivec2 &size, int lod, float *heights);
    // void getLightRange(const vm::ivec3 &worldPosition, const vm::ivec3 &size, int lod, uint8_t *skyLights, uint8_t *aos);

    //

    /* float *getChunkHeightfield(const vm::ivec2 &worldPositionXZ, int lod);
    unsigned char *getChunkSkylight(const vm::ivec2 &worldPosition, int lod);
    unsigned char *getChunkAo(const vm::ivec2 &worldPosition, int lod); */

    //

    uint8_t *createChunkGrass(const vm::ivec2 &worldPositionXZ, const int lod, const int numGrassInstances);
    uint8_t *createChunkVegetation(const vm::ivec2 &worldPositionXZ, const int lod, const int numVegetationInstances);
    uint8_t *createMobSplat(const vm::ivec2 &worldPositionXZ, const int lod);
    
    //
    
    ChunkResult *createChunkMesh(const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray, int generateFlags);
    // uint8_t *createLiquidChunkMesh(const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray);
    OctreeContext getChunkSeedOctree(const vm::ivec2 &worldPosition, int lod, int chunkSize);

    //

    void setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix);
    void setClipRange(const vm::vec2 &min, const vm::vec2 &max);

    //

    float getTemperature(const vm::vec2 &worldPosition, const int &lod);
    float getHumidity(const vm::vec2 &worldPosition, const int &lod);

    //
    
    /* template<typename PositionType>
    bool tryLock(const PositionType &chunkPosition, int lod, GenerateFlags flags) {
        Mutex *chunkLock = getChunkLock(chunkPosition, lod, flags);
        return chunkLock->try_lock();
    }
    template<typename PositionType>
    void unlock(const PositionType &chunkPosition, int lod) {
        Mutex *chunkLock = getChunkLock(chunkPosition, lod);
        chunkLock->unlock();
    } */

    //

    void createChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray, int generateFlags);
    // void createLiquidChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray);
    void createChunkGrassAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int numGrassInstances);
    void createChunkVegetationAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int numVegetationInstances);

    //

    // void getHeightfieldRangeAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const vm::ivec2 &sizeXZ, int lod, float *heights, int priority);
    // void getLightRangeAsync(uint32_t id, const vm::ivec3 &worldPosition, const vm::ivec3 &size, int lod, uint8_t *skylights, uint8_t *aos, int priority);

    //

    void createMobSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority);

    //

    // float signedDistanceToBox(float sx, float sy, float sz, float px, float py, float pz);
    // float signedDistanceToSphere(float cx, float cy, float cz, float r, float px, float py, float pz);
    
    //
    
    void trackerUpdateAsync(uint32_t id, Tracker *tracker, const vm::vec3 &position, int priority);
};

#endif // _INSTANCE_H_
