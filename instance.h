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

//

class Tracker;
class OctreeNodeSpec;
class OctreeNode;

//

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

    NoiseField getNoise(int bx, int by);
    uint8_t getBiome(int bx, int bz);
    
    void getHeightField(int bx, int bz, int lod, Heightfield *heightfield);
    void getHeightFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, Heightfield *heightfieldSeams);
    Heightfield getHeightField(int bx, int bz);
    
    void getWaterField(int bx, int bz, int lod, float *waterField);
    float getWaterField(int bx, int bz);

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

    uint8_t *createGrassSplat(const vm::ivec2 &worldPositionXZ, const int lod);
    uint8_t *createVegetationSplat(const vm::ivec2 &worldPositionXZ, const int lod);
    uint8_t *createMobSplat(const vm::ivec2 &worldPositionXZ, const int lod);
    
    //
    
    uint8_t *createTerrainChunkMesh(const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray);
    uint8_t *createLiquidChunkMesh(const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray);

    //

    void setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix);
    void setClipRange(const vm::vec2 &min, const vm::vec2 &max);

    //

    float getTemperature(const vm::vec2 &worldPosition, const int &lod);
    float getHumidity(const vm::vec2 &worldPosition, const int &lod);

    //
    
    template<typename PositionType>
    bool tryLock(const PositionType &chunkPosition, int lod, GenerateFlags flags) {
        Mutex *chunkLock = getChunkLock(chunkPosition, lod, flags);
        return chunkLock->try_lock();
    }
    template<typename PositionType>
    void unlock(const PositionType &chunkPosition, int lod) {
        Mutex *chunkLock = getChunkLock(chunkPosition, lod);
        chunkLock->unlock();
    }

    //

    void createTerrainChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray);
    // void createLiquidChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray);

    //

    // void getHeightfieldRangeAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const vm::ivec2 &sizeXZ, int lod, float *heights, int priority);
    // void getLightRangeAsync(uint32_t id, const vm::ivec3 &worldPosition, const vm::ivec3 &size, int lod, uint8_t *skylights, uint8_t *aos, int priority);

    //

    void createGrassSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority);
    void createVegetationSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority);
    void createMobSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority);

    //

    // float signedDistanceToBox(float sx, float sy, float sz, float px, float py, float pz);
    // float signedDistanceToSphere(float cx, float cy, float cz, float r, float px, float py, float pz);
    
    //
    
    void trackerUpdateAsync(uint32_t id, Tracker *tracker, const vm::vec3 &position, int priority);
};

#endif // _INSTANCE_H_
