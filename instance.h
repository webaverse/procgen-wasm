#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <ctime>
#include <string.h>
#include <memory>
#include "context.h"
#include "mesh.h"
#include "task.h"
#include "vector.h"
#include "noises.h"
#include "tracker.h"
#include "promise.h"
#include "./generation/generator.h"

//

class Tracker;
class OctreeNodeSpec;
class OctreeNode;
class OctreeContext;

//

class ChunkResult
{
public:
    uint8_t *terrainMeshBuffer = nullptr;
    uint8_t *waterMeshBuffer = nullptr;
    uint8_t *treeInstancesBuffer = nullptr;
    uint8_t *bushInstancesBuffer = nullptr;
    uint8_t *rockInstancesBuffer = nullptr;
    uint8_t *stoneInstancesBuffer = nullptr;
    uint8_t *grassInstancesBuffer = nullptr;
    uint8_t *poiInstancesBuffer = nullptr;
    uint8_t *heightfieldsBuffer = nullptr;

    void free(PGInstance *inst);
};


class RenderingInfo
{
public:
    vm::vec3 worldPosition;
    vm::vec3 cameraPosition;
    Quat cameraQuaternion;
    std::array<float, 16> projectionMatrix;
};

class PGInstance
{
public:
    std::unique_ptr<vm::box3> clipRange;

    Generator generator;
    RenderingInfo renderingInfo;

    PGInstance(int seed, int chunkSize);
    ~PGInstance();

    void setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix);
    void setClipRange(const vm::vec2 &min, const vm::vec2 &max);

    OctreeContext getChunkSeedOctree(
        const vm::ivec2 &worldPosition,
        int minLod,
        int maxLod,
        int chunkSize);

    uint8_t *createMobSplat(const vm::ivec2 &worldPositionXZ, const int lod);

    void createChunkMesh(
        ChunkResult *result,
        const vm::ivec2 &worldPosition,
        int lod,
        const std::array<int, 2> &lodArray,
        int generateFlags,
        int numVegetationInstances,
        int numRockInstances,
        int numGrassInstances,
        int numPoiInstances);


    void createChunkMeshAsync(
        uint32_t id,
        const vm::ivec2 &worldPosition,
        int lod,
        const std::array<int, 2> &lodArray,
        int generateFlags,
        int numVegetationInstances,
        int numRockInstances,
        int numGrassInstances,
        int numPoiInstances);

    void createMobSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority);

    uint8_t *createBarrierMesh(
        const vm::ivec2 &worldPosition,
        int minLod,
        int maxLod);

    void createBarrierMeshAsync(
        uint32_t id,
        const vm::ivec2 &worldPosition,
        int minLod,
        int maxLod);

    void trackerUpdateAsync(
        uint32_t id,
        Tracker *tracker,
        const vm::vec3 &position,
        int minLod,
        int maxLod,
        int lod1Range,
        int priority);
};

#endif // _INSTANCE_H_
