#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <ctime>
#include <string.h>
#include <memory>

#include "libs/vector.h"
#include "generation/noises.h"
#include "generation/heightfield-generator.h"
#include "generation/instance-generator.h"
#include "polygonization/polygonizer.h"
#include "polygonization/mesh.h"
#include "task/task.h"
#include "task/promise.h"
#include "task/octree.h"

class ChunkResult
{
public:
    uint8_t *terrainMeshBuffer;
    uint8_t *waterMeshBuffer;
    uint8_t *treeInstancesBuffer;
    uint8_t *bushInstancesBuffer;
    uint8_t *rockInstancesBuffer;
    uint8_t *stoneInstancesBuffer;
    uint8_t *grassInstancesBuffer;
    uint8_t *poiInstancesBuffer;
    uint8_t *heightfieldsBuffer;

    void free();
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
    HeightfieldGenerator heightfieldGenerator;
    VegetationGenerator vegetationGenerator;
    PoiGenerator poiGenerator;
    Polygonizer polygonizer;
    RenderingInfo renderingInfo;

    PGInstance(int seed, int chunkSize);
    ~PGInstance();

    void setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix);

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
