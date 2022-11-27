#include "instance.h"
#include "procgen.h"
#include "biomes.h"
#include "task/tracker.h"
#include "libs/vector.h"
#include "utils/util.h"
#include "libs/MurmurHash3.h"
#include <emscripten.h>

constexpr int CHUNK_RANGE = 1;

void ChunkResult::free(PGInstance *inst)
{
    std::free(terrainMeshBuffer);
    std::free(waterMeshBuffer);
    std::free(treeInstancesBuffer);
    std::free(bushInstancesBuffer);
    std::free(rockInstancesBuffer);
    std::free(stoneInstancesBuffer);
    std::free(grassInstancesBuffer);
    std::free(poiInstancesBuffer);
    std::free(heightfieldsBuffer);
    std::free(this);
}

// constructor/destructor
PGInstance::PGInstance(int seed, int chunkSize) : heightfieldGenerator(seed, chunkSize)

{
    // std::cout << "new pg instance " << seed << " " << chunkSize << std::endl;
}
PGInstance::~PGInstance() {}

uint8_t *PGInstance::createMobSplat(const vm::ivec2 &worldPositionXZ, const int lod)
{
    std::vector<float> ps;
    std::vector<float> qs;
    std::vector<float> instances;
    unsigned int count = 0;

    const int chunkSize = heightfieldGenerator.getChunkSize();

    int minX = worldPositionXZ.x / chunkSize * chunkSize;
    int minZ = worldPositionXZ.y / chunkSize * chunkSize;

    float seed = heightfieldGenerator.noises.mobNoise.in2D(minX, minZ);
    unsigned int seedInt;
    memcpy(&seedInt, &seed, sizeof(unsigned int));
    std::mt19937 rng(seedInt);

    const int maxNumMobs = 2;
    const float mobRate = 0.4;
    for (int i = 0; i < maxNumMobs; i++)
    {
        float dx = (float)rng() / (float)0xFFFFFFFF * (float)chunkSize;
        float dz = (float)rng() / (float)0xFFFFFFFF * (float)chunkSize;

        float ax = (float)minX + dx;
        float az = (float)minZ + dz;

        float noiseValue = heightfieldGenerator.noises.mobNoise.in2D(ax, az);

        if (noiseValue < mobRate)
        {
            const float height = heightfieldGenerator.getHeight(ax, az);

            ps.push_back(ax);
            ps.push_back(height);
            ps.push_back(az);

            Quat q = Quat().setFromAxisAngle(Vec{0, 1, 0}, rng() * 2.0f * M_PI);
            qs.push_back(q.x);
            qs.push_back(q.y);
            qs.push_back(q.z);
            qs.push_back(q.w);

            instances.push_back((float)rng() / (float)0xFFFFFFFF);

            count++;
        }
    }
    // serialize
    {
        const size_t size = sizeof(uint32_t) +
                            sizeof(float) * ps.size() +
                            sizeof(uint32_t) +
                            sizeof(float) * qs.size() +
                            sizeof(uint32_t) +
                            sizeof(float) * instances.size();

        uint8_t *buffer = (uint8_t *)malloc(size);
        int index = 0;

        ((uint32_t *)(buffer + index))[0] = ps.size();
        index += sizeof(uint32_t);
        memcpy(buffer + index, ps.data(), sizeof(float) * ps.size());
        index += sizeof(float) * ps.size();

        ((uint32_t *)(buffer + index))[0] = qs.size();
        index += sizeof(uint32_t);
        memcpy(buffer + index, qs.data(), sizeof(float) * qs.size());
        index += sizeof(float) * qs.size();

        ((uint32_t *)(buffer + index))[0] = instances.size();
        index += sizeof(uint32_t);
        memcpy(buffer + index, instances.data(), sizeof(float) * instances.size());
        index += sizeof(float) * instances.size();

        return buffer;
    }
}
//

inline int getLodInt(int lod)
{
    return 1 << (lod - 1);
}
inline int getLodRange(int lod, int chunkSize)
{
    return getLodInt(lod) * chunkSize;
}
OctreeContext PGInstance::getChunkSeedOctree(
    const vm::ivec2 &worldPosition,
    // int lod,
    int minLod,
    int maxLod, // we will sample a 3x3 of this lod
    int chunkSize)
{
    const int maxLodInt = getLodInt(maxLod);
    const int maxLodRange = getLodRange(maxLod, chunkSize);

    constexpr int maxNumSplits = 3;

    std::vector<vm::ivec2> maxLodChunkPositions;
    for (int dz = -1; dz <= 1; dz++)
    {
        for (int dx = -1; dx <= 1; dx++)
        {
            vm::ivec2 baseNode{
                (int)std::floor(
                    (float)(((float)worldPosition.x) / (float)maxLodRange) + (float)dx) *
                    maxLodInt,
                (int)std::floor(
                    (float)(((float)worldPosition.y) / (float)maxLodRange) + (float)dz) *
                    maxLodInt};

            // insert the node if it does not exist
            auto iter = std::find(
                maxLodChunkPositions.begin(),
                maxLodChunkPositions.end(),
                baseNode);
            if (iter == maxLodChunkPositions.end())
            {
                maxLodChunkPositions.push_back(baseNode);
            }
            else
            {
                std::cerr << "ERROR: duplicate node found: " << baseNode.x << ", " << baseNode.y << " : " << dx << " " << dz << std::endl;
                abort();
            }
        }
    }

    // compute splits
    std::vector<std::pair<vm::ivec2, int>> lodSplits;
    for (size_t i = 0; i < maxLodChunkPositions.size(); i++)
    {
        const vm::ivec2 &baseNode = maxLodChunkPositions[i];

        float chunkSeed = heightfieldGenerator.noises.numSplitsNoise.in2D(baseNode.x, baseNode.y);
        unsigned int seedInt = *(unsigned int *)&chunkSeed;
        std::mt19937 rng(seedInt);
        std::uniform_real_distribution<float> dis(0.f, 1.f);

        uint32_t numSplits = (uint32_t)(dis(rng) * (float)maxNumSplits);
        for (uint32_t i = 0; i < numSplits; i++)
        {
            uint32_t splitLodDX = (uint32_t)(dis(rng) * (float)maxLodInt);
            uint32_t splitLodDZ = (uint32_t)(dis(rng) * (float)maxLodInt);
            uint32_t splitLod = (uint32_t)(dis(rng) * (float)(maxLod - minLod) + minLod);

            int splitLodDXInt = baseNode.x + (int)splitLodDX;
            int splitLodDZInt = baseNode.y + (int)splitLodDZ;
            int splitLodInt = getLodInt(splitLod);

            lodSplits.push_back(
                std::make_pair(
                    vm::ivec2{
                        splitLodDXInt,
                        splitLodDZInt},
                    splitLodInt));
        }
    }

    OctreeContext octreeContext;
    constructSeedTree(
        octreeContext,
        maxLodChunkPositions,
        maxLodInt,
        lodSplits);
    return octreeContext;
}

void generateGridHeightfield(
    const std::vector<Heightfield> &heightfields,
    HeightfieldGeometry &heightfieldGeometry,
    int chunkSize)
{
    heightfieldGeometry.heightfieldImage.resize(chunkSize * chunkSize);

    int dstIndex = 0;
    const int chunkSizeP2 = chunkSize + 2;
    for (int y = 0; y < chunkSize; y++)
    {
        for (int x = 0; x < chunkSize; x++)
        {
            int dx = x + 1;
            int dy = y + 1;
            const int srcIndex = dx + dy * chunkSizeP2;
            const Heightfield &heightfield = heightfields[srcIndex];

            heightfieldGeometry.heightfieldImage[dstIndex++] = vm::vec4{
                heightfield.height,
                heightfield.liquidHeight,
                0,
                0};
        }
    }
}

void PGInstance::createChunkMesh(
    ChunkResult *result,
    const vm::ivec2 &worldPosition,
    int lod,
    const std::array<int, 2> &lodArray,
    int generateFlags,
    int numVegetationInstances,
    int numRockInstances,
    int numGrassInstances,
    int numPoiInstances)
{
    const int chunkSize = heightfieldGenerator.getChunkSize();

    // heightfield
    std::vector<Heightfield> heightfields;

    if (
        (generateFlags & GF_TERRAIN) |
        (generateFlags & GF_WATER) |
        (generateFlags & GF_VEGETATION) |
        (generateFlags & GF_ROCK) |
        (generateFlags & GF_GRASS) |
        (generateFlags & GF_POI) |
        (generateFlags & GF_HEIGHTFIELD))
    {
        heightfields = heightfieldGenerator.getHeightfields(worldPosition.x, worldPosition.y, lod, lodArray);
        polygonizer.calculateSurfaceNormals(heightfields, lod, lodArray, chunkSize);
        heightfieldGenerator.applyMaterials(worldPosition.x, worldPosition.y, lod, lodArray, heightfields);
    }

    // terrain
    if (generateFlags & GF_TERRAIN)
    {
        TerrainGeometry terrainGeometry;

        polygonizer.generateTerrainGeometry(
            worldPosition,
            lod,
            lodArray,
            chunkSize,
            heightfields,
            terrainGeometry);

        result->terrainMeshBuffer = terrainGeometry.getBuffer();
    }
    else
    {
        result->terrainMeshBuffer = nullptr;
    }

    // water
    if (generateFlags & GF_WATER)
    {
        WaterGeometry waterGeometry;

        const std::vector<Waterfield> &waterfields = *((std::vector<Waterfield> *)&heightfields);

        polygonizer.generateWaterGeometry(
            worldPosition,
            lod,
            lodArray,
            chunkSize,
            waterfields,
            waterGeometry);

        result->waterMeshBuffer = waterGeometry.getBuffer();
    }
    else
    {
        result->waterMeshBuffer = nullptr;
    }

    // vegetation
    if (generateFlags & GF_VEGETATION)
    {
        VegetationGeometry treeGeometry;
        VegetationGeometry bushGeometry;

        const uint8_t TREE_INSTANCE = (uint8_t)VEGETATION::TREE;

        const PushInstancesFunction pushInstancesFunction = [&](const float &ax, const float &az, const float &rot, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield, std::mt19937 &rng, std::uniform_real_distribution<float> &dis)
        {
            instanceGenerator.pushSplatInstances(
                ax,
                az,
                rot,
                treeGeometry,
                instanceId,
                heightfieldSampler);

            instanceGenerator.pushSubSplatInstances(
                ax,
                az,
                rot,
                bushGeometry,
                instanceId,
                heightfieldSampler,
                NUM_BUSHES_AROUND_TREE,
                rng,
                dis);
        };

        instanceGenerator.generateInstances<TREE_INSTANCE, VegetationGeometry>(
            worldPosition,
            lod,
            chunkSize,
            numVegetationInstances,
            MAX_NUM_VEGGIES_PER_CHUNK,
            heightfields,
            heightfieldGenerator.noises,
            pushInstancesFunction);

        result->treeInstancesBuffer = treeGeometry.getBuffer();
        result->bushInstancesBuffer = bushGeometry.getBuffer();
    }
    else
    {
        result->treeInstancesBuffer = nullptr;
        result->bushInstancesBuffer = nullptr;
    }

    // // rocks
    // if (generateFlags & GF_ROCK)
    // {
    //     SplatInstanceGeometry rockGeometry;
    //     SplatInstanceGeometry stoneGeometry;

    //     const uint8_t ROCK_INSTANCE = (uint8_t)LAYER::MINERALS;

    //     const PushInstancesFunction pushInstancesFunction = [&](const float &ax, const float &az, const float &rot, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield, std::mt19937 &rng, std::uniform_real_distribution<float> &dis)
    //     {
    //         instanceGenerator.pushSplatInstances(
    //             ax,
    //             az,
    //             rot,
    //             rockGeometry,
    //             instanceId,
    //             heightfieldSampler);

    //         instanceGenerator.pushSubSplatInstances(
    //             ax,
    //             az,
    //             rot,
    //             stoneGeometry,
    //             instanceId,
    //             heightfieldSampler,
    //             NUM_STONES_AROUND_ROCK,
    //             rng,
    //             dis);
    //     };

    //     instanceGenerator.generateInstances<ROCK_INSTANCE, SplatInstanceGeometry>(
    //         worldPosition,
    //         lod,
    //         chunkSize,
    //         numRockInstances,
    //         MAX_NUM_ROCKS_PER_CHUNK,
    //         heightfields,
    //         heightfieldGenerator.noises,
    //         pushInstancesFunction);

    //     result->rockInstancesBuffer = rockGeometry.getBuffer();
    //     result->stoneInstancesBuffer = stoneGeometry.getBuffer();
    // }
    // else
    // {
        result->rockInstancesBuffer = nullptr;
        result->stoneInstancesBuffer = nullptr;
    // }

    // // grass
    // if (generateFlags & GF_GRASS)
    // {
    //     GrassGeometry grassGeometry;
    //     const uint8_t GRASS_INSTANCE = (uint8_t)VEGETATION::GRASS;

    //     const PushInstancesFunction pushInstancesFunction = [&](const float &ax, const float &az, const float &rot, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield, std::mt19937 &rng, std::uniform_real_distribution<float> &dis)
    //     {
    //         instanceGenerator.pushGrassInstances(
    //             ax,
    //             az,
    //             rot,
    //             grassGeometry,
    //             instanceId,
    //             heightfieldSampler,
    //             heightfield,
    //             heightfieldGenerator.noises);
    //     };

    //     // TODO: get the cutoff value from JS
    //     const int GRASS_LOD_CUTOFF = 2;
    //     if (lod <= 2)
    //     {
    //         instanceGenerator.generateInstances<GRASS_INSTANCE, GrassGeometry>(
    //             worldPosition,
    //             lod,
    //             chunkSize,
    //             numGrassInstances,
    //             MAX_NUM_GRASSES_PER_CHUNK,
    //             heightfields,
    //             heightfieldGenerator.noises,
    //             pushInstancesFunction);
    //     }

    //     result->grassInstancesBuffer = grassGeometry.getBuffer();
    // }
    // else
    // {
        result->grassInstancesBuffer = nullptr;
    // }

    // poi
    if (generateFlags & GF_POI)
    {
        PoiGeometry poiGeometry;

        instanceGenerator.generatePoiInstances(
            worldPosition,
            lod,
            chunkSize,
            numPoiInstances,
            heightfields,
            heightfieldGenerator.noises,
            poiGeometry);

        result->poiInstancesBuffer = poiGeometry.getBuffer();
    }
    else
    {
        result->poiInstancesBuffer = nullptr;
    }

    // poi
    if (generateFlags & GF_HEIGHTFIELD)
    {
        HeightfieldGeometry heightfieldGeometry;

        generateGridHeightfield(
            heightfields,
            heightfieldGeometry,
            chunkSize);
        result->heightfieldsBuffer = heightfieldGeometry.getBuffer();
    }
    else
    {
        result->heightfieldsBuffer = nullptr;
    }
}
void PGInstance::createMobSplatAsync(
    uint32_t id,
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int priority)
{
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPositionXZ.y};
    Task *mobSplatTask = new Task(id, worldPositionF, lod, priority, [this, promise, worldPositionXZ, lod]() -> void
                                  {
        void *result = createMobSplat(worldPositionXZ, lod);
        promise->resolve(result); });
    ProcGen::taskQueue.pushTask(mobSplatTask);
}

//

uint8_t *PGInstance::createBarrierMesh(
    const vm::ivec2 &worldPosition,
    int minLod,
    int maxLod)
{
    const int chunkSize = heightfieldGenerator.getChunkSize();

    OctreeContext octreeContext = getChunkSeedOctree(
        worldPosition,
        minLod,
        maxLod,
        chunkSize);

    BarrierGeometry barrierGeometry;

    polygonizer.generateBarrierGeometry(
        worldPosition,
        chunkSize,
        octreeContext,
        barrierGeometry);

    uint8_t *result = barrierGeometry.getBuffer();

    return result;
}

void PGInstance::createBarrierMeshAsync(
    uint32_t id,
    const vm::ivec2 &worldPosition,
    int minLod,
    int maxLod)
{
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // this is so the chunk center is roughly the center of the one used by createBarrierMesh
    const int chunkSize = heightfieldGenerator.getChunkSize();

    const int maxLodRange = getLodRange(maxLod, chunkSize);
    vm::vec3 basePositionF{
        // center of the middle 3x3 chunks
        std::floor(
            (float)(((float)worldPosition.x) / (float)maxLodRange)) *
                maxLodRange +
            (maxLodRange / 2.f),
        (float)(-WORLD_BASE_HEIGHT) + ((float)MIN_WORLD_HEIGHT + (float)MAX_WORLD_HEIGHT) / 2.f,
        std::floor(
            (float)(((float)worldPosition.y + (float)maxLodRange / 2.f) / (float)maxLodRange)) *
                maxLodRange +
            (maxLodRange / 2.f),
    };
    const int maxLodP1 = maxLod + 1;
    Task *terrainTask = new Task(id, basePositionF, maxLodP1, [this, promise, worldPosition, minLod, maxLod]() -> void
                                 {
        uint8_t *result = createBarrierMesh(
            worldPosition,
            minLod,
            maxLod
        );
        if (!promise->resolve(result)) {
            free(result);
        } });
    ProcGen::taskQueue.pushTask(terrainTask);
}

void PGInstance::setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix)
{
    renderingInfo.worldPosition = worldPosition;
    renderingInfo.cameraPosition = cameraPosition;
    renderingInfo.cameraQuaternion = cameraQuaternion;
    renderingInfo.projectionMatrix = projectionMatrix;

    ProcGen::taskQueue.setCamera(worldPosition, cameraPosition, cameraQuaternion, projectionMatrix);
}

//

void PGInstance::createChunkMeshAsync(
    uint32_t id,
    const vm::ivec2 &worldPosition,
    int lod,
    const std::array<int, 2> &lodArray,
    int generateFlags,
    int numVegetationInstances,
    int numRockInstances,
    int numGrassInstances,
    int numPoiInstances)
{
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x,
        (float)(-WORLD_BASE_HEIGHT) + ((float)MIN_WORLD_HEIGHT + (float)MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y};
    std::array<int, 2> lodArray2{
        lodArray[0],
        lodArray[1]};
    Task *terrainTask = new Task(id, worldPositionF, lod, [this,
                                                           // result,
                                                           promise, worldPosition, lod, lodArray2, generateFlags, numVegetationInstances, numRockInstances, numGrassInstances, numPoiInstances]()
                                 {
        ChunkResult *result = (ChunkResult *)malloc(sizeof(ChunkResult));

        createChunkMesh(
            result,
            worldPosition,
            lod,
            lodArray2,
            generateFlags,
            numVegetationInstances,
            numRockInstances,
            numGrassInstances,
            numPoiInstances
        );
        if (!promise->resolve(result)) {
            result->free(this);
        } });
    ProcGen::taskQueue.pushTask(terrainTask);
}
// 2d caches

void PGInstance::trackerUpdateAsync(
    uint32_t id,
    Tracker *tracker,
    const vm::vec3 &position,
    int minLod,
    int maxLod,
    int lod1Range,
    int priority)
{
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // std::cout << "tracker update async " << priority << std::endl;

    Task *trackerUpdateTask = new Task(id, priority, [this, promise, tracker, position, minLod, maxLod, lod1Range]() -> void
                                       {
        const TrackerUpdate &trackerUpdate = tracker->update(this, position, minLod, maxLod, lod1Range);
        uint8_t *buffer = trackerUpdate.getBuffer();
        // std::cout << "trakcer update buffer address" << (void *)buffer << std::endl;
        if (!promise->resolve(buffer)) {
            free(buffer);
        } });
    ProcGen::taskQueue.pushTask(trackerUpdateTask);
}