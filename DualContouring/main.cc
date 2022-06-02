#include "main.h"
#include "octree.h"

// namespace ChunkMesh
// {
//     OctreeNode *chunkRoot;
//     OctreeNode *chunkWithLod;
//     std::vector<OctreeNode *> neighbouringChunks;
//     OctreeNode *seamRoot;
// };

struct AllocationsMetrics
{
    uint32_t TotalAllocated = 0;
    uint32_t TotalFreed = 0;
    uint32_t CurrentUsage() { return TotalAllocated - TotalFreed; }
};

static AllocationsMetrics sAllocationsMetrics;

void *operator new(size_t size)
{
    sAllocationsMetrics.TotalAllocated += size;
    return std::malloc(size);
}
void operator delete(void *memory, size_t size)
{
    sAllocationsMetrics.TotalFreed += size;
    free(memory);
}

static void PrintMemoryUsage(){
    std::cout << "Memory Usage : " << sAllocationsMetrics.CurrentUsage() << std::endl;
}

namespace DualContouring
{
    // chunk settings
    int chunkSize = 16;
    FastNoise *fastNoise = nullptr;

    // storing the octrees that we would delete after mesh construction
    std::vector<OctreeNode *> neighbourNodes;

    // storing the octree roots here for search
    std::unordered_map<uint64_t, OctreeNode *> chunksListHashMap;
    std::unordered_map<uint64_t, Chunk> chunksNoiseHashMap;

    void initialize(int newChunkSize, int seed)
    {
        chunkSize = newChunkSize;

        std::mt19937 rng(seed);
        fastNoise = new FastNoise(rng());
        fastNoise->SetFrequency(0.02);
        fastNoise->SetFractalOctaves(4);
    }
    uint8_t *constructOutputBuffer(VertexBuffer &vertexBuffer)
    {
        return vertexBuffer.getBuffer();
    }
    Chunk &getChunk(const vm::ivec3 &min)
    {
        uint64_t minHash = hashOctreeMin(min);

        const auto &iter = chunksNoiseHashMap.find(minHash);
        if (iter == chunksNoiseHashMap.end())
        {
            chunksNoiseHashMap.emplace(std::make_pair(minHash, Chunk(min)));
        }

        Chunk &chunkNoise = chunksNoiseHashMap.find(minHash)->second;
        return chunkNoise;
    }
    float *getChunkHeightField(float x, float y, float z)
    {
        const vm::ivec3 min = vm::ivec3(x, y, z);
        Chunk &chunkNoise = getChunk(min);
        std::vector<float> &heightfieldRaw = chunkNoise.cachedHeightField;
        float *heightfieldOut = (float *)malloc(chunkSize * chunkSize * chunkSize * sizeof(float));
        int gridSize = chunkSize + 3;
        for (int x = 0; x < chunkSize; x++)
        {
            for (int z = 0; z < chunkSize; z++)
            {
                for (int y = 0; y < chunkSize; y++)
                {
                    int lx = x + 1;
                    int ly = y + 1;
                    int lz = z + 1;

                    int outIndex = x + z * chunkSize + y * chunkSize * chunkSize;
                    int inIndex = lx + lz * gridSize + ly * gridSize * gridSize;
                    heightfieldOut[outIndex] = heightfieldRaw[inIndex];
                }
            }
        }
        return heightfieldOut;
    }

    void clearChunkRoot(float x, float y, float z)
    {
        // we destroy the chunk root separately because we might need it for LOD switch if it's already generated
        const vm::ivec3 octreeMin = vm::ivec3(x, y, z);
        OctreeNode *chunkRoot = getChunkRootFromHashMap(octreeMin, chunksListHashMap);
        if (!chunkRoot)
        {
            return;
        }
        // sort and remove duplicates
        std::sort(neighbourNodes.begin(), neighbourNodes.end());
        neighbourNodes.erase(std::unique(neighbourNodes.begin(), neighbourNodes.end()), neighbourNodes.end());
        for (int i = 0; i < neighbourNodes.size(); i++)
        {
            const OctreeNode *node = neighbourNodes[i];
            if (node->drawInfo)
            {
                delete node->drawInfo;
            }
            delete node;
        }
        neighbourNodes.clear();
        removeOctreeFromHashMap(octreeMin, chunksListHashMap);
        destroyOctree(chunkRoot);
    }

    uint8_t *createChunkMesh(float x, float y, float z, const int lod)
    {
        VertexBuffer vertexBuffer;

        const vm::ivec3 octreeMin = vm::ivec3(x, y, z);
        // OctreeNode *chunkRoot = getChunkRootFromHashMap(octreeMin, chunksListHashMap);

        Chunk &chunk = getChunk(octreeMin);
        ChunkOctree chunkOctree(chunk, octreeMin, chunkSize, lod);
        if (!chunkOctree.root)
        {
            // printf("Chunk Has No Data\n");
            return nullptr;
        }

        generateMeshFromOctree(chunkOctree.root, lod, false, vertexBuffer);

        std::vector<OctreeNode *> seamNodes = generateSeamNodes(chunk, chunkOctree, neighbourNodes);
        OctreeNode *seamRoot = constructOctreeUpwards(seamRoot, seamNodes, chunkOctree.root->min, chunkOctree.root->size * 2);
        generateMeshFromOctree(seamRoot, lod, true, vertexBuffer);

        // mesh is not valid
        if (vertexBuffer.indices.size() == 0)
        {
            // printf("Generated Mesh Is Not Valid\n");
            return nullptr;
        }

        addChunkRootToHashMap(chunkOctree.root, chunksListHashMap);

        return constructOutputBuffer(vertexBuffer);
    }

    bool drawSphereDamage(const float &x, const float &y, const float &z, const float radius, float *outPositions, unsigned int *outPositionsCount, float *outDamages)
    {
        unsigned int maxPositionsCount = *outPositionsCount;
        *outPositionsCount = 0;

        bool drew = false;
        std::set<uint64_t> seenHashes;
        for (float dx = -1; dx <= 1; dx += 2)
        {
            for (float dz = -1; dz <= 1; dz += 2)
            {
                for (float dy = -1; dy <= 1; dy += 2)
                {
                    float ax = x + dx * radius;
                    float ay = y + dy * radius;
                    float az = z + dz * radius;
                    vm::ivec3 min = vm::ivec3(std::floor(ax / (float)chunkSize), std::floor(ay / (float)chunkSize), std::floor(az / (float)chunkSize)) * chunkSize;
                    uint64_t minHash = hashOctreeMin(min);
                    if (seenHashes.find(minHash) == seenHashes.end())
                    {
                        seenHashes.insert(minHash);

                        Chunk &chunkNoise = getChunk(min);
                        if (chunkNoise.addSphereDamage(ax, ay, az, radius))
                        {
                            if (*outPositionsCount < maxPositionsCount)
                            {
                                int gridSize = chunkSize + 3;
                                int damageBufferSize = gridSize * gridSize * gridSize;
                                memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.data(), sizeof(float) * damageBufferSize);

                                outPositions[(*outPositionsCount) * 3] = min.x;
                                outPositions[(*outPositionsCount) * 3 + 1] = min.y;
                                outPositions[(*outPositionsCount) * 3 + 2] = min.z;

                                (*outPositionsCount)++;
                            }

                            drew = true;
                        }
                    }
                }
            }
        }
        return drew;
    }

    bool eraseSphereDamage(const float &x, const float &y, const float &z, const float radius, float *outPositions, unsigned int *outPositionsCount, float *outDamages)
    {
        unsigned int maxPositionsCount = *outPositionsCount;
        *outPositionsCount = 0;

        bool drew = false;
        std::set<uint64_t> seenHashes;
        for (float dx = -1; dx <= 1; dx += 2)
        {
            for (float dz = -1; dz <= 1; dz += 2)
            {
                for (float dy = -1; dy <= 1; dy += 2)
                {
                    float ax = x + dx * radius;
                    float ay = y + dy * radius;
                    float az = z + dz * radius;
                    vm::ivec3 min = vm::ivec3(std::floor(ax / (float)chunkSize), std::floor(ay / (float)chunkSize), std::floor(az / (float)chunkSize)) * chunkSize;
                    uint64_t minHash = hashOctreeMin(min);
                    if (seenHashes.find(minHash) == seenHashes.end())
                    {
                        seenHashes.insert(minHash);

                        Chunk &chunkNoise = getChunk(min);
                        if (chunkNoise.removeSphereDamage(ax, ay, az, radius))
                        {
                            if (*outPositionsCount < maxPositionsCount)
                            {
                                int gridSize = chunkSize + 3;
                                int damageBufferSize = gridSize * gridSize * gridSize;
                                memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.data(), sizeof(float) * damageBufferSize);

                                outPositions[(*outPositionsCount) * 3] = min.x;
                                outPositions[(*outPositionsCount) * 3 + 1] = min.y;
                                outPositions[(*outPositionsCount) * 3 + 2] = min.z;

                                (*outPositionsCount)++;
                            }

                            drew = true;
                        }
                    }
                }
            }
        }
        return drew;
    }

    bool drawCubeDamage(
        float x, float y, float z,
        float qx, float qy, float qz, float qw,
        float sx, float sy, float sz,
        float *outPositions,
        unsigned int *outPositionsCount,
        float *outDamages)
    {
        unsigned int maxPositionsCount = *outPositionsCount;
        *outPositionsCount = 0;

        Matrix m(Vec{x, y, z}, Quat{qx, qy, qz, qw}, Vec{sx, sy, sz});

        bool drew = false;
        std::set<uint64_t> seenHashes;
        for (float dx = -1.f; dx <= 1.f; dx += 2.f)
        {
            for (float dz = -1.f; dz <= 1.f; dz += 2.f)
            {
                for (float dy = -1.f; dy <= 1.f; dy += 2.f)
                {
                    Vec p = (Vec(dx, dy, dz) * 0.5).applyMatrix(m);
                    float ax = p.x;
                    float ay = p.y;
                    float az = p.z;
                    vm::ivec3 min = vm::ivec3(std::floor(ax / (float)chunkSize), std::floor(ay / (float)chunkSize), std::floor(az / (float)chunkSize)) * chunkSize;
                    uint64_t minHash = hashOctreeMin(min);
                    if (seenHashes.find(minHash) == seenHashes.end())
                    {
                        seenHashes.insert(minHash);

                        Chunk &chunkNoise = getChunk(min);
                        if (chunkNoise.addCubeDamage(
                                x, y, z,
                                qx, qy, qz, qw,
                                sx, sy, sz))
                        {
                            if (*outPositionsCount < maxPositionsCount)
                            {
                                int gridSize = chunkSize + 3;
                                int damageBufferSize = gridSize * gridSize * gridSize;
                                memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.data(), sizeof(float) * damageBufferSize);

                                outPositions[(*outPositionsCount) * 3] = min.x;
                                outPositions[(*outPositionsCount) * 3 + 1] = min.y;
                                outPositions[(*outPositionsCount) * 3 + 2] = min.z;

                                (*outPositionsCount)++;
                            }

                            drew = true;
                        }
                    }
                }
            }
        }
        return drew;
    }

    bool eraseCubeDamage(
        float x, float y, float z,
        float qx, float qy, float qz, float qw,
        float sx, float sy, float sz,
        float *outPositions,
        unsigned int *outPositionsCount,
        float *outDamages)
    {
        unsigned int maxPositionsCount = *outPositionsCount;
        *outPositionsCount = 0;

        Matrix m(Vec{x, y, z}, Quat{qx, qy, qz, qw}, Vec{sx, sy, sz});

        bool drew = false;
        std::set<uint64_t> seenHashes;
        for (float dx = -1.f; dx <= 1.f; dx += 2.f)
        {
            for (float dz = -1.f; dz <= 1.f; dz += 2.f)
            {
                for (float dy = -1.f; dy <= 1.f; dy += 2.f)
                {
                    Vec p = (Vec(dx, dy, dz) * 0.5).applyMatrix(m);
                    float ax = p.x;
                    float ay = p.y;
                    float az = p.z;
                    vm::ivec3 min = vm::ivec3(std::floor(ax / (float)chunkSize), std::floor(ay / (float)chunkSize), std::floor(az / (float)chunkSize)) * chunkSize;
                    uint64_t minHash = hashOctreeMin(min);
                    if (seenHashes.find(minHash) == seenHashes.end())
                    {
                        seenHashes.insert(minHash);

                        Chunk &chunkNoise = getChunk(min);
                        if (chunkNoise.removeCubeDamage(
                                x, y, z,
                                qx, qy, qz, qw,
                                sx, sy, sz))
                        {
                            if (*outPositionsCount < maxPositionsCount)
                            {
                                int gridSize = chunkSize + 3;
                                int damageBufferSize = gridSize * gridSize * gridSize;
                                memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.data(), sizeof(float) * damageBufferSize);

                                outPositions[(*outPositionsCount) * 3] = min.x;
                                outPositions[(*outPositionsCount) * 3 + 1] = min.y;
                                outPositions[(*outPositionsCount) * 3 + 2] = min.z;

                                (*outPositionsCount)++;
                            }

                            drew = true;
                        }
                    }
                }
            }
        }
        return drew;
    }

    void injectDamage(const float &x, const float &y, const float &z, float *damageBuffer)
    {
        const vm::ivec3 min = vm::ivec3(x, y, z);
        uint64_t minHash = hashOctreeMin(min);
        const auto &iter = chunksNoiseHashMap.find(minHash);
        if (iter != chunksNoiseHashMap.end())
        {
            Chunk &chunkNoise = iter->second;
            chunkNoise.injectDamage(damageBuffer);
        }
        else
        {
            std::vector<float> cachedHeightField;

            int gridSize = chunkSize + 3;
            std::vector<float> cachedSdf(gridSize * gridSize * gridSize);
            memcpy(cachedSdf.data(), damageBuffer, cachedSdf.size() * sizeof(float));

            Chunk chunkNoise(min, std::move(cachedHeightField), std::move(cachedSdf));
            chunkNoise.initHeightField();

            chunksNoiseHashMap.emplace(std::make_pair(minHash, std::move(chunkNoise)));
        }
    }
}