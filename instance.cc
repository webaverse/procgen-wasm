#include "instance.h"
#include "procgen.h"
#include "octree.h"
// #include "lock.h"
#include "biomes.h"
#include "tracker.h"
#include "vector.h"
#include "util.h"
#include <emscripten.h>
// #include "peek.h"

constexpr int CHUNK_RANGE = 1;

// constructor/destructor
PGInstance::PGInstance(int seed, int chunkSize) :
    seed(seed),
    chunkSize(chunkSize),
    noises(seed)
    // cachedNoiseField(this),
    // cachedBiomesField(this),
    // cachedHeightField(this),
    // cachedWaterField(this),
    // cachedSkylightField(this),
    // cachedAoField(this),
    // cachedCaveField(this),
    // cachedSdf(this),
    // cachedWaterSdf(this)
    // cachedDamageSdf(this)
{
    // std::cout << "new pg instance " << seed << " " << chunkSize << std::endl;
}
PGInstance::~PGInstance() {}

// chunk

/* // locks
Mutex *PGInstance::getChunkLock(const vm::ivec2 &worldPos, const int lod, const int flags) {
    Mutex *chunkLock;
    uint64_t minLodHash = hashOctreeMin(worldPos);
    {
        std::unique_lock<Mutex> lock(locksMutex);
        chunkLock = &chunkLocks2D[minLodHash];
    }
    return chunkLock;
}
Mutex *PGInstance::getChunkLock(const vm::ivec3 &worldPos, const int lod) {
    Mutex *chunkLock;
    uint64_t minLodHash = hashOctreeMin(worldPos);
    {
        std::unique_lock<Mutex> lock(locksMutex);
        chunkLock = &chunkLocks3D[minLodHash];
    }
    return chunkLock;
}

// field ranges
void PGInstance::getHeightfieldRange(const vm::ivec2 &worldPositionXZ, const vm::ivec2 &size, int lod, float *heights) {
    
    for (int z = 0; z < size.y; z++)
    {
        for (int x = 0; x < size.x; x++)
        {
            int index2D = x + z * size.x;

            int ax = worldPositionXZ.x + x;
            int az = worldPositionXZ.y + z;
            heights[index2D] = cachedHeightField.get(ax, az).heightField;
        }
    }
}
void PGInstance::getLightRange(const vm::ivec3 &worldPosition, const vm::ivec3 &size, int lod, uint8_t *skylights, uint8_t *aos) {
    for (int z = 0; z < size.z; z++)
    {
        for (int y = 0; y < size.y; y++)
        {
            for (int x = 0; x < size.x; x++)
            {
                int dstIndex = x + y * size.x + z * size.x * size.y;

                int ax = worldPosition.x + x;
                int ay = worldPosition.y + y;
                int az = worldPosition.z + z;

                skylights[dstIndex] = cachedSkylightField.get(ax, ay, az);
                aos[dstIndex] = cachedAoField.get(ax, ay, az);
            }
        }
    }
} */

/* // fields
float *PGInstance::getChunkHeightfield(const vm::ivec2 &worldPositionXZ, int lod) {
    const int &size = chunkSize;
    
    float *heights = (float *)malloc(sizeof(float) * chunkSize * chunkSize);
    
    for (int z = 0; z < size; z++)
    {
        for (int x = 0; x < size; x++)
        {
            int index2D = x + z * size;

            int ax = worldPositionXZ.x + x;
            int az = worldPositionXZ.y + z;
            heights[index2D] = cachedHeightField.get(ax, az).heightField;
        }
    }
    
    return heights;
}
unsigned char *PGInstance::getChunkSkylight(const vm::ivec3 &worldPosition, int lod) {
    const int &size = chunkSize;

    unsigned char *skylights = (unsigned char *)malloc(sizeof(unsigned char) * size * size * size);

    for (int z = 0; z < size; z++)
    {
        for (int y = 0; y < size; y++)
        {
            for (int x = 0; x < size; x++)
            {
                int dstIndex = x + y * size + z * size * size;

                // int lx = x + 1;
                // int ly = y + 1;
                // int lz = z + 1;
                // int srcIndex = lx + lz * gridPoints + ly * gridPoints * gridPoints; // note: output is y-first, but storage is z-first

                int ax = worldPosition.x + x;
                int ay = worldPosition.y + y;
                int az = worldPosition.z + z;
                // int index = getIndex(ax, ay);
                skylights[dstIndex] = cachedSkylightField.get(ax, ay, az);
            }
        }
    }

    return skylights;
}
unsigned char *PGInstance::getChunkAo(const vm::ivec3 &worldPosition, int lod) {
    const int &size = chunkSize;

    unsigned char *aos = (unsigned char *)malloc(sizeof(unsigned char) * size * size * size);

    for (int z = 0; z < size; z++)
    {
        for (int y = 0; y < size; y++)
        {
            for (int x = 0; x < size; x++)
            {
                int dstIndex = x + y * size + z * size * size;
                // int srcIndex = x + z * size + y * size * size; // note: output is y-first, but storage is z-first

                int ax = worldPosition.x + x;
                int ay = worldPosition.y + y;
                int az = worldPosition.y + z;
                // int index = getIndex(ax, ay);
                aos[dstIndex] = cachedAoField.get(ax, ay, az);
            }
        }
    }

    return aos;
} */

// splats
uint8_t *PGInstance::createGrassSplat(const vm::ivec2 &worldPositionXZ, const int lod)
{
    std::vector<float> ps;
    std::vector<float> qs;
    std::vector<float> instances;
    unsigned int count = 0;

    int size = chunkSize * lod;

    // accumulate
    // Chunk2D &chunk = getChunk(worldPositionXZ, lod, GF_HEIGHTFIELD);
    int minX = worldPositionXZ.x / size * size;
    int minZ = worldPositionXZ.y / size * size;

    float seed = noises.grassNoise.in2D(minX, minZ);
    unsigned int seedInt;
    memcpy(&seedInt, &seed, sizeof(unsigned int));
    std::mt19937 rng(seedInt);

    const int maxNumGrasses = 4 * 1024;
    // ps.resize(maxNumGrasses * 3);
    // qs.resize(maxNumGrasses * 4);
    // instances.resize(maxNumGrasses);
    for (int i = 0; i < maxNumGrasses; i++)
    {
        float dx = (float)rng() / (float)0xFFFFFFFF * (float)size;
        float dz = (float)rng() / (float)0xFFFFFFFF * (float)size;

        float ax = (float)minX + dx;
        float az = (float)minZ + dz;

        // XXX this can be optimized to not fetch the full spec
        const Heightfield &heightfield = getHeightField(ax, az);
        const float &height = heightfield.heightField;

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
uint8_t *PGInstance::createVegetationSplat(const vm::ivec2 &worldPositionXZ, const int lod)
{
    std::vector<float> ps;
    std::vector<float> qs;
    std::vector<float> instances;
    unsigned int count = 0;

    // Chunk2D &chunk = getChunk(worldPositionXZ, lod, GF_HEIGHTFIELD);

    int minX = worldPositionXZ.x / chunkSize * chunkSize;
    int minZ = worldPositionXZ.y / chunkSize * chunkSize;

    float seed = noises.vegetationNoise.in2D(minX, minZ);
    unsigned int seedInt;
    memcpy(&seedInt, &seed, sizeof(unsigned int));
    std::mt19937 rng(seedInt);

    const int maxNumVeggies = 128;
    const float veggieRate = 0.35;
    for (int i = 0; i < maxNumVeggies; i++)
    {
        float dx = (float)rng() / (float)0xFFFFFFFF * (float)chunkSize;
        float dz = (float)rng() / (float)0xFFFFFFFF * (float)chunkSize;

        float ax = (float)minX + dx * lod;
        float az = (float)minZ + dz * lod;

        float noiseValue = noises.vegetationNoise.in2D(ax, az);

        if (noiseValue < veggieRate)
        {
            int index = 0;

            const Heightfield &heightfield = getHeightField(ax, az);
            const float &height = heightfield.heightField;

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
uint8_t *PGInstance::createMobSplat(const vm::ivec2 &worldPositionXZ, const int lod)
{
    std::vector<float> ps;
    std::vector<float> qs;
    std::vector<float> instances;
    unsigned int count = 0;

    // Chunk2D &chunk = getChunk(worldPositionXZ, lod, GF_HEIGHTFIELD);

    int minX = worldPositionXZ.x / chunkSize * chunkSize;
    int minZ = worldPositionXZ.y / chunkSize * chunkSize;

    float seed = noises.mobNoise.in2D(minX, minZ);
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

        float noiseValue = noises.mobNoise.in2D(ax, az);

        if (noiseValue < mobRate)
        {
            const Heightfield &heightfield = getHeightField(ax, az);
            const float &height = heightfield.heightField;

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

/* std::vector<vm::ivec3> getChunkRangeInclusive(const vm::ivec3 &worldPosition, int minChunkDelta, int maxChunkDelta, int chunkSize) {
    std::vector<vm::ivec3> result;
    for (int dy = minChunkDelta; dy <= maxChunkDelta; dy++)
    {
        for (int dz = minChunkDelta; dz <= maxChunkDelta; dz++)
        {
            for (int dx = minChunkDelta; dx <= maxChunkDelta; dx++)
            {
                result.push_back(vm::ivec3{
                    worldPosition.x + dx * chunkSize,
                    worldPosition.y + dy * chunkSize,
                    worldPosition.z + dz * chunkSize
                });
            }
        }
    }
    return result;
}
std::vector<vm::ivec2> getChunkRangeInclusive(const vm::ivec2 &worldPosition, int minChunkDelta, int maxChunkDelta, int chunkSize) {
    std::vector<vm::ivec2> result;
    for (int dz = minChunkDelta; dz <= maxChunkDelta; dz++)
    {
        for (int dx = minChunkDelta; dx <= maxChunkDelta; dx++)
        {
            result.push_back(vm::ivec2{
                worldPosition.x + dx * chunkSize,
                worldPosition.y + dz * chunkSize
            });
        }
    }
    return result;
} */

//

class Geometry {
public:
    std::vector<float> positions;
    std::vector<float> normals;
    std::vector<float> uvs;
    std::vector<uint32_t> indices;

    // Geometry() {}
};
Geometry createPlaneGeometry(int width, int height, int widthSegments, int heightSegments, const std::vector<Heightfield> &heightfield) {
    Geometry geometry;

    /* this.parameters = {
        width: width,
        height: height,
        widthSegments: widthSegments,
        heightSegments: heightSegments
    }; */

    // const int width_half = width / 2;
    // const int height_half = height / 2;

    const int gridX = widthSegments;
    const int gridY = heightSegments;

    const int gridX1 = gridX + 1;
    const int gridY1 = gridY + 1;

    const int segment_width = width / gridX;
    const int segment_height = height / gridY;

    //

    // const indices = [];
    // const vertices = [];
    // const normals = [];
    // const uvs = [];

    int index = 0;
    for (int iy = 0; iy < gridY1; iy ++) {

        // const int y = iy * segment_height - height_half;
        const int y = iy * segment_height;

        for (int ix = 0; ix < gridX1; ix++) {

            // const int x = ix * segment_width - width_half;
            const int x = ix * segment_width;

            /* if (index >= heightfield.size()) {
                std::cerr << "height overflow " << index << " " << heightfield.size() << " " << gridX1 << " " << gridY1 << std::endl;
                abort();
            } */
            const Heightfield &localHeightfield = heightfield[index];
            const float &height = localHeightfield.heightField;

            geometry.positions.push_back(x);
            geometry.positions.push_back(height);
            geometry.positions.push_back(-y);

            geometry.normals.push_back(0);
            geometry.normals.push_back(1);
            geometry.normals.push_back(0);

            geometry.uvs.push_back((float)ix / (float)gridX);
            geometry.uvs.push_back(1.f - ((float)iy / (float)gridY));

            index++;

        }

    }

    for (int iy = 0; iy < gridY; iy++) {

        for (int ix = 0; ix < gridX; ix++) {

            const int a = ix + gridX1 * iy;
            const int b = ix + gridX1 * (iy + 1);
            const int c = (ix + 1) + gridX1 * (iy + 1);
            const int d = (ix + 1) + gridX1 * iy;

            // geometry.indices.push_back(a);
            // geometry.indices.push_back(b);
            // geometry.indices.push_back(d);
            // geometry.indices.push_back(b);
            // geometry.indices.push_back(c);
            // geometry.indices.push_back(d);

            geometry.indices.push_back(a);
            geometry.indices.push_back(d);
            geometry.indices.push_back(b);
            geometry.indices.push_back(b);
            geometry.indices.push_back(d);
            geometry.indices.push_back(c);

        }

    }

    // this.setIndex( indices );
    // this.setAttribute( 'position', new Float32BufferAttribute( vertices, 3 ) );
    // this.setAttribute( 'normal', new Float32BufferAttribute( normals, 3 ) );
    // this.setAttribute( 'uv', new Float32BufferAttribute( uvs, 2 ) );

    return geometry;
}
void generateHeightfieldMesh(const vm::ivec2 &worldPosition, int lod, int chunkSize, const std::vector<Heightfield> &heightfield, TerrainVertexBuffer &vertexBuffer) {
    // std::cout << "generate plane geometry: " << worldPosition.x << " " << worldPosition.y << " " << lod << " " << chunkSize << std::endl;
    const int worldSize = chunkSize * lod;
    Geometry planeGeometry = createPlaneGeometry(worldSize, worldSize, chunkSize, chunkSize, heightfield);

    for (size_t i = 0, j = 0; i < planeGeometry.positions.size(); i += 3, j += 2) { 
      vertexBuffer.positions.push_back(vm::vec3{
        planeGeometry.positions[i] + (float)worldPosition.x,
        planeGeometry.positions[i + 1] - (float)WORLD_BASE_HEIGHT,
        planeGeometry.positions[i + 2] + (float)worldPosition.y
      });

      vertexBuffer.normals.push_back(vm::vec3{
        planeGeometry.normals[i],
        planeGeometry.normals[i + 1],
        planeGeometry.normals[i + 2]
      });

      // size_t j = i / 3 * 2;
      vertexBuffer.biomesUvs1.push_back(std::array<UV, 2>{
        planeGeometry.uvs[j],
        planeGeometry.uvs[j + 1]
      });
    }
    for (size_t i = 0; i < planeGeometry.indices.size(); i++) {
      vertexBuffer.indices.push_back(planeGeometry.indices[i]);
    }
}
uint8_t *PGInstance::createTerrainChunkMesh(const vm::ivec2 &worldPosition, int lod) {
    const int chunkSizeP1 = chunkSize + 1;
    std::vector<Heightfield> heightfields(chunkSizeP1 * chunkSizeP1);
    getHeightField(worldPosition.x, worldPosition.y, lod, heightfields.data());

    // TerrainDCContext vertexContext;
    // generateMeshFromOctree<TerrainDCContext, false>(chunkOctree.root, vertexContext);
    // generateMeshFromOctree<TerrainDCContext, true>(chunkOctree.seamRoot, vertexContext);

    // auto &vertexBuffer = vertexContext.vertexBuffer;
    TerrainVertexBuffer vertexBuffer;
    generateHeightfieldMesh(worldPosition, lod, chunkSize, heightfields, vertexBuffer);
    
    /* if (vertexBuffer.indices.size() == 0)
    {
        // printf("Generated Mesh Is Not Valid\n");
        return nullptr;
    } */

    // const vm::ivec3 chunkMax = worldPosition + (chunkSize * lod);
    // setPeeks<TerrainDCContext>(this, worldPosition, chunkMax, lod, vertexBuffer.peeks, PEEK_FACE_INDICES.array);

return vertexBuffer.getBuffer();
}
uint8_t *PGInstance::createLiquidChunkMesh(const vm::ivec2 &worldPosition, int lod)
{
    const int chunkSizeP1 = chunkSize + 1;
    std::vector<float> waterfields(chunkSizeP1 * chunkSizeP1);
    getWaterField(worldPosition.x, worldPosition.y, lod, waterfields.data());

    // ChunkOctree<LiquidDCContext> chunkOctree(this, worldPosition, lodArray);
    // if (!chunkOctree.root)
    // {
    //     // printf("Chunk Has No Data\n");
    //     return nullptr;
    // }
    // WaterDCContext vertexContext;
    // generateMeshFromOctree<LiquidDCContext, false>(chunkOctree.root, vertexContext);
    // generateMeshFromOctree<LiquidDCContext, true>(chunkOctree.seamRoot, vertexContext);

    // auto &vertexBuffer = vertexContext.vertexBuffer;
    LiquidVertexBuffer vertexBuffer;
    if (vertexBuffer.indices.size() == 0)
    {
        // printf("Generated Mesh Is Not Valid\n");
        return nullptr;
    }

    return vertexBuffer.getBuffer();
}

//

/* bool PGInstance::drawSphereDamage(const float &x, const float &y, const float &z,
                                  const float &radius, float *outPositions, unsigned int *outPositionsCount,
                                  const int &lod)
{
    return damageBuffers.damage(vm::vec3{x,y,z}, radius, outPositions, outPositionsCount, lod);
}

bool PGInstance::eraseSphereDamage(const float &x, const float &y, const float &z,
                                   const float radius, float *outPositions, unsigned int *outPositionsCount, float *outDamages,
                                   const int &lod)
{
    unsigned int maxPositionsCount = *outPositionsCount;
    *outPositionsCount = 0;

    bool drew = false;
    std::set<uint64_t> seenHashes;

    // chunk min of the hit point
    vm::ivec3 chunkMin = chunkMinForPosition(vm::ivec3{(int)x, (int)y, (int)z}, lod);

    for (float dx = -1; dx <= 1; dx += 1)
    {
        for (float dz = -1; dz <= 1; dz += 1)
        {
            for (float dy = -1; dy <= 1; dy += 1)
            {
               vm::ivec3 min = chunkMin + vm::ivec3{(int)dx, (int)dy, (int)dz} * chunkSize;

                uint64_t minHash = hashOctreeMin(min);
                if (seenHashes.find(minHash) == seenHashes.end())
                {
                    seenHashes.insert(minHash);

                    // Chunk3D &chunkNoise = getChunk(min, lod, GF_SDF);
                    if (removeSphereDamage(x, y, z, radius))
                    {
                        if (*outPositionsCount < maxPositionsCount)
                        {
                            // int gridSize = chunkSize + 3 + lod;
                            // int damageBufferSize = gridSize * gridSize * gridSize;
                            // memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.value.data(), sizeof(float) * damageBufferSize);

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
} */

/* bool PGInstance::drawCubeDamage(
    float x, float y, float z,
    float qx, float qy, float qz, float qw,
    float sx, float sy, float sz,
    float *outPositions,
    unsigned int *outPositionsCount,
    float *outDamages,
    const int &lod)
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
                vm::ivec3 min = vm::ivec3{
                                    (int)std::floor(ax / (float)chunkSize),
                                    (int)std::floor(ay / (float)chunkSize),
                                    (int)std::floor(az / (float)chunkSize)} *
                                chunkSize;
                uint64_t minHash = hashOctreeMin(min);
                if (seenHashes.find(minHash) == seenHashes.end())
                {
                    seenHashes.insert(minHash);

                    // Chunk3D &chunkNoise = getChunk(min, lod, GF_SDF);
                    if (addCubeDamage(
                            x, y, z,
                            qx, qy, qz, qw,
                            sx, sy, sz))
                    {
                        if (*outPositionsCount < maxPositionsCount)
                        {
                            // int gridSize = chunkSize + 3 + lod;
                            // int damageBufferSize = gridSize * gridSize * gridSize;
                            // memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.value.data(), sizeof(float) * damageBufferSize);

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

bool PGInstance::eraseCubeDamage(
    float x, float y, float z,
    float qx, float qy, float qz, float qw,
    float sx, float sy, float sz,
    float *outPositions,
    unsigned int *outPositionsCount,
    float *outDamages,
    const int &lod)
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
                vm::ivec3 min = vm::ivec3{
                                    (int)std::floor(ax / (float)chunkSize),
                                    (int)std::floor(ay / (float)chunkSize),
                                    (int)std::floor(az / (float)chunkSize)} *
                                chunkSize;
                uint64_t minHash = hashOctreeMin(min);
                if (seenHashes.find(minHash) == seenHashes.end())
                {
                    seenHashes.insert(minHash);

                    // Chunk3D &chunkNoise = getChunk(min, lod, GF_SDF);
                    if (removeCubeDamage(
                            x, y, z,
                            qx, qy, qz, qw,
                            sx, sy, sz))
                    {
                        if (*outPositionsCount < maxPositionsCount)
                        {
                            // int gridSize = chunkSize + 3 + lod;
                            // int damageBufferSize = gridSize * gridSize * gridSize;
                            // memcpy(outDamages + (*outPositionsCount) * damageBufferSize, chunkNoise.cachedSdf.value.data(), sizeof(float) * damageBufferSize);

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
} */

void PGInstance::setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix) {
    this->worldPosition = worldPosition;
    this->cameraPosition = cameraPosition;
    this->cameraQuaternion = cameraQuaternion;
    this->projectionMatrix = projectionMatrix;
    
    ProcGen::taskQueue.setCamera(worldPosition, cameraPosition, cameraQuaternion, projectionMatrix);
}

void PGInstance::setClipRange(const vm::vec2 &min, const vm::vec2 &max)
{
    clipRange.reset(new vm::box3{
        vm::vec3{min.x, min.y},
        vm::vec3{max.x, max.y}});
}

//

void PGInstance::createTerrainChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod)
{
    // std::cout << "create terrain async " << worldPosition.x << " " << worldPosition.y << std::endl;
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // int lod = lodArray[0];
    // std::vector<int> lodVector(lodArray, lodArray + 8);

    vm::vec3 worldPositionF{
        (float)worldPosition.x + (float)lod / 2.f,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y + (float)lod / 2.f
    };
    // vm::vec3 sizeF{
    //     (float)lod / 2.f,
    //     (float)lod / 2.f,
    //     (float)lod / 2.f
    // };
    Task *terrainTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod
    ]() -> void {
        uint8_t *result = createTerrainChunkMesh(worldPosition, lod);
        if (!promise->resolve(result)) {
            // XXX cleanup
            // XXX also cleanup the other async calls
        }
    });
    ProcGen::taskQueue.pushTask(terrainTask);
}
void PGInstance::createLiquidChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod)
{
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // int lod = lodArray[0];
    // std::vector<int> lodVector(lodArray, lodArray + 8);

    vm::vec3 worldPositionF{
        (float)worldPosition.x + (float)lod / 2.f,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y + (float)lod / 2.f
    };
    // vm::vec3 sizeF{
    //     (float)lod / 2.f,
    //     (float)lod / 2.f,
    //     (float)lod / 2.f
    // };
    Task *liquidTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod
    ]() -> void {
        uint8_t *result = createLiquidChunkMesh(worldPosition, lod);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(liquidTask);
}

// get ranges
/* void PGInstance::getHeightfieldRangeAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const vm::ivec2 &sizeXZ, int lod, float *heights, int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x + (float)sizeXZ.x / 2.f,
        0.f,
        (float)worldPositionXZ.y + (float)sizeXZ.y / 2.f
    };
    vm::vec3 sizeF{
        (float)sizeXZ.x / 2.f,
        0.f,
        (float)sizeXZ.y / 2.f
    };
    Task *heightfieldRangeTask = new Task(id, worldPositionF, sizeF, priority, [
        this,
        promise,
        worldPositionXZ,
        sizeXZ,
        lod,
        heights
    ]() -> void {
        getHeightfieldRange(worldPositionXZ, sizeXZ, lod, heights);
        void *result = nullptr;
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(heightfieldRangeTask);
}
void PGInstance::getLightRangeAsync(uint32_t id, const vm::ivec3 &worldPosition, const vm::ivec3 &size, int lod, uint8_t *skylights, uint8_t *aos, int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x + size.x / 2.f,
        (float)worldPosition.y + size.y / 2.f,
        (float)worldPosition.z + size.z / 2.f
    };
    vm::vec3 worldSizeF{
        (float)size.x / 2.f,
        (float)size.y / 2.f,
        (float)size.z / 2.f
    };
    Task *lightRangeTask = new Task(id, worldPositionF, worldSizeF, priority, [
        this,
        promise,
        worldPosition,
        size,
        lod,
        skylights,
        aos
    ]() -> void {
        getLightRange(worldPosition, size, lod, skylights, aos);
        void *result = nullptr;
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(lightRangeTask);
} */

/* // get chunk attributes
void PGInstance::getChunkHeightfieldAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, int lod, int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // std::cout << "get chunk heightfield priority " << priority << std::endl;

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x,
        0.f,
        (float)worldPositionXZ.y
    };
    Task *heightfieldTask = new Task(id, worldPositionF, lod, priority, [
        this,
        promise,
        worldPositionXZ,
        lod
    ]() -> void {
        void *result = getChunkHeightfield(worldPositionXZ, lod);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(heightfieldTask);
}
void PGInstance::getChunkSkylightAsync(uint32_t id, const vm::ivec3 &worldPosition, int lod, int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x,
        (float)worldPosition.y,
        (float)worldPosition.z
    };
    Task *skylightTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod
    ]() -> void {
        void *result = getChunkSkylight(worldPosition, lod);
        promise->resolve(result); });
    ProcGen::taskQueue.pushTask(skylightTask);
}
void PGInstance::getChunkAoAsync(uint32_t id, const vm::ivec3 &worldPosition, int lod, int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);
    
    vm::vec3 worldPositionF{
        (float)worldPosition.x,
        (float)worldPosition.y,
        (float)worldPosition.z
    };
    Task *aoTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod
    ]() -> void {
        void *result = getChunkAo(worldPosition, lod);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(aoTask);
} */

void PGInstance::createGrassSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x + (float)lod / 2.f,
        0.f,
        (float)worldPositionXZ.y + (float)lod / 2.f
    };
    // vm::vec3 sizeF{
    //     (float)lod / 2.f,
    //     0.f,
    //     (float)lod / 2.f
    // };
    /* std::cout << "grass splat world position" <<
        worldPositionF.x << " " <<
        worldPositionF.y << " " <<
        worldPositionF.z << std::endl; */
    Task *grassSplatTask = new Task(id, worldPositionF, lod, priority, [
        this,
        promise,
        worldPositionXZ,
        lod
    ]() -> void {
        uint8_t *result = createGrassSplat(worldPositionXZ, lod);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(grassSplatTask);
}
void PGInstance::createVegetationSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x + (float)lod / 2.f,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPositionXZ.y + (float)lod / 2.f
    };
    // vm::vec3 sizeF{
    //     (float)lod / 2.f,
    //     0.f,
    //     (float)lod / 2.f
    // };
    Task *vegetationSplatTask = new Task(id, worldPositionF, lod, priority, [
        this,
        promise,
        worldPositionXZ,
        lod
    ]() -> void {
        void *result = createVegetationSplat(worldPositionXZ, lod);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(vegetationSplatTask);
}
void PGInstance::createMobSplatAsync(uint32_t id, const vm::ivec2 &worldPositionXZ, const int lod, const int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x + (float)lod / 2.f,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPositionXZ.y + (float)lod / 2.f
    };
    Task *mobSplatTask = new Task(id, worldPositionF, lod, priority, [
        this,
        promise,
        worldPositionXZ,
        lod
    ]() -> void {
        void *result = createMobSplat(worldPositionXZ, lod);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(mobSplatTask);
}

// 2d caches

NoiseField PGInstance::getNoise(int bx, int bz) {
    // const int &size = chunkSize;
    // const vm::ivec2 &min = chunk->min;
    // const int &lod = chunk->lod;
    
    /* NoiseField noiseField;
    noiseField.temperature.resize(size * size);
    noiseField.humidity.resize(size * size);
    noiseField.ocean.resize(size * size);
    noiseField.river.resize(size * size);
    for (int z = 0; z < size; z++)
    {
        for (int x = 0; x < size; x++)
        { */
            // int index = x + z * size;
            // int ax = x;
            // int az = y;

            float tNoise = (float)noises.temperatureNoise.in2D(bx, bz);
            // noiseField.temperature[index] = tNoise;

            float hNoise = (float)noises.humidityNoise.in2D(bx, bz);
            // noiseField.humidity[index] = hNoise;

            float oNoise = (float)noises.oceanNoise.in2D(bx, bz);
            // noiseField.ocean[index] = oNoise;

            float rNoise = (float)noises.riverNoise.in2D(bx, bz);
            // noiseField.river[index] = rNoise;

            return NoiseField{
                tNoise,
                hNoise,
                oNoise,
                rNoise
            };
        /* }
    }

    return noiseField; */
}
uint8_t PGInstance::getBiome(int bx, int bz) {
    // const int &size = chunkSize;
    // const auto &cachedNoiseField = chunk->cachedNoiseField;
    
    // std::vector<uint8_t> biomesField(size * size);
    /* for (int z = 0; z < size; z++)
    {
        for (int x = 0; x < size; x++)
        { */
            // int index = x + z * size;
            // unsigned char biome = cachedBiomesField.get(x, z);
            unsigned char biome = 0xFF;

            const auto &noise = getNoise(bx, bz);
            float temperatureNoise = noise.temperature;
            float humidityNoise = noise.humidity;
            float oceanNoise = noise.ocean;
            float riverNoise = noise.river;

            if (oceanNoise < (80.0f / 255.0f))
            {
                biome = (unsigned char)BIOME::biOcean;
            }
            if (biome == 0xFF)
            {
                const float range = 0.022f;
                if (riverNoise > 0.5f - range && riverNoise < 0.5f + range)
                {
                    biome = (unsigned char)BIOME::biRiver;
                }
            }
            if (std::pow(temperatureNoise, 1.3f) < ((4.0f * 16.0f) / 255.0f))
            {
                if (biome == (unsigned char)BIOME::biOcean)
                {
                    biome = (unsigned char)BIOME::biFrozenOcean;
                }
                else if (biome == (unsigned char)BIOME::biRiver)
                {
                    biome = (unsigned char)BIOME::biFrozenRiver;
                }
            }
            if (biome == 0xFF)
            {
                float temperatureNoise2 = vm::clamp(std::pow(temperatureNoise, 1.3f), 0.f, 1.f);
                float humidityNoise2 = vm::clamp(std::pow(humidityNoise, 1.3f), 0.f, 1.f);

                int t = (int)std::floor(temperatureNoise2 * 16.0f);
                int h = (int)std::floor(humidityNoise2 * 16.0f);
                biome = (unsigned char)BIOMES_TEMPERATURE_HUMIDITY[t + 16 * h];
            }
            return biome;
        /* }
    } */
    // return biomesField;
}
void PGInstance::getHeightField(int bx, int bz, int lod, Heightfield *heightfield) {
    const int chunkSizeP1 = chunkSize + 1;
    for (int z = 0; z < chunkSizeP1; z++)
    {
        for (int x = 0; x < chunkSizeP1; x++)
        {
            int index2D = x + z * chunkSizeP1;
            Heightfield &localHeightfield = heightfield[index2D];
            // int ax = x;
            // int az = z;
            
            // int lx = x - 1;
            // int lz = z - 1;
            // int index2D2 = lx + lz * size;
            // bool isInRange = lx >= 0 && lx < size && lz >= 0 && lz < size;

            localHeightfield = getHeightField(bx + x * lod, bz + z * lod);
        }
    }
    // return heightfield;
}
Heightfield PGInstance::getHeightField(int bx, int bz) {
    Heightfield localHeightfield;

    std::unordered_map<unsigned char, unsigned int> biomeCounts(numBiomes);
    int numSamples = 0;
    for (int dz = -chunkSize/2; dz < chunkSize/2; dz++)
    {
        for (int dx = -chunkSize/2; dx < chunkSize/2; dx++)
        {
            int ax = bx + dx;
            int az = bz + dz;
            unsigned char b = getBiome(ax, az);

            biomeCounts[b]++;
            numSamples++;
        }
    }

    std::vector<unsigned char> seenBiomes;
    for (auto kv : biomeCounts)
    {
        unsigned char b = kv.first;
        seenBiomes.push_back(b);
    }
    
    // sort by increasing occurence count of the biome
    std::sort(
        seenBiomes.begin(),
        seenBiomes.end(),
        [&](unsigned char b1, unsigned char b2) -> bool {
            return biomeCounts[b1] > biomeCounts[b2];
        }
    );

    for (size_t i = 0; i < 4; i++)
    {
        if (i < seenBiomes.size())
        {
            localHeightfield.biomesVectorField[i] = seenBiomes[i];
            localHeightfield.biomesWeightsVectorField[i] = (float)biomeCounts[seenBiomes[i]] / (float)numSamples * 255.0f;
        }
        else
        {
            localHeightfield.biomesVectorField[i] = 0;
            localHeightfield.biomesWeightsVectorField[i] = 0;
        }
    }

    float elevationSum = 0.f;
    vm::vec2 fWorldPosition{(float)bx, (float)bz};
    for (auto const &iter : biomeCounts)
    {
        elevationSum += iter.second * getComputedBiomeHeight(iter.first, fWorldPosition);
    }

    float elevation = elevationSum / (float)numSamples;
    localHeightfield.heightField = elevation;

    return localHeightfield;
}
void PGInstance::getWaterField(int bx, int bz, int lod, float *waterfield) {
    const int chunkSizeP1 = chunkSize + 1;
    for (int z = 0; z < chunkSizeP1; z++)
    {
        for (int x = 0; x < chunkSizeP1; x++)
        {
            int index2D = x + z * chunkSizeP1;
            float &localWaterfield = waterfield[index2D];

            float value = 0;
            for (int dz = -chunkSize/2; dz < chunkSize/2; dz++)
            {
                for (int dx = -chunkSize/2; dx < chunkSize/2; dx++)
                {
                    int ax = bx + dx;
                    int az = bz + dz;
                    unsigned char b = getBiome(ax, az);

                    if (isWaterBiome(b)) {
                        value++;
                    }
                }
            }
            localWaterfield = value;
        }
    }
}
float PGInstance::getWaterField(int bx, int bz) {
    float value = 0;
    for (int dz = -chunkSize/2; dz < chunkSize/2; dz++)
    {
        for (int dx = -chunkSize/2; dx < chunkSize/2; dx++)
        {
            int ax = bx + dx;
            int az = bz + dz;
            unsigned char b = getBiome(ax, az);

            if (isWaterBiome(b)) {
                value++;
            }
        }
    }
    return value;
}

// 3d caches

/* uint8_t PGInstance::initSkylightField(PGInstance *inst, int x, int y, int z) {
    constexpr uint8_t maxSkyLightu8 = 8;
    constexpr int maxSkyLighti = (int)maxSkyLightu8;
    constexpr float maxSkyLightf = (float)maxSkyLighti;

    skylight = std::min(std::max(skylight, (uint8_t)0), maxSkyLightu8);

    // XXX should flood fill the light

    return skylight;
}
uint8_t PGInstance::initAoField(PGInstance *inst, int x, int y, int z) {
    const int &size = chunkSize;

    unsigned char numOpens = 0;
    for (int dy = -1; dy <= 1; dy++)
    {
        for (int dz = -1; dz <= 1; dz++)
        {
            for (int dx = -1; dx <= 1; dx++)
            {
                // int sdfIndex = (lx + dx) + (lz + dz) * gridPoints + (ly + dy) * gridPoints * gridPoints;
                numOpens += (unsigned char)(inst->cachedSdf.get(x + dx, y + dy, z + dz) >= 0.f);
            }
        }
    }

    return numOpens;
} */
/* float PGInstance::initCaveField(PGInstance *inst, int x, int y, int z) {
    return ProcGen::getComputedCaveNoise(x, y, z);
} */
float randomFromPoint(int x, int y, int z) {
    uint64_t hash = hashOctreeMin(vm::ivec3{x, y, z});
    uint32_t hash32 = (uint32_t)hash ^ (uint32_t)(hash >> 32);
    float f = (float)hash32 / (float)UINT32_MAX;
    return f;
}
/* float PGInstance::initSdf(PGInstance *inst, int x, int y, int z) {
    float height = inst->cachedHeightField.get(x, z).heightField;

    float heightValue = (float)y - height;
    heightValue = std::min(
        std::max(
            heightValue,
            (float)-1),
        (float)1);

    float caveValue = inst->cachedCaveField.get(x, y, z);
    float f = heightValue + caveValue * 1.1f + randomFromPoint(x, y, z) * 0.0001f;

    return f;
} */
/* float PGInstance::initWaterSdf(PGInstance *inst, int x, int y, int z) {
    const float fSize = (float)chunkSize;

    float waterValue = -inst->cachedWaterField.get(x, z) / fSize;

    float heightValue = (float)y - waterBaseHeight;
    heightValue = std::min(
        std::max(
            heightValue,
            -1.f
        ),
        1.f
    );
    
    float value = std::max(waterValue, heightValue);

    return value;
} */

//

// biomes
/* void PGInstance::getCachedBiome2D(const vm::ivec2 &worldPosition, vm::ivec4 &biome, vm::vec4 &biomeWeights, std::array<UV, 2> &biomeUvs1, std::array<UV, 2> &biomeUvs2) {
    const auto &heightfield = cachedHeightField.get(worldPosition.x, worldPosition.y);
    biome.x = heightfield.biomesVectorField[0];
    biome.y = heightfield.biomesVectorField[1];
    biome.z = heightfield.biomesVectorField[2];
    biome.w = heightfield.biomesVectorField[3];

    biomeWeights.x = heightfield.biomesWeightsVectorField[0];
    biomeWeights.y = heightfield.biomesWeightsVectorField[1];
    biomeWeights.z = heightfield.biomesWeightsVectorField[2];
    biomeWeights.w = heightfield.biomesWeightsVectorField[3];

    biomeUvs1[0] = BIOME_UVS[(int)biome.x];
    biomeUvs1[1] = BIOME_UVS[(int)biome.y];
    biomeUvs2[0] = BIOME_UVS[(int)biome.z];
    biomeUvs2[1] = BIOME_UVS[(int)biome.w];
}
inline void shiftOverrideBiome(vm::ivec4 &biome, vm::vec4 &biomeWeights, std::array<UV, 2> &biomeUvs1, std::array<UV, 2> &biomeUvs2, BIOME b) {
    // move the biomes to make room
    biome.w = biome.z;
    biome.z = biome.y;
    biome.y = biome.x;
    biome.x = (unsigned char)b;

    biomeWeights.x = 1;
    biomeWeights.y = 0;
    biomeWeights.z = 0;
    biomeWeights.w = 0;

    biomeUvs1[0] = BIOME_UVS[(int)b];
    biomeUvs1[1] = {0, 0};
    biomeUvs2[0] = {0, 0};
    biomeUvs2[1] = {0, 0};
}
void PGInstance::getCachedInterpolatedBiome3D(const vm::vec3 &worldPosition, vm::ivec4 &biome, vm::vec4 &biomeWeights, std::array<UV, 2> &biomeUvs1, std::array<UV, 2> &biomeUvs2) {
    vm::ivec3 iWorldPosition{(int)worldPosition.x, (int)worldPosition.y, (int)worldPosition.z};
    vm::ivec2 iWorldPositionXZ{(int)worldPosition.x, (int)worldPosition.z};

    const int &x = iWorldPosition.x;
    const int &y = iWorldPosition.y;
    const int &z = iWorldPosition.z;

    getCachedBiome2D(iWorldPositionXZ, biome, biomeWeights, biomeUvs1, biomeUvs2);

    float heightValue = cachedHeightField.get(x, z).heightField;
    float sdfValue = cachedSdf.get(x, y, z);

    bool neighborHeightsValid = true;
    for (int dx = -1; dx <= 1; dx += 2)
    {
        for (int dz = -1; dz <= 1; dz += 2)
        {
            int lx2 = x + dx;
            int lz2 = z + dz;
            float heightValue = cachedHeightField.get(lx2, lz2).heightField;
            if (y + 3 > heightValue)
            {
                neighborHeightsValid = false;
                break;
            }
        }
        if (!neighborHeightsValid)
        {
            break;
        }
    }

    if (neighborHeightsValid)
    {
        if (y < heightValue - 12)
        {
            shiftOverrideBiome(biome, biomeWeights, biomeUvs1, biomeUvs2, BIOME::teStone);
        }
        else if (y < heightValue - 2)
        {
            shiftOverrideBiome(biome, biomeWeights, biomeUvs1, biomeUvs2, BIOME::teDirt);
        }
    }
} */

// lighting
/* void PGInstance::getCachedInterpolatedLight(const vm::vec3 &worldPosition, uint8_t &skylight, uint8_t &ao) {
    vm::ivec3 iWorldPosition{(int)worldPosition.x, (int)worldPosition.y, (int)worldPosition.z};

    const int &x = iWorldPosition.x;
    const int &y = iWorldPosition.y;
    const int &z = iWorldPosition.z;

    skylight = cachedSkylightField.get(x, y, z);
    ao = cachedAoField.get(x, y, z);
} */

// heightfield

/* float heightfieldLens(const Heightfield &heightfield) {
    return heightfield.heightField;
}
float PGInstance::getCachedInterpolatedHeightfield(const vm::vec2 &worldPosition, const int lod) {
    return bilinearMap<
        decltype(cachedHeightField),
        float,
        decltype(heightfieldLens)
    >(
        worldPosition,
        lod,
        cachedHeightField,
        heightfieldLens
    );
} */

//

// signed distance field function for a box at the origin
// returns negative for points inside the box, zero at the box's surface, and positive for points outside the box
// sx sy sz is the size of the box. the box goes from -sx/2 to sx/2, -sy/2 to sy/2, -sz/2 to sz/2
// px py pz is the point to check
/* float PGInstance::signedDistanceToBox(float sx, float sy, float sz, float px, float py, float pz)
{
    float dx = std::abs(px) - sx / 2;
    float dy = std::abs(py) - sy / 2;
    float dz = std::abs(pz) - sz / 2;
    float d = std::max(std::max(dx, dy), dz);
    return d;
}

// signed distance to sphere
// returns negative for points inside the sphere, zero at the sphere's surface, and positive for points outside the sphere
// cx, cy, cz is the center of the sphere. r is the radius. px, py, pz is the point to check
float PGInstance::signedDistanceToSphere(float cx, float cy, float cz, float r, float px, float py, float pz)
{
    float dx = px - cx;
    float dy = py - cy;
    float dz = pz - cz;
    float d = sqrt(dx * dx + dy * dy + dz * dz);
    return d - r;
} */

// biomes
float PGInstance::getComputedBiomeHeight(unsigned char b, const vm::vec2 &worldPosition) {
    const Biome &biome = BIOMES[b];
    float ax = worldPosition.x;
    float az = worldPosition.y;

    float biomeHeight = biome.baseHeight +
        noises.elevationNoise1.in2D(ax * biome.amps[0][0], az * biome.amps[0][0]) * biome.amps[0][1] +
        noises.elevationNoise2.in2D(ax * biome.amps[1][0], az * biome.amps[1][0]) * biome.amps[1][1] +
        noises.elevationNoise3.in2D(ax * biome.amps[2][0], az * biome.amps[2][0]) * biome.amps[2][1];
    return biomeHeight;
}

void PGInstance::trackerUpdateAsync(uint32_t id, Tracker *tracker, const vm::vec3 &position, int priority) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // std::cout << "tracker update async " << priority << std::endl;

    Task *trackerUpdateTask = new Task(id, priority, [
        this,
        promise,
        tracker,
        position
    ]() -> void {
        const TrackerUpdate &trackerUpdate = tracker->update(position);
        uint8_t *buffer = trackerUpdate.getBuffer();
        // std::cout << "trakcer update buffer address" << (void *)buffer << std::endl;
        if (!promise->resolve(buffer)) {
          // XXX clean up
        }
    });
    ProcGen::taskQueue.pushTask(trackerUpdateTask);
}