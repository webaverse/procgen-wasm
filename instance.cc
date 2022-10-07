#include "instance.h"
#include "procgen.h"
#include "octree.h"
// #include "lock.h"
#include "biomes.h"
#include "tracker.h"
#include "vector.h"
#include "util.h"
#include "MurmurHash3.h"
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
            const float height = getHeight(ax, az);

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

/* class Geometry {
public:
    std::vector<float> positions;
    std::vector<float> normals;
    // std::vector<float> uvs;
    std::vector<uint32_t> indices;

    // Geometry() {}
}; */
void normalizeNormals(std::vector<vm::vec3> &normals) {
    for (size_t i = 0, il = normals.size(); i < il; i++) {
        Vec vec{
            normals[i].x,
            normals[i].y,
            normals[i].z
        };
        vec.normalize();
        // normals.setXYZ(i, vec.x, vec.y, vec.z);
        normals[i].x = vec.x;
        normals[i].y = vec.y;
        normals[i].z = vec.z;
    }
}
void computeVertexNormals(std::vector<vm::vec3> &positions, std::vector<vm::vec3> &normals, std::vector<uint32_t> &indices) {
    // const index = this.index;
    // const positionAttribute = this.getAttribute( 'position' );

    // reset existing normals to zero
    // for ( let i = 0, il = normalAttribute.count; i < il; i ++ ) {
    //   normalAttribute.setXYZ( i, 0, 0, 0 );
    // }
    // std::fill(normals.begin(), normals.end(), 0);

    // const pA = new Vector3(), pB = new Vector3(), pC = new Vector3();
    // const nA = new Vector3(), nB = new Vector3(), nC = new Vector3();
    // const cb = new Vector3(), ab = new Vector3();

    // indexed elements
    for (size_t i = 0, il = indices.size(); i < il; i += 3) {
        // const uint32_t vA = index.getX(i + 0);
        // const uint32_t vB = index.getX(i + 1);
        // const uint32_t vC = index.getX(i + 2);

        const uint32_t &vA = indices[i];
        const uint32_t &vB = indices[i + 1];
        const uint32_t &vC = indices[i + 2];

        // pA.fromBufferAttribute( positionAttribute, vA );
        // pB.fromBufferAttribute( positionAttribute, vB );
        // pC.fromBufferAttribute( positionAttribute, vC );
        Vec pA{
            positions[vA].x,
            positions[vA].y,
            positions[vA].z
        };
        Vec pB{
            positions[vB].x,
            positions[vB].y,
            positions[vB].z
        };
        Vec pC{
            positions[vC].x,
            positions[vC].y,
            positions[vC].z
        };

        // cb.subVectors( pC, pB );
        // ab.subVectors( pA, pB );
        // cb.cross( ab );
        Vec cb = pC - pB;
        Vec ab = pA - pB;
        // cb.cross(ab);
        cb ^= ab;

        // nA.fromBufferAttribute( normalAttribute, vA );
        // nB.fromBufferAttribute( normalAttribute, vB );
        // nC.fromBufferAttribute( normalAttribute, vC );
        Vec nA{
            normals[vA].x,
            normals[vA].y,
            normals[vA].z
        };
        Vec nB{
            normals[vB].x,
            normals[vB].y,
            normals[vB].z
        };
        Vec nC{
            normals[vC].x,
            normals[vC].y,
            normals[vC].z
        };

        // nA.add( cb );
        // nB.add( cb );
        // nC.add( cb );
        nA += cb;
        nB += cb;
        nC += cb;

        // normalAttribute.setXYZ( vA, nA.x, nA.y, nA.z );
        // normalAttribute.setXYZ( vB, nB.x, nB.y, nB.z );
        // normalAttribute.setXYZ( vC, nC.x, nC.y, nC.z );
        normals[vA].x = nA.x;
        normals[vA].y = nA.y;
        normals[vA].z = nA.z;
        normals[vB].x = nB.x;
        normals[vB].y = nB.y;
        normals[vB].z = nB.z;
        normals[vC].x = nC.x;
        normals[vC].y = nC.y;
        normals[vC].z = nC.z;
    }

    normalizeNormals(normals);
}
/* void fillVec3(std::vector<vm::vec3> &array, const vm::vec3 &v) {
    for (size_t i = 0, il = array.size(); i < il; i++) {
        array[i] = v;
    }
} */

//

enum class WindingDirection {
    CCW,
    CW
};
enum class ComputeNormals {
    NO,
    YES
};

template<
    typename T,
    typename G,
    WindingDirection windingDirection,
    ComputeNormals computeNormals
>
void createPlaneGeometry(
    int width,
    int height,
    int widthSegments,
    int heightSegments,
    int rowSize,
    const std::vector<T> &heightfields,
    G &geometry
) {
    const int &gridX = widthSegments; // equals chunkSize - 1
    const int &gridY = heightSegments; // equals chunkSize - 1

    const int gridX1 = gridX + 1; // equals chunkSize
    const int gridY1 = gridY + 1; // equals chunkSize

    const int segment_width = width / gridX; // equals lod
    const int segment_height = height / gridY; // equals lod

    //

    auto pushPoint = [&](int x, int y) -> void {
        const int dx = x + 1;
        const int dy = y + 1;
        const int index = dx + dy * rowSize;
        const T &v0 = heightfields[index];

        const int ax = x * segment_width;
        const int ay = y * segment_height;

        // position
        const float height = v0.getHeight();

        const MaterialsArray &materials = v0.materials;
        const MaterialsWeightsArray &materialWeights = v0.materialsWeights;

        geometry.positions.push_back(vm::vec3{
            (float)ax,
            height,
            (float)ay
        });

        geometry.materials.push_back(vm::ivec4{
            materials[0],
            materials[1],
            materials[2],
            materials[3]
        });

        geometry.materialsWeights.push_back(vm::vec4{
            materialWeights[0],
            materialWeights[1],
            materialWeights[2],
            materialWeights[3]
        });

        // normal
        vm::vec3 normal;
        if (computeNormals == ComputeNormals::YES) {
            const int Lx = (x - 1) + 1;
            const int Ly = y + 1;
            const int Lindex = Lx + Ly * rowSize;
            const float Lheight = heightfields[Lindex].getHeight();

            const int Rx = (x + 1) + 1;
            const int Ry = y + 1;
            const int Rindex = Rx + Ry * rowSize;
            const float Rheight = heightfields[Rindex].getHeight();

            const int Ux = x + 1;
            const int Uy = (y - 1) + 1;
            const int Uindex = Ux + Uy * rowSize;
            const float Uheight = heightfields[Uindex].getHeight();

            const int Dx = x + 1;
            const int Dy = (y + 1) + 1;
            const int Dindex = Dx + Dy * rowSize;
            const float Dheight = heightfields[Dindex].getHeight();

            normal = vm::normalize(
                vm::vec3{
                    Lheight - Rheight,
                    2.0f,
                    Uheight - Dheight
                }
            );
        } else {
            normal = vm::vec3{
                0,
                1,
                0
            };
        }
        geometry.normals.push_back(normal);
        
        // metadata
        geometry.pushPointMetadata(v0);
    };

    // positions
    for (int y = 0; y < gridY1; y++) {
        for (int x = 0; x < gridX1; x++) {
            pushPoint(x, y);
        }
    }

    auto pushTriangle = [&](int ra, int rb, int rc, int wa, int wb, int wc) -> void {
        const T &hfA = heightfields[ra];
        const T &hfB = heightfields[rb];
        const T &hfC = heightfields[rc];
        if (T::acceptIndices(hfA, hfB, hfC)) {
            if (windingDirection == WindingDirection::CCW) {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wb);
                geometry.indices.push_back(wc);
            } else {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wc);
                geometry.indices.push_back(wb);
            }
        }
    };

    // indices
    for (int y = 0; y < gridY; y++) {
        for (int x = 0; x < gridX; x++) {
            const int dx = x + 1;
            const int dy = y + 1;
            
            const int ra = dx + rowSize * dy;
            const int rb = dx + rowSize * (dy + 1);
            const int rc = (dx + 1) + rowSize * (dy + 1);
            const int rd = (dx + 1) + rowSize * dy;

            const int wa = x + gridX1 * y;
            const int wb = x + gridX1 * (y + 1);
            const int wc = (x + 1) + gridX1 * (y + 1);
            const int wd = (x + 1) + gridX1 * y;

            pushTriangle(ra, rb, rd, wa, wb, wd);
            pushTriangle(rb, rc, rd, wb, wc, wd);

        }

    }
}
template<
    typename T,
    typename G,
    WindingDirection windingDirection,
    ComputeNormals computeNormals
>
void createPlaneSeamsGeometry(
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    int rowSize,
    const std::vector<T> &heightfields,
    G &geometry
) {
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = chunkSize * lod / bottomLod;
    const int gridWidthP1 = gridWidth + 1;
    const int gridWidthP3 = gridWidth + 3;

    const int gridHeight = chunkSize * lod / rightLod;
    const int gridHeightP1 = gridHeight + 1;
    const int gridHeightP3 = gridHeight + 3;

    const int heightfieldsCenterDataWriteOffset = chunkSize * chunkSize;
    const int heightfieldsCenterDataReadOffset = rowSize * rowSize;

    //

    auto pushBottomPoint = [&](int x, int y) -> void {
        const int ax = x * bottomLod;
        const int ay = (y + chunkSize) * lod;

        const int dx = x + 1;
        const int dy = y + 1;
        const int index =
            heightfieldsCenterDataReadOffset +
            (dy * gridWidthP3) +
            dx;
        const T &v0 = heightfields[index];

        // position
        const float height = v0.getHeight();

        const MaterialsArray &materials = v0.materials;
        const MaterialsWeightsArray &materialWeights = v0.materialsWeights;

        geometry.positions.push_back(vm::vec3{
            (float)ax,
            height,
            (float)ay
        });

        geometry.materials.push_back(vm::ivec4{
            materials[0],
            materials[1],
            materials[2],
            materials[3]
        });

        geometry.materialsWeights.push_back(vm::vec4{
            materialWeights[0],
            materialWeights[1],
            materialWeights[2],
            materialWeights[3]
        });

        // normal
        vm::vec3 normal;
        if (computeNormals == ComputeNormals::YES) {
            const int Lx = dx - 1;
            const int Ly = dy;
            const int Lindex =
                heightfieldsCenterDataReadOffset +
                (Ly * gridWidthP3) +
                Lx;
            const float Lheight = heightfields[Lindex].getHeight();

            const int Rx = dx + 1;
            const int Ry = dy;
            const int Rindex =
                heightfieldsCenterDataReadOffset +
                (Ry * gridWidthP3) +
                Rx;
            const float Rheight = heightfields[Rindex].getHeight();

            const int Ux = dx;
            const int Uy = dy - 1;
            const int Uindex =
                heightfieldsCenterDataReadOffset +
                (Uy * gridWidthP3) +
                Ux;
            const float Uheight = heightfields[Uindex].getHeight();

            const int Dx = dx;
            const int Dy = dy + 1;
            const int Dindex =
                heightfieldsCenterDataReadOffset +
                (Dy * gridWidthP3) +
                Dx;
            const float Dheight = heightfields[Dindex].getHeight();

            normal = vm::normalize(
                vm::vec3{
                    Lheight - Rheight,
                    2.0f,
                    Uheight - Dheight
                }
            );
        } else {
            normal = vm::vec3{
                0,
                1,
                0
            };
        }
        geometry.normals.push_back(normal);

        // metadata
        geometry.pushPointMetadata(v0);
    };
    auto pushRightPoint = [&](int x, int y) -> void {
        const int ax = (x + chunkSize) * lod;
        const int ay = y * rightLod;

        const int dx = x + 1;
        const int dy = y + 1;
        const int index =
            heightfieldsCenterDataReadOffset +
            (3 * gridWidthP3) +
            (dx * gridHeightP3) +
            dy;
        const T &v0 = heightfields[index];

        // position
        const float height = v0.getHeight();

        const MaterialsArray &materials = v0.materials;
        const MaterialsWeightsArray &materialWeights = v0.materialsWeights;

        geometry.positions.push_back(vm::vec3{
            (float)ax,
            height,
            (float)ay
        });

        geometry.materials.push_back(vm::ivec4{
            materials[0],
            materials[1],
            materials[2],
            materials[3]
        });

        geometry.materialsWeights.push_back(vm::vec4{
            materialWeights[0],
            materialWeights[1],
            materialWeights[2],
            materialWeights[3]
        });
        
        // normal
        vm::vec3 normal;
        if (computeNormals == ComputeNormals::YES) {
            const int Lx = dx - 1;
            const int Ly = dy;
            const int Lindex =
                heightfieldsCenterDataReadOffset +
                (3 * gridWidthP3) +
                (Lx * gridHeightP3) +
                Ly;
            const float Lheight = heightfields[Lindex].getHeight();

            const int Rx = dx + 1;
            const int Ry = dy;
            const int Rindex =
                heightfieldsCenterDataReadOffset +
                (3 * gridWidthP3) +
                (Rx * gridHeightP3) +
                Ry;
            const float Rheight = heightfields[Rindex].getHeight();

            const int Ux = dx;
            const int Uy = dy - 1;
            const int Uindex =
                heightfieldsCenterDataReadOffset +
                (3 * gridWidthP3) +
                (Ux * gridHeightP3) +
                Uy;
            const float Uheight = heightfields[Uindex].getHeight();

            const int Dx = dx;
            const int Dy = dy + 1;
            const int Dindex =
                heightfieldsCenterDataReadOffset +
                (3 * gridWidthP3) +
                (Dx * gridHeightP3) +
                Dy;
            const float Dheight = heightfields[Dindex].getHeight();

            normal = vm::normalize(
                vm::vec3{
                    Lheight - Rheight,
                    2.0f,
                    Uheight - Dheight
                }
            );
        } else {
            normal = vm::vec3{
                0,
                1,
                0
            };
        }
        geometry.normals.push_back(normal);

        // metadata
        geometry.pushPointMetadata(v0);
    };

    // positions
    const int heightfieldsBottomDataWriteOffset =
        heightfieldsCenterDataWriteOffset +
        gridWidthP1;
    const int heightfieldsBottomDataReadOffset =
        heightfieldsCenterDataReadOffset +
        3 * gridWidthP3;
    // bottom
    {
        const int y = 0;
        for (int x = 0; x < gridWidthP1; x++) {
            pushBottomPoint(x, y);
        }
    }
    // right
    {
        const int x = 0;
        for (int y = 0; y < gridHeightP1; y++) {
            pushRightPoint(x, y);
        }
    }

    //

    const int chunkSizeM1 = chunkSize - 1;
    const int gridX = chunkSizeM1;
    const int gridY = chunkSizeM1;

    const int gridX1 = gridX + 1; // equals chunkSize
    const int gridY1 = gridY + 1; // equals chunkSize

    auto pushTriangle = [&](int ra, int rb, int rc, int wa, int wb, int wc) {
        const T &hfA = heightfields[ra];
        const T &hfB = heightfields[rb];
        const T &hfC = heightfields[rc];
        if (T::acceptIndices(hfA, hfB, hfC)) {
            if (windingDirection == WindingDirection::CCW) {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wb);
                geometry.indices.push_back(wc);
            } else {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wc);
                geometry.indices.push_back(wb);
            }
        }
    };

    // indices
    // bottom
    if (bottomLod == lod) {
        const int innerPointY = chunkSize - 1;
        for (int innerPointX = 0; innerPointX < chunkSize; innerPointX++) {
            const int outerPointX = innerPointX;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            const int outerPointDx = outerPointX + 1;

            // inner
            const int ra = innerPointDx + rowSize * innerPointDy;
            const int rd = (innerPointDx + 1) + rowSize * innerPointDy;
            // outer
            const int rb = heightfieldsCenterDataReadOffset + outerPointDx;
            const int rc = heightfieldsCenterDataReadOffset + (outerPointDx + 1);

            const int wa = innerPointX + chunkSize * innerPointY;
            const int wd = (innerPointX + 1) + chunkSize * innerPointY;
            const int wb = heightfieldsCenterDataWriteOffset + outerPointX;
            const int wc = heightfieldsCenterDataWriteOffset + (outerPointX + 1);

            pushTriangle(ra, rb, rc, wa, wb, wc);
            if (innerPointX != (chunkSize - 1)) { // only single triangle in corner
                pushTriangle(ra, rc, rd, wa, wc, wd);
            }
        }
    } else if (bottomLod > lod) {
        const int innerPointY = chunkSize - 1;
        for (int innerPointX = 0; innerPointX < chunkSize; innerPointX++) {
            const int outerPointX = innerPointX / 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            const int outerPointDx = outerPointX + 1;

            if (innerPointX % 2 == 0) {
                // inner
                const int ra = innerPointDx + rowSize * innerPointDy;
                const int rd = (innerPointDx + 1) + rowSize * innerPointDy;
                // outer
                const int rb = heightfieldsCenterDataReadOffset + outerPointDx;
                const int rc = heightfieldsCenterDataReadOffset + (outerPointDx + 1);

                const int wa = innerPointX + chunkSize * innerPointY;
                const int wd = (innerPointX + 1) + chunkSize * innerPointY;
                const int wb = heightfieldsCenterDataWriteOffset + outerPointX;
                const int wc = heightfieldsCenterDataWriteOffset + (outerPointX + 1);

                pushTriangle(ra, rb, rc, wa, wb, wc);
                pushTriangle(ra, rc, rd, wa, wc, wd);
            } else {
                if (innerPointX != (chunkSize - 1)) {
                    // inner
                    const int ra = innerPointDx + rowSize * innerPointDy;
                    const int rd = (innerPointDx + 1) + rowSize * innerPointDy;
                    // outer
                    const int rb = heightfieldsCenterDataReadOffset + outerPointDx;
                    const int rc = heightfieldsCenterDataReadOffset + (outerPointDx + 1);

                    const int wa = innerPointX + chunkSize * innerPointY;
                    const int wd = (innerPointX + 1) + chunkSize * innerPointY;
                    const int wb = heightfieldsCenterDataWriteOffset + outerPointX;
                    const int wc = heightfieldsCenterDataWriteOffset + (outerPointX + 1);

                    pushTriangle(ra, rc, rd, wa, wc, wd);
                }
            }
        }
    } else if (bottomLod < lod) {
        const int innerPointY = chunkSize - 1;
        for (int innerPointX = 0; innerPointX < chunkSize; innerPointX++) {
            int outerPointX = innerPointX * 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            int outerPointDx = outerPointX + 1;

            if (innerPointX != (chunkSize - 1)) {
                {
                    // inner
                    const int ra = innerPointDx + rowSize * innerPointDy;
                    const int rd = (innerPointDx + 1) + rowSize * innerPointDy;
                    // outer
                    const int rb = heightfieldsCenterDataReadOffset + outerPointDx;
                    const int rc = heightfieldsCenterDataReadOffset + (outerPointDx + 1);

                    const int wa = innerPointX + chunkSize * innerPointY;
                    const int wd = (innerPointX + 1) + chunkSize * innerPointY;
                    const int wb = heightfieldsCenterDataWriteOffset + outerPointX;
                    const int wc = heightfieldsCenterDataWriteOffset + (outerPointX + 1);

                    pushTriangle(ra, rb, rc, wa, wb, wc);
                    pushTriangle(ra, rc, rd, wa, wc, wd);
                }
                outerPointX++;
                outerPointDx++;
                {
                    // inner
                    const int ra = innerPointDx + rowSize * innerPointDy;
                    const int rd = (innerPointDx + 1) + rowSize * innerPointDy;
                    // outer
                    const int rb = heightfieldsCenterDataReadOffset + outerPointDx;
                    const int rc = heightfieldsCenterDataReadOffset + (outerPointDx + 1);

                    const int wa = innerPointX + chunkSize * innerPointY;
                    const int wd = (innerPointX + 1) + chunkSize * innerPointY;
                    const int wb = heightfieldsCenterDataWriteOffset + outerPointX;
                    const int wc = heightfieldsCenterDataWriteOffset + (outerPointX + 1);

                    pushTriangle(rb, rc, rd, wb, wc, wd);
                }
            } else {
                // inner
                const int ra = innerPointDx + rowSize * innerPointDy;
                // outer
                const int rb = heightfieldsCenterDataReadOffset + outerPointDx;
                const int rc = heightfieldsCenterDataReadOffset + (outerPointDx + 1);
                const int rd = heightfieldsCenterDataReadOffset + (outerPointDx + 2);

                const int wa = innerPointX + chunkSize * innerPointY;
                const int wb = heightfieldsCenterDataWriteOffset + outerPointX;
                const int wc = heightfieldsCenterDataWriteOffset + (outerPointX + 1);
                const int wd = heightfieldsCenterDataWriteOffset + (outerPointX + 2);

                pushTriangle(ra, rb, rc, wa, wb, wc);
                pushTriangle(ra, rc, rd, wa, wc, wd);
            }
        }
    }
    // right
    if (rightLod == lod) {
        const int innerPointX = chunkSize - 1;
        const int outerPointX = innerPointX;
        for (int innerPointY = 0; innerPointY < chunkSize; innerPointY++) {
            const int outerPointY = innerPointY;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            const int outerPointDy = outerPointY + 1;

            // inner
            const int ra = innerPointDx + rowSize * innerPointDy;
            const int rb = innerPointDx + rowSize * (innerPointDy + 1);
            // outer
            const int rd = heightfieldsBottomDataReadOffset + outerPointDy;
            const int rc = heightfieldsBottomDataReadOffset + (outerPointDy + 1);

            const int wa = innerPointX + chunkSize * innerPointY;
            const int wb = innerPointX + chunkSize * (innerPointY + 1);
            const int wd = heightfieldsBottomDataWriteOffset + outerPointY;
            const int wc = heightfieldsBottomDataWriteOffset + (outerPointY + 1);

            pushTriangle(ra, rc, rd, wa, wc, wd);
            if (innerPointY != (chunkSize - 1)) { // only single triangle in corner
                pushTriangle(ra, rb, rc, wa, wb, wc);
            }
        }
    } else if (rightLod > lod) {
        const int innerPointX = chunkSize - 1;
        const int outerPointX = innerPointX;
        for (int innerPointY = 0; innerPointY < chunkSize; innerPointY++) {
            const int outerPointY = innerPointY / 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            const int outerPointDy = outerPointY + 1;

            if (innerPointY % 2 == 0) {
                // inner
                const int ra = innerPointDx + rowSize * innerPointDy;
                const int rb = innerPointDx + rowSize * (innerPointDy + 1);
                // outer
                const int rd = heightfieldsBottomDataReadOffset + outerPointDy;
                const int rc = heightfieldsBottomDataReadOffset + (outerPointDy + 1);

                const int wa = innerPointX + chunkSize * innerPointY;
                const int wb = innerPointX + chunkSize * (innerPointY + 1);
                const int wd = heightfieldsBottomDataWriteOffset + outerPointY;
                const int wc = heightfieldsBottomDataWriteOffset + (outerPointY + 1);

                pushTriangle(ra, rc, rd, wa, wc, wd);
                pushTriangle(ra, rb, rc, wa, wb, wc);
            } else {
                if (innerPointY != (chunkSize - 1)) {
                    // inner
                    const int ra = innerPointDx + rowSize * innerPointDy;
                    const int rb = innerPointDx + rowSize * (innerPointDy + 1);
                    // outer
                    const int rd = heightfieldsBottomDataReadOffset + outerPointDy;
                    const int rc = heightfieldsBottomDataReadOffset + (outerPointDy + 1);

                    const int wa = innerPointX + chunkSize * innerPointY;
                    const int wb = innerPointX + chunkSize * (innerPointY + 1);
                    const int wd = heightfieldsBottomDataWriteOffset + outerPointY;
                    const int wc = heightfieldsBottomDataWriteOffset + (outerPointY + 1);

                    pushTriangle(ra, rb, rc, wa, wb, wc);
                }
            }
        }
    } else if (rightLod < lod) {
        const int innerPointX = chunkSize - 1;
        const int outerPointX = innerPointX;
        for (int innerPointY = 0; innerPointY < chunkSize; innerPointY++) {
            int outerPointY = innerPointY * 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            int outerPointDy = outerPointY + 1;

            if (innerPointY != (chunkSize - 1)) {
                {
                    // inner
                    const int ra = innerPointDx + rowSize * innerPointDy;
                    const int rb = innerPointDx + rowSize * (innerPointDy + 1);
                    // outer
                    const int rd = heightfieldsBottomDataReadOffset + outerPointDy;
                    const int rc = heightfieldsBottomDataReadOffset + (outerPointDy + 1);

                    const int wa = innerPointX + chunkSize * innerPointY;
                    const int wb = innerPointX + chunkSize * (innerPointY + 1);
                    const int wd = heightfieldsBottomDataWriteOffset + outerPointY;
                    const int wc = heightfieldsBottomDataWriteOffset + (outerPointY + 1);

                    pushTriangle(ra, rc, rd, wa, wc, wd);
                    pushTriangle(ra, rb, rc, wa, wb, wc);
                }
                outerPointY++;
                outerPointDy++;
                {
                    // inner
                    const int ra = innerPointDx + rowSize * innerPointDy;
                    const int rb = innerPointDx + rowSize * (innerPointDy + 1);
                    // outer
                    const int rd = heightfieldsBottomDataReadOffset + outerPointDy;
                    const int rc = heightfieldsBottomDataReadOffset + (outerPointDy + 1);

                    const int wa = innerPointX + chunkSize * innerPointY;
                    const int wb = innerPointX + chunkSize * (innerPointY + 1);
                    const int wd = heightfieldsBottomDataWriteOffset + outerPointY;
                    const int wc = heightfieldsBottomDataWriteOffset + (outerPointY + 1);

                    pushTriangle(rb, rc, rd, wb, wc, wd);
                }
            } else {
                // inner
                const int ra = innerPointDx + rowSize * innerPointDy;
                // outer
                const int rd = heightfieldsBottomDataReadOffset + outerPointDy;
                const int rc = heightfieldsBottomDataReadOffset + (outerPointDy + 1);
                const int rb = heightfieldsBottomDataReadOffset + (outerPointDy + 2);

                const int wa = innerPointX + chunkSize * innerPointY;
                const int wd = heightfieldsBottomDataWriteOffset + outerPointY;
                const int wc = heightfieldsBottomDataWriteOffset + (outerPointY + 1);
                const int wb = heightfieldsBottomDataWriteOffset + (outerPointY + 2);

                pushTriangle(ra, rb, rc, wa, wb, wc);
                pushTriangle(ra, rc, rd, wa, wc, wd);
            }
        }
    }
}

//

template<typename G>
void offsetGeometry(G &geometry, const vm::ivec2 &worldPosition, float height = -(float)WORLD_BASE_HEIGHT) {
    for (size_t i = 0; i < geometry.positions.size(); i++) {
        vm::vec3 &p = geometry.positions[i]; 
        p.x += (float)worldPosition.x;
        p.y += height;
        p.z += (float)worldPosition.y;
    }
}

//

void generateHeightfieldCenterMesh(
    int lod,
    int chunkSize,
    const std::vector<Heightfield> &heightfields,
    TerrainGeometry &geometry
) {
    const int worldSize = chunkSize * lod;
    const int worldSizeM1 = worldSize - lod;
    const int chunkSizeM1 = chunkSize - 1;
    const int rowSize = chunkSize + 2;
    createPlaneGeometry<
        Heightfield,
        TerrainGeometry,
        WindingDirection::CCW,
        ComputeNormals::YES
    >(
        worldSizeM1,
        worldSizeM1,
        chunkSizeM1,
        chunkSizeM1,
        rowSize,
        heightfields,
        geometry
    );
}
void generateHeightfieldSeamsMesh(
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    const std::vector<Heightfield> &heightfields,
    TerrainGeometry &geometry
) {
    const int rowSize = chunkSize + 2;
    createPlaneSeamsGeometry<
        Heightfield,
        TerrainGeometry,
        WindingDirection::CCW,
        ComputeNormals::YES
    >(
        lod,
        lodArray,
        chunkSize,
        rowSize,
        heightfields,
        geometry
    );
}

//

void generateWaterfieldCenterMesh(
    int lod,
    int chunkSize,
    const std::vector<Waterfield> &waterfields,
    WaterGeometry &geometry
) {
    const int worldSize = chunkSize * lod;
    const int worldSizeM1 = worldSize - lod;
    const int chunkSizeM1 = chunkSize - 1;
    const int rowSize = chunkSize + 2;
    createPlaneGeometry<
        Waterfield,
        WaterGeometry,
        WindingDirection::CCW,
        ComputeNormals::NO
    >(
        worldSizeM1,
        worldSizeM1,
        chunkSizeM1,
        chunkSizeM1,
        rowSize,
        waterfields,
        geometry
    );
}
void generateWaterfieldSeamsMesh(
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    const std::vector<Waterfield> &waterfields,
    WaterGeometry &geometry
) {
    const int rowSize = chunkSize + 2;
    createPlaneSeamsGeometry<
        Waterfield,
        WaterGeometry,
        WindingDirection::CCW,
        ComputeNormals::NO
    >(
        lod,
        lodArray,
        chunkSize,
        rowSize,
        waterfields,
        geometry
    );
}

//

template<typename G>
void createBoxGeometry(float width, float height, float depth, int widthSegments, int heightSegments, int depthSegments, G &geometry) {
    // segments

    // widthSegments = Math.floor(widthSegments);
    // heightSegments = Math.floor(heightSegments);
    // depthSegments = Math.floor(depthSegments);

    // buffers

    // const indices = [];
    // const vertices = [];
    // const normals = [];
    // const uvs = [];

    // helper variables

    int numberOfVertices = 0;
    // let groupStart = 0;

    auto buildPlane = [&](std::function<vm::vec3(float z, float y, float x)> &vectorFn, float udir, float vdir, float width, float height, float depth, int gridX, int gridY) -> void {
        const float segmentWidth = width / gridX;
        const float segmentHeight = height / gridY;

        const float widthHalf = width / 2.f;
        const float heightHalf = height / 2.f;
        const float depthHalf = depth / 2.f;

        const int gridX1 = gridX + 1;
        const int gridY1 = gridY + 1;

        int vertexCounter = 0;
        // let groupCount = 0;

        // const vector = new Vector3();
        vm::vec3 vector;

        // generate vertices, normals and uvs

        for (int iy = 0; iy < gridY1; iy++) {

            const int y = iy * segmentHeight - heightHalf;

            for (int ix = 0; ix < gridX1; ix++) {

                const int x = ix * segmentWidth - widthHalf;

                // set values to correct vector component

                vector = vectorFn(x * udir, y * vdir, depthHalf);
                // vector[ u ] = x * udir;
                // vector[ v ] = y * vdir;
                // vector[ w ] = depthHalf;

                // now apply vector to vertex buffer

                // vertices.push(vector.x, vector.y, vector.z);
                geometry.positions.push_back(vector);

                // set values to correct vector component

                vector = vectorFn(0, 0, depth > 0 ? 1 : - 1);
                // vector[ u ] = 0;
                // vector[ v ] = 0;
                // vector[ w ] = depth > 0 ? 1 : - 1;

                // now apply vector to normal buffer

                // normals.push(vector.x, vector.y, vector.z);
                geometry.normals.push_back(vector);

                // uvs

                // uvs.push((float)ix / (float)gridX);
                // uvs.push(1.f - ((float)iy / (float)gridY));
                /* geometry.uvs.push_back(vm::vec2{
                    (float)ix / (float)gridX,
                    1.f - ((float)iy / (float)gridY)
                }); */
                geometry.uvs.push_back(vm::vec2{
                    (float)ix * width,
                    (1.f - (float)iy) * height
                });

                // counters

                vertexCounter += 1;

            }

        }

        // indices

        // 1. you need three indices to draw a single face
        // 2. a single segment consists of two faces
        // 3. so we need to generate six (2*3) indices per segment

        for (int iy = 0; iy < gridY; iy++) {

            for (int ix = 0; ix < gridX; ix++) {

                const int a = numberOfVertices + ix + gridX1 * iy;
                const int b = numberOfVertices + ix + gridX1 * (iy + 1);
                const int c = numberOfVertices + (ix + 1) + gridX1 * (iy + 1);
                const int d = numberOfVertices + (ix + 1) + gridX1 * iy;

                // faces

                geometry.indices.push_back(a);
                geometry.indices.push_back(b);
                geometry.indices.push_back(d);
                geometry.indices.push_back(b);
                geometry.indices.push_back(c);
                geometry.indices.push_back(d);

                // increase counter

                // groupCount += 6;

            }

        }

        // add a group to the geometry. this will ensure multi material support

        // scope.addGroup( groupStart, groupCount, materialIndex );

        // calculate new start value for groups

        // groupStart += groupCount;

        // update total number of vertices

        numberOfVertices += vertexCounter;
    };

    // build each side of the box geometry

    std::function<vm::vec3(float z, float y, float x)> pushVertices1 = [&](float z, float y, float x) -> vm::vec3 {
        return vm::vec3{
            x,
            y,
            z
        };
    };
    buildPlane(pushVertices1, -1, -1, depth, height, width, depthSegments, heightSegments); // px
    std::function<vm::vec3(float z, float y, float x)> pushVertices2 = [&](float z, float y, float x) -> vm::vec3 {
        return vm::vec3{
            x,
            y,
            z
        };
    };
    buildPlane(pushVertices2, 1, -1, depth, height, -width, depthSegments, heightSegments); // nx
    // std::function<vm::vec3(float z, float y, float x)> pushVertices3 = [&](float x, float z, float y) -> vm::vec3 {
    //     return vm::vec3{
    //         x,
    //         y,
    //         z
    //     };
    // };
    // buildPlane(pushVertices3, 1, 1, width, depth, height, widthSegments, depthSegments); // py
    // std::function<vm::vec3(float z, float y, float x)> pushVertices4 = [&](float x, float z, float y) -> vm::vec3 {
    //     return vm::vec3{
    //         x,
    //         y,
    //         z
    //     };
    // };
    // buildPlane(pushVertices4, 1, -1, width, depth, -height, widthSegments, depthSegments); // ny
    std::function<vm::vec3(float z, float y, float x)> pushVertices5 = [&](float x, float y, float z) -> vm::vec3 {
        return vm::vec3{
            x,
            y,
            z
        };
    };
    buildPlane(pushVertices5, 1, -1, width, height, depth, widthSegments, heightSegments); // pz
    std::function<vm::vec3(float z, float y, float x)> pushVertices6 = [&](float x, float y, float z) -> vm::vec3 {
        return vm::vec3{
            x,
            y,
            z
        };
    };
    buildPlane(pushVertices6, -1, -1, width, height, -depth, widthSegments, heightSegments); // nz

    // build geometry

    // this.setIndex( indices );
    // this.setAttribute( 'position', new Float32BufferAttribute( vertices, 3 ) );
    // this.setAttribute( 'normal', new Float32BufferAttribute( normals, 3 ) );
    // this.setAttribute( 'uv', new Float32BufferAttribute( uvs, 2 ) );
}

//

/* template<typename G>
G &mergeGeometry(G &dst, G &a, G &b) {
    dst.positions.reserve(a.positions.size() + b.positions.size());
    dst.normals.reserve(a.normals.size() + b.normals.size());
    dst.indices.reserve(a.indices.size() + b.indices.size());
    
    for (size_t i = 0; i < a.positions.size(); i++) {
        dst.positions.push_back(a.positions[i]);
    }
    for (size_t i = 0; i < a.normals.size(); i++) {
        dst.normals.push_back(a.normals[i]);
    }
    for (size_t i = 0; i < a.indices.size(); i++) {
        dst.indices.push_back(a.indices[i]);
    }

    uint32_t indexOffset = (uint32_t)a.positions.size();
    for (size_t i = 0; i < b.positions.size(); i++) {
        dst.positions.push_back(b.positions[i]);
    }
    for (size_t i = 0; i < b.normals.size(); i++) {
        dst.normals.push_back(b.normals[i]);
    }
    for (size_t i = 0; i < b.indices.size(); i++) {
        dst.indices.push_back(indexOffset + b.indices[i]);
    }

    return dst;
} */
template<typename G>
G &mergeGeometries(G &dst, std::vector<G> &geometries) {
    size_t numPositions = 0;
    for (size_t i = 0; i < geometries.size(); i++) {
        numPositions += geometries[i].positions.size();
    }

    size_t numNormals = 0;
    for (size_t i = 0; i < geometries.size(); i++) {
        numNormals += geometries[i].normals.size();
    }
    
    size_t numIndices = 0;
    for (size_t i = 0; i < geometries.size(); i++) {
        numIndices += geometries[i].indices.size();
    }

    size_t numUvs = 0;
    for (size_t i = 0; i < geometries.size(); i++) {
        numUvs += geometries[i].uvs.size();
    }

    size_t numPositions2D = 0;
    for (size_t i = 0; i < geometries.size(); i++) {
        numPositions2D += geometries[i].positions2D.size();
    }

    dst.positions.reserve(numPositions);
    dst.normals.reserve(numNormals);
    dst.indices.reserve(numIndices);
    
    for (size_t i = 0; i < geometries.size(); i++) {
        G &g = geometries[i];
        for (size_t j = 0; j < g.positions.size(); j++) {
            dst.positions.push_back(g.positions[j]);
        }
        for (size_t j = 0; j < g.normals.size(); j++) {
            dst.normals.push_back(g.normals[j]);
        }
        size_t positionOffset = dst.positions.size() / 3;
        for (size_t j = 0; j < g.indices.size(); j++) {
            dst.indices.push_back(positionOffset + g.indices[j]);
        }
        for (size_t j = 0; j < g.uvs.size(); j++) {
            dst.uvs.push_back(g.uvs[j]);
        }
        for (size_t j = 0; j < g.positions2D.size(); j++) {
            dst.positions2D.push_back(g.positions2D[j]);
        }
    }

    return dst;
}
void setPositions2D(BarrierGeometry &geometry, const vm::ivec2 position2D) {
    geometry.positions2D.reserve(geometry.positions.size());
    for (size_t i = 0; i < geometry.positions.size(); i++) {
        geometry.positions2D.push_back(position2D);
    }
}

//

void generateBarrierGeometry(
    const vm::ivec2 &worldPosition,
    int chunkSize,
    OctreeContext &octreeContext,
    BarrierGeometry &geometry
) {
    std::vector<OctreeNodePtr> seedLeafNodes = octreeContext.getLeafNodes();

    // geometry
    std::vector<BarrierGeometry> geometries;
    for (size_t i = 0; i < seedLeafNodes.size(); i++) {
        OctreeNodePtr node = seedLeafNodes[i];

        BarrierGeometry g;

        /* {
            std::cout << "main leaf node: ";
            std::cout << node->min.x << "," << node->min.y << "/" << node->lod << ", ";
            std::cout << std::endl;
        }
        {
            std::cout << "seed leaf nodes: ";
            for (auto iter = seedLeafNodes.begin(); iter != seedLeafNodes.end(); iter++) {
                OctreeNodePtr node = *iter;
                std::cout << node->min.x << "," << node->min.y << "/" << node->lod << ", ";
            }
            std::cout << std::endl;
        } */

        const vm::ivec2 &nodePosition = node->min;
        const int &nodeLod = node->lod;

        /* float barrierMinHeight = std::numeric_limits<float>::infinity();
        float barrierMaxHeight = -std::numeric_limits<float>::infinity();
        for (int dz = 0; dz < chunkSize; dz++) {
            for (int dx = 0; dx < chunkSize; dx++) {
                const float height = inst->getHeight(
                    nodePosition.x + dx * nodeLod,
                    nodePosition.y + dz * nodeLod
                ) - (float)WORLD_BASE_HEIGHT;
                if (height < barrierMinHeight) {
                    barrierMinHeight = height;
                }
                if (height > barrierMaxHeight) {
                    barrierMaxHeight = height;
                }
            }
        }
        barrierMinHeight = std::floor(barrierMinHeight);
        barrierMaxHeight = std::ceil(barrierMaxHeight); */
        constexpr int barrierMinHeight = -WORLD_BASE_HEIGHT + MIN_WORLD_HEIGHT;
        constexpr int barrierMaxHeight = -WORLD_BASE_HEIGHT + MAX_WORLD_HEIGHT;

        int width = nodeLod;
        constexpr int height = barrierMaxHeight - barrierMinHeight;
        int depth = nodeLod;
        createBoxGeometry(
            width,
            height,
            depth,
            1,
            1,
            1,
            g
        );
        vm::ivec2 worldOffset{
            width / 2 + worldPosition.x,
            depth / 2 + worldPosition.y
        };
        offsetGeometry(
            g,
            worldOffset,
            height / 2.f + barrierMinHeight
        );
        setPositions2D(g, node->min / chunkSize);

        geometries.push_back(std::move(g));
    }
    mergeGeometries(geometry, geometries);

    // leaf nodes
    geometry.leafNodes = seedLeafNodes;

    // leaf node index
    {
        vm::ivec2 leafNodesMin{ // in chunks space
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max()
        };
        vm::ivec2 leafNodesMax{ // in chunks space
            std::numeric_limits<int>::min(),
            std::numeric_limits<int>::min()
        };
        for (size_t i = 0; i < seedLeafNodes.size(); i++) {
            OctreeNodePtr node = seedLeafNodes[i];
            vm::ivec2 nodeChunkPosition = node->min; // in chunks space
            const int &nodeLod = node->lod;

            int minX = nodeChunkPosition.x;
            int minZ = nodeChunkPosition.y;
            int maxX = nodeChunkPosition.x + nodeLod;
            int maxZ = nodeChunkPosition.y + nodeLod;

            leafNodesMin.x = std::min(leafNodesMin.x, minX);
            leafNodesMin.y = std::min(leafNodesMin.y, minZ);
            leafNodesMax.x = std::max(leafNodesMax.x, maxX);
            leafNodesMax.y = std::max(leafNodesMax.y, maxZ);
        }
        geometry.leafNodesMin = leafNodesMin; // in chunks space
        geometry.leafNodesMax = leafNodesMax; // in chunks space

        int w = leafNodesMax.x - leafNodesMin.x; // in chunks space
        int h = leafNodesMax.y - leafNodesMin.y; // in chunks space
        geometry.leafNodesIndex.resize(w * h);

        for (size_t i = 0; i < seedLeafNodes.size(); i++) {
            OctreeNodePtr node = seedLeafNodes[i];
            vm::ivec2 nodeChunkPosition = node->min; // in chunks space
            const int &nodeLod = node->lod;

            int minX = nodeChunkPosition.x;
            int minZ = nodeChunkPosition.y;
            int maxX = nodeChunkPosition.x + nodeLod;
            int maxZ = nodeChunkPosition.y + nodeLod;

            for (int z = minZ; z < maxZ; z++) {
                for (int x = minX; x < maxX; x++) {
                    int dx = x - leafNodesMin.x;
                    int dz = z - leafNodesMin.y;
                    int index = dx + dz * w;
                    geometry.leafNodesIndex[index] = i;
                }
            }
        }
    }
}

//

void generateTerrainGeometry(
    const vm::ivec2 &worldPosition,
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    const std::vector<Heightfield> &heightfields,
    TerrainGeometry &geometry
) {
    generateHeightfieldCenterMesh(lod, chunkSize, heightfields, geometry);
    generateHeightfieldSeamsMesh(lod, lodArray, chunkSize, heightfields, geometry);
    offsetGeometry(geometry, worldPosition);
    // computeVertexNormals(geometry.positions, geometry.normals, geometry.indices);
}

//

void generateWaterGeometry(
    const vm::ivec2 &worldPosition,
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    const std::vector<Waterfield> &waterfields,
    WaterGeometry &geometry
) {
    generateWaterfieldCenterMesh(lod, chunkSize, waterfields, geometry);
    generateWaterfieldSeamsMesh(lod, lodArray, chunkSize, waterfields, geometry);
    offsetGeometry(geometry, worldPosition);
}

//

class HeightfieldSampler {
public:
    const vm::vec2 &worldPositionXZ;
    const int &lod;
    const int chunkSize;
    const int chunkSizeP2;
    const std::vector<Heightfield> &heightfields;
    
    HeightfieldSampler(
        const vm::vec2 &worldPositionXZ,
        const int &lod,
        const int chunkSize,
        const std::vector<Heightfield> &heightfields
    ) :
        worldPositionXZ(worldPositionXZ),
        lod(lod),
        chunkSize(chunkSize),
        chunkSizeP2(chunkSize + 2),
        heightfields(heightfields)
        {}
    float getHeight(float x, float z) {
        vm::vec2 location{
            x - worldPositionXZ.x,
            z - worldPositionXZ.y
        };
        location /= (float)lod;

        float result = bilinear<HeightfieldSampler, float>(location, chunkSize, *this);
        return result;
    }
    float get(int x, int z) {
        const int dx = x + 1;
        const int dz = z + 1;
        const int index = dx + dz * chunkSizeP2;
        const Heightfield &heightfield = heightfields[index];
        const float &height = heightfield.heightField;
        return height;
    }
};

//

void generateVegetationInstances(
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int chunkSize,
    const int numVegetationInstances,
    const std::vector<Heightfield> &heightfields,
    Noises &noises,
    VegetationGeometry &vegetationGeometry
) {
    constexpr int maxNumVeggiesPerChunk = 8;
    constexpr float maxVeggieRate = 0.35;
    // const float veggieRate = maxVeggieRate / (float)(lod * lod);
    const float veggieRate = maxVeggieRate / (float)lod;
    // const float veggieRate = maxVeggieRate;

    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y
    };
    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields
    );
    
    for (int dz = 0; dz < lod; dz++) {
        for (int dx = 0; dx < lod; dx++) {
            int chunkMinX = baseMinX + dx * chunkSize;
            int chunkMinZ = baseMinZ + dz * chunkSize;

            float chunkSeed = noises.vegetationSeedNoise.in2D(chunkMinX, chunkMinZ);
            unsigned int seedInt = *(unsigned int *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < maxNumVeggiesPerChunk; i++) {
                float noiseValue = dis(rng);
                float chunkOffsetX = dis(rng) * (float)chunkSize;
                float chunkOffsetZ = dis(rng) * (float)chunkSize;
                float rot = dis(rng) * 2.0f * M_PI;
                int instanceId = (int)std::round(dis(rng) * (float)(numVegetationInstances - 1));

                if (noiseValue < veggieRate) {
                    auto iterPair = vegetationGeometry.instances.emplace(
                        std::make_pair(instanceId, SplatInstance{})
                    );
                    auto iter = iterPair.first;
                    const bool &inserted = iterPair.second;
                    SplatInstance &instance = iter->second;
                    if (inserted) {
                        instance.instanceId = instanceId;
                    }

                    float ax = (float)chunkMinX + chunkOffsetX;
                    float az = (float)chunkMinZ + chunkOffsetZ;
                    const float height = heightfieldSampler.getHeight(ax, az) - (float)WORLD_BASE_HEIGHT;

                    instance.ps.push_back(ax);
                    instance.ps.push_back(height);
                    instance.ps.push_back(az);

                    Quat q = Quat().setFromAxisAngle(Vec{0, 1, 0}, rot);
                    instance.qs.push_back(q.x);
                    instance.qs.push_back(q.y);
                    instance.qs.push_back(q.z);
                    instance.qs.push_back(q.w);
                }
            }
        }
    }
}

//

void generateGrassInstances(
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int chunkSize,
    const int numGrassInstances,
    const std::vector<Heightfield> &heightfields,
    Noises &noises,
    GrassGeometry &grassGeometry
) {
    constexpr int maxNumGrassesPerChunk = 2048;
    constexpr float grassRate = 0.5;
    // const float grassRate = maxGrassRate / (float)(lod * lod);
    const float grassThrowRate = 1.f / (float)lod;
    // const float grassRate = maxGrassRate;

    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y
    };
    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields
    );
    
    for (int dz = 0; dz < lod; dz++) {
        for (int dx = 0; dx < lod; dx++) {
            int chunkMinX = baseMinX + dx * chunkSize;
            int chunkMinZ = baseMinZ + dz * chunkSize;

            float chunkSeed = noises.grassSeedNoise.in2D(chunkMinX, chunkMinZ);
            unsigned int seedInt = *(unsigned int *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < maxNumGrassesPerChunk; i++) {
                float chunkOffsetX = dis(rng) * (float)chunkSize;
                float chunkOffsetZ = dis(rng) * (float)chunkSize;
                float rot = dis(rng) * 2.0f * M_PI;
                int instanceId = (int)std::round(dis(rng) * (float)(numGrassInstances - 1));
                float throwNoise = dis(rng);
                
                float ax = (float)chunkMinX + chunkOffsetX;
                float az = (float)chunkMinZ + chunkOffsetZ;
                float noiseValue = noises.grassNoise.in2D(ax, az);

                if (noiseValue < grassRate && throwNoise <= grassThrowRate) {
                    auto iterPair = grassGeometry.instances.emplace(
                        std::make_pair(instanceId, SplatInstance{})
                    );
                    auto iter = iterPair.first;
                    const bool &inserted = iterPair.second;
                    SplatInstance &instance = iter->second;
                    if (inserted) {
                        instance.instanceId = instanceId;
                    }
                    
                    const float height = heightfieldSampler.getHeight(ax, az) - (float)WORLD_BASE_HEIGHT;

                    instance.ps.push_back(ax);
                    instance.ps.push_back(height);
                    instance.ps.push_back(az);

                    Quat q = Quat().setFromAxisAngle(Vec{0, 1, 0}, rot);
                    instance.qs.push_back(q.x);
                    instance.qs.push_back(q.y);
                    instance.qs.push_back(q.z);
                    instance.qs.push_back(q.w);
                }
            }
        }
    }
}

//

void generatePoiInstances(
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int chunkSize,
    const int numPoiInstances,
    const std::vector<Heightfield> &heightfields,
    Noises &noises,
    PoiGeometry &poiGeometry
) {
    constexpr int maxNumPoisPerChunk = 16;
    constexpr float poiRate = 0.3;
    // const float poiRate = maxPoiRate / (float)(lod * lod);
    const float poiThrowRate = 1.f / (float)lod;
    // const float poiRate = maxPoiRate;

    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y
    };
    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields
    );
    
    for (int dz = 0; dz < lod; dz++) {
        for (int dx = 0; dx < lod; dx++) {
            int chunkMinX = baseMinX + dx * chunkSize;
            int chunkMinZ = baseMinZ + dz * chunkSize;

            float chunkSeed = noises.poiSeedNoise.in2D(chunkMinX, chunkMinZ);
            unsigned int seedInt = *(unsigned int *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < maxNumPoisPerChunk; i++) {
                float noiseValue = dis(rng);
                float chunkOffsetX = dis(rng) * (float)chunkSize;
                float chunkOffsetZ = dis(rng) * (float)chunkSize;
                int instanceId = (int)std::round(dis(rng) * (float)(numPoiInstances - 1));

                if (noiseValue < poiRate) {
                    float ax = (float)chunkMinX + chunkOffsetX;
                    float az = (float)chunkMinZ + chunkOffsetZ;
                    const float height = heightfieldSampler.getHeight(ax, az) - (float)WORLD_BASE_HEIGHT;

                    poiGeometry.ps.push_back(ax);
                    poiGeometry.ps.push_back(height);
                    poiGeometry.ps.push_back(az);

                    poiGeometry.instances.push_back(instanceId);
                }
            }
        }
    }
}

//

/* void generateBarrierGeometry(
    const vm::ivec2 &worldPosition,
    int lod,
    int chunkSize,
    OctreeContext &octreeContext,
    BarrierGeometry &geometry
) {
    generateBarrierMesh(worldPosition, lod, chunkSize, octreeContext, geometry);
    // offsetGeometry(geometry, worldPosition);
    // computeVertexNormals(geometry.positions, geometry.normals, geometry.indices);
} */

//

/* uint8_t getMostCommonBiome(const std::vector<Heightfield> &heightfields) {
    std::unordered_map<unsigned char, unsigned int> biomeCounts(numBiomes);
    for (const auto &hf : heightfields) {
        biomeCounts[hf.biomesVectorField[0]]++;
    }

    std::vector<unsigned char> seenBiomes;
    for (auto kv : biomeCounts)
    {
        unsigned char b = kv.first;
        seenBiomes.push_back(b);
    }
    
    // sort by increasing occurence count of the biome
    auto iter = std::max_element(
        seenBiomes.begin(),
        seenBiomes.end(),
        [&](unsigned char b1, unsigned char b2) -> bool {
            return biomeCounts[b1] > biomeCounts[b2];
        }
    );
    return *iter;
} */
inline int getLodInt(int lod) {
    return 1 << (lod - 1);
}
inline int getLodRange(int lod, int chunkSize) {
    return getLodInt(lod) * chunkSize;
}
OctreeContext PGInstance::getChunkSeedOctree(
    const vm::ivec2 &worldPosition,
    // int lod,
    int minLod,
    int maxLod, // we will sample a 3x3 of this lod
    int chunkSize
) {
    // constexpr int minLod = 1;
    // constexpr int maxLod = 6; // we will sample a 3x3 of this lod

    const int maxLodInt = getLodInt(maxLod);
    const int maxLodRange = getLodRange(maxLod, chunkSize);
    /* const int maxLodP1 = maxLod + 1;
    const int maxLodP1Range = (1 << (maxLodP1 - 1)) * chunkSize;
    vm::ivec2 maxLodP1Center{
        (int)std::floor((float)worldPosition.x / (float)maxLodRange) * maxLodRange,
        (int)std::floor((float)worldPosition.y / (float)maxLodRange) * maxLodRange
    }; */

    constexpr int maxNumSplits = 3;

    std::vector<vm::ivec2> maxLodChunkPositions;
    for (int dz = -1; dz <= 1; dz++) {
        for (int dx = -1; dx <= 1; dx++) {
            vm::ivec2 baseNode{
                (int)std::floor(
                    (float)(((float)worldPosition.x) / (float)maxLodRange) + (float)dx
                ) * maxLodInt,
                (int)std::floor(
                    (float)(((float)worldPosition.y) / (float)maxLodRange) + (float)dz
                ) * maxLodInt
            };

            // insert the node if it does not exist
            auto iter = std::find(
                maxLodChunkPositions.begin(),
                maxLodChunkPositions.end(),
                baseNode
            );
            if (iter == maxLodChunkPositions.end()) {
                maxLodChunkPositions.push_back(baseNode);
            } else {
                std::cerr << "ERROR: duplicate node found: " <<
                    baseNode.x << ", " << baseNode.y << " : " << dx << " " << dz <<
                    std::endl;
                abort();
            }
        }
    }

    // compute splits
    std::vector<std::pair<vm::ivec2, int>> lodSplits;
    for (size_t i = 0; i < maxLodChunkPositions.size(); i++) {
        const vm::ivec2 &baseNode = maxLodChunkPositions[i];

        float chunkSeed = noises.numSplitsNoise.in2D(baseNode.x, baseNode.y);
        unsigned int seedInt = *(unsigned int *)&chunkSeed;
        std::mt19937 rng(seedInt);
        std::uniform_real_distribution<float> dis(0.f, 1.f);

        uint32_t numSplits = (uint32_t)(dis(rng) * (float)maxNumSplits);
        for (uint32_t i = 0; i < numSplits; i++) {
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
                        splitLodDZInt
                    },
                    splitLodInt
                )
            );
        }
    }

    /* {
        float splitsLodNoise = (float)noises.numSplitsNoise.in2D(maxLodCenter.x, maxLodCenter.y);
        for (uint32_t i = 0; i < numSplitsHash; i++) {
            uint32_t splitLodDX;
            MurmurHash3_x86_32(&splitsLodNoise, sizeof(splitsLodNoise), 0, &splitLodDX);
            splitLodDX = (uint32_t)((float)splitLodDX * (float)maxLodRange / (float)0xFFFFFFFFu);
            splitsLodNoise++;

            uint32_t splitLodDZ;
            MurmurHash3_x86_32(&splitsLodNoise, sizeof(splitsLodNoise), 1, &splitLodDZ);
            splitLodDZ = (uint32_t)((float)splitLodDZ * (float)maxLodRange / (float)0xFFFFFFFFu);
            splitsLodNoise++;

            uint32_t splitLod;
            MurmurHash3_x86_32(&splitsLodNoise, sizeof(splitsLodNoise), 1, &splitLod);
            splitLod = (uint32_t)((float)minLod + (float)splitLod * (float)(maxLod - minLod) / (float)0xFFFFFFFFu);
            splitsLodNoise++;

            int splitLodDXInt = maxLodCenter.x + (int)splitLodDX;
            int splitLodDZInt = maxLodCenter.y + (int)splitLodDZ;
            int splitLodInt = 1 << ((int)splitLod - 1);

            lodSplits.push_back(
                std::make_pair(
                    vm::ivec2{
                        splitLodDXInt,
                        splitLodDZInt
                    },
                    splitLodInt
                )
            );
        }
    } */

    OctreeContext octreeContext;
    constructSeedTree(
        octreeContext,
        maxLodChunkPositions,
        maxLodInt,
        lodSplits
    );
    return octreeContext;
}
std::vector<Heightfield> PGInstance::getHeightfields(
    int x,
    int z,
    int lod,
    const std::array<int, 2> &lodArray
) {
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = chunkSize * lod / bottomLod;
    const int gridWidthP1 = gridWidth + 1;
    const int gridWidthP3 = gridWidth + 3;
    const int gridWidthP3T3 = gridWidthP3 * 3;
    
    const int gridHeight = chunkSize * lod / rightLod;
    const int gridHeightP1 = gridHeight + 1;
    const int gridHeightP3 = gridHeight + 3;
    const int gridHeightP3T3 = gridHeightP3 * 3;

    const int chunkSizeP2 = chunkSize + 2;

    std::vector<Heightfield> heightfields(
        (chunkSizeP2 * chunkSizeP2) + // center
        (gridWidthP3T3 + gridHeightP3T3) // seams
    );

    getHeightFieldCenter(x, z, lod, heightfields);
    getHeightFieldSeams(x, z, lod, lodArray, chunkSizeP2, heightfields);

    return heightfields;
}
enum GenerateFlags {
    GF_NONE = 0,
    GF_TERRAIN = 1 << 0,
    GF_WATER = 1 << 1,
    GF_VEGETATION = 1 << 2,
    GF_GRASS = 1 << 3,
    GF_POI = 1 << 4
};
ChunkResult *PGInstance::createChunkMesh(
    const vm::ivec2 &worldPosition,
    int lod,
    const std::array<int, 2> &lodArray,
    int generateFlags,
    int numVegetationInstances,
    int numGrassInstances,
    int numPoiInstances
) {
    ChunkResult *result = (ChunkResult *)malloc(sizeof(ChunkResult));

    // heightfield
    std::vector<Heightfield> heightfields;
    if (
        (generateFlags & GF_TERRAIN) |
        (generateFlags & GF_WATER) |
        (generateFlags & GF_VEGETATION) |
        (generateFlags & GF_GRASS) |
        (generateFlags & GF_POI)
    ) {
        heightfields = getHeightfields(worldPosition.x, worldPosition.y, lod, lodArray);
    }

    // terrain
    if (generateFlags & GF_TERRAIN) {
        TerrainGeometry terrainGeometry;

        generateTerrainGeometry(
            worldPosition,
            lod,
            lodArray,
            chunkSize,
            heightfields,
            terrainGeometry
        );
        result->terrainMeshBuffer = terrainGeometry.getBuffer();
    } else {
        result->terrainMeshBuffer = nullptr;
    }

    // water
    if (generateFlags & GF_WATER) {
        WaterGeometry waterGeometry;

        const std::vector<Waterfield> &waterfields = *((std::vector<Waterfield> *)&heightfields);

        generateWaterGeometry(
            worldPosition,
            lod,
            lodArray,
            chunkSize,
            waterfields,
            waterGeometry
        );
        result->waterMeshBuffer = waterGeometry.getBuffer();
    } else {
        result->waterMeshBuffer = nullptr;
    }

    // vegetation
    if (generateFlags & GF_VEGETATION) {
        VegetationGeometry vegetationGeometry;

        generateVegetationInstances(
            worldPosition,
            lod,
            chunkSize,
            numVegetationInstances,
            heightfields,
            noises,
            vegetationGeometry
        );
        result->vegetationInstancesBuffer = vegetationGeometry.getBuffer();
    } else {
        result->vegetationInstancesBuffer = nullptr;
    }

    // grass
    if (generateFlags & GF_GRASS) {
        GrassGeometry grassGeometry;

        generateGrassInstances(
            worldPosition,
            lod,
            chunkSize,
            numGrassInstances,
            heightfields,
            noises,
            grassGeometry
        );
        result->grassInstancesBuffer = grassGeometry.getBuffer();
    } else {
        result->grassInstancesBuffer = nullptr;
    }

    // poi
    if (generateFlags & GF_POI) {
        PoiGeometry poiGeometry;

        generatePoiInstances(
            worldPosition,
            lod,
            chunkSize,
            numPoiInstances,
            heightfields,
            noises,
            poiGeometry
        );
        result->poiInstancesBuffer = poiGeometry.getBuffer();
    } else {
        result->poiInstancesBuffer = nullptr;
    }

    return result;
}
void PGInstance::createMobSplatAsync(
    uint32_t id,
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int priority
) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPositionXZ.x,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPositionXZ.y
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

//

uint8_t *PGInstance::createBarrierMesh(
    const vm::ivec2 &worldPosition,
    int minLod,
    int maxLod
) {
    OctreeContext octreeContext = getChunkSeedOctree(
        worldPosition,
        minLod,
        maxLod,
        chunkSize
    );

    BarrierGeometry barrierGeometry;

    generateBarrierGeometry(
        worldPosition,
        chunkSize,
        octreeContext,
        barrierGeometry
    );

    uint8_t *result = barrierGeometry.getBuffer();
    return result;
}
void PGInstance::createBarrierMeshAsync(
    uint32_t id,
    const vm::ivec2 &worldPosition,
    int minLod,
    int maxLod
) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    // this is so the chunk center is roughly the center of the one used by createBarrierMesh
    const int maxLodRange = getLodRange(maxLod, chunkSize);
    vm::vec3 basePositionF{ // center of the middle 3x3 chunks
        std::floor(
            (float)(((float)worldPosition.x) / (float)maxLodRange)
        ) * maxLodRange + (maxLodRange / 2.f),
        (float)(-WORLD_BASE_HEIGHT) + ((float)MIN_WORLD_HEIGHT + (float)MAX_WORLD_HEIGHT) / 2.f,
        std::floor(
            (float)(((float)worldPosition.y + (float)maxLodRange / 2.f) / (float)maxLodRange)
        ) * maxLodRange + (maxLodRange / 2.f),
    };
    const int maxLodP1 = maxLod + 1;
    Task *terrainTask = new Task(id, basePositionF, maxLodP1, [
        this,
        promise,
        worldPosition,
        minLod,
        maxLod
    ]() -> void {
        uint8_t *result = createBarrierMesh(
            worldPosition,
            minLod,
            maxLod
        );
        if (!promise->resolve(result)) {
            free(result);
        }
    });
    ProcGen::taskQueue.pushTask(terrainTask);
}

/* uint8_t *PGInstance::createLiquidChunkMesh(const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray)
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
} */

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

void PGInstance::createChunkMeshAsync(
    uint32_t id,
    const vm::ivec2 &worldPosition,
    int lod,
    const std::array<int, 2> &lodArray,
    int generateFlags,
    int numVegetationInstances,
    int numGrassInstances,
    int numPoiInstances
) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x,
        (float)(-WORLD_BASE_HEIGHT) + ((float)MIN_WORLD_HEIGHT + (float)MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y
    };
    std::array<int, 2> lodArray2{
        lodArray[0],
        lodArray[1]
    };
    Task *terrainTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod,
        lodArray2,
        generateFlags,
        numVegetationInstances,
        numGrassInstances,
        numPoiInstances
    ]() -> void {
        ChunkResult *result = createChunkMesh(
            worldPosition,
            lod,
            lodArray2,
            generateFlags,
            numVegetationInstances,
            numGrassInstances,
            numPoiInstances
        );
        if (!promise->resolve(result)) {
            result->free();
        }
    });
    ProcGen::taskQueue.pushTask(terrainTask);
}
/* void PGInstance::createLiquidChunkMeshAsync(uint32_t id, const vm::ivec2 &worldPosition, int lod, const std::array<int, 2> &lodArray)
{
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x + (float)lod / 2.f,
        (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y + (float)lod / 2.f
    };
    std::array<int, 2> lodArray2{
        lodArray[0],
        lodArray[1]
    };
    Task *liquidTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod,
        lodArray2
    ]() -> void {
        uint8_t *result = createLiquidChunkMesh(worldPosition, lod, lodArray2);
        promise->resolve(result);
    });
    ProcGen::taskQueue.pushTask(liquidTask);
} */
/* void PGInstance::createChunkVegetationAsync(uint32_t id, const vm::ivec2 &worldPosition, const int lod, const int numInstances) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x,
        (float)(-WORLD_BASE_HEIGHT) + ((float)MIN_WORLD_HEIGHT + (float)MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y
    };
    Task *vegetationSplatTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod,
        numInstances
    ]() -> void {
        void *result = createChunkVegetation(worldPosition, lod, numInstances);
        if (!promise->resolve(result)) {
            free(result);
        }
    });
    ProcGen::taskQueue.pushTask(vegetationSplatTask);
} */

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

/* void PGInstance::createChunkGrassAsync(uint32_t id, const vm::ivec2 &worldPosition, const int lod, const int numInstances) {
    std::shared_ptr<Promise> promise = ProcGen::resultQueue.createPromise(id);

    vm::vec3 worldPositionF{
        (float)worldPosition.x,
        (float)(-WORLD_BASE_HEIGHT) + ((float)MIN_WORLD_HEIGHT + (float)MAX_WORLD_HEIGHT) / 2.f,
        (float)worldPosition.y
    };
    Task *grassSplatTask = new Task(id, worldPositionF, lod, [
        this,
        promise,
        worldPosition,
        lod,
        numInstances
    ]() -> void {
        std::array<int, 2> lodArray;
        lodArray[0] = lod;
        lodArray[1] = lod;
        std::vector<Heightfield> heightfields = getHeightfields(worldPosition.x, worldPosition.y, lod, lodArray);
        uint8_t *result = createChunkGrass(worldPosition, lod, heightfields, numInstances);
        if (!promise->resolve(result)) {
          free(result);
        }
    });
    ProcGen::taskQueue.pushTask(grassSplatTask);
} */

// 2d caches

NoiseField PGInstance::getNoise(float bx, float bz) {
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
}
uint8_t PGInstance::getBiome(float bx, float bz) {
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
}

//

void PGInstance::getHeightFieldCenter(int bx, int bz, int lod, std::vector<Heightfield> &heightfields) {
    const int chunkSizeP2 = chunkSize + 2;
    for (int dz = 0; dz < chunkSizeP2; dz++) {
        for (int dx = 0; dx < chunkSizeP2; dx++) {
            const int index = dx + dz * chunkSizeP2;
            Heightfield &localHeightfield = heightfields[index];
            
            const int x = dx - 1;
            const int z = dz - 1;
            const int ax = bx + x * lod;
            const int az = bz + z * lod;
            localHeightfield = getHeightField(ax, az);
        }
    }
}
void PGInstance::getHeightFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, int rowSize, std::vector<Heightfield> &heightfields) {
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = chunkSize * lod / bottomLod;
    const int gridWidthP3 = gridWidth + 3;

    const int gridHeight = chunkSize * lod / rightLod;
    const int gridHeightP3 = gridHeight + 3;

    const int heightfieldsCenterDataOffset = rowSize * rowSize;

    // bottom
    int index = heightfieldsCenterDataOffset;
    {
        for (int dz = 0; dz < 3; dz++) {
            for (int dx = 0; dx < gridWidthP3; dx++) {
                const int x = dx - 1;
                const int z = dz - 1;

                Heightfield &localHeightfieldSeam = heightfields[index];
                
                const int ax = bx + x * bottomLod;
                const int az = bz + chunkSize * lod + z * bottomLod;
                localHeightfieldSeam = getHeightField(ax, az);

                index++;
            }
        }
    }
    // right
    {
        for (int dx = 0; dx < 3; dx++) {
            for (int dz = 0; dz < gridHeightP3; dz++) {
                const int x = dx - 1;
                const int z = dz - 1;
                
                Heightfield &localHeightfieldSeam = heightfields[index];

                const int ax = bx + chunkSize * lod + x * rightLod;
                const int az = bz + z * rightLod;
                localHeightfieldSeam = getHeightField(ax, az);

                index++;
            }
        }
    }
}

Heightfield PGInstance::getHeightField(float bx, float bz) {
    Heightfield localHeightfield;

    vm::vec2 fWorldPosition{bx, bz};

    // height
    const float halfChunkSizeF = (float)chunkSize / 2.f;
    const float maxDistance = std::sqrt(halfChunkSizeF);

    // water
    constexpr int waterRange = 4;
    const float maxWaterDistance = (float)std::sqrt((float)waterRange * (float)waterRange);
    constexpr float baseWaterFactor = 0.25;

    // acc height
    std::vector<float> biomeCounts(numBiomes);
    float totalHeightFactors = 0;
    // acc material
    std::vector<MaterialWeightAccumulator> materialWeightAccumulators(numMaterials);
    float totalMaterialFactors = 0;
    // acc water
    float sumWaterFactor = 0;
    float totalWaterFactors = 0;
    // loop
    for (float dz = -halfChunkSizeF; dz <= halfChunkSizeF; dz++) {
        for (float dx = -halfChunkSizeF; dx <= halfChunkSizeF; dx++) {
            float distance = std::sqrt(dx*dx + dz*dz);

            float ax = bx + dx;
            float az = bz + dz;
            unsigned char b = 0x0;
            bool hasBiome = false;

            if (distance < maxDistance) {
                if (!hasBiome) {
                    b = getBiome(ax, az);
                    hasBiome = true;
                }

                float heightFactor = 1.f - (distance / maxDistance);

                biomeCounts[b] += heightFactor;
                totalHeightFactors += heightFactor;
            }

            if (distance < maxWaterDistance) {
                if (!hasBiome) {
                    b = getBiome(ax, az);
                    hasBiome = true;
                }

                float waterFactor = 1.f - (distance / maxWaterDistance);
                waterFactor =
                    baseWaterFactor +
                    (1.f - baseWaterFactor) * waterFactor;

                if (isWaterBiome(b)) {
                    sumWaterFactor += waterFactor;
                }
                totalWaterFactors += waterFactor;
            }
        }
    }
    sumWaterFactor /= totalWaterFactors;

    // postprocess height
    {
        std::vector<uint8_t> seenBiomes;
        for (size_t i = 0; i < biomeCounts.size(); i++)
        {
            const uint8_t biome = (uint8_t)i;
            const float &biomeWeight = biomeCounts[biome];
            if (biomeWeight > 0.f) {
                seenBiomes.push_back(biome);
            }
        }
        
        // sort by increasing occurence count of the biome
        std::sort(
            seenBiomes.begin(),
            seenBiomes.end(),
            [&](uint8_t b1, uint8_t b2) -> bool {
                return biomeCounts[b1] > biomeCounts[b2];
            }
        );

        // letting the biome weight fit in an unsigned char
        const float biomeWeightFitter = 255.f;

        for (size_t i = 0; i < 4; i++)
        {
            if (i < seenBiomes.size())
            {
                const uint8_t &biome = seenBiomes[i];
                localHeightfield.biomesVectorField[i] = biome;
                localHeightfield.biomesWeightsVectorField[i] = biomeCounts[biome] / totalHeightFactors * biomeWeightFitter;
            }
            else
            {
                localHeightfield.biomesVectorField[i] = 0;
                localHeightfield.biomesWeightsVectorField[i] = 0;
            }
        }

        float elevationSum = 0.f;
        vm::vec2 fWorldPosition{bx, bz};
        for (size_t i = 0; i < biomeCounts.size(); i++)
        {
            const uint8_t biome = (uint8_t)i;
            const float &biomeWeight = biomeCounts[biome];
            if (biomeWeight > 0.f) {
                elevationSum += biomeWeight * getComputedBiomeHeight(biome, fWorldPosition);
            }
        }

        float elevation = elevationSum / totalHeightFactors;
        localHeightfield.heightField = elevation;

        // materials

        getComputedMaterials(localHeightfield, materialWeightAccumulators, totalMaterialFactors, fWorldPosition);

        std::vector<uint8_t> seenMaterials;
        for (size_t i = 0; i < materialWeightAccumulators.size(); i++)
        {
            const uint8_t material = (uint8_t)i;
            const bool &isSeen = materialWeightAccumulators[material].getSeen();
            if(isSeen) {
                seenMaterials.push_back(material);
            }
        }
        
        // sort by the overall weight of the material
        // std::sort(
        //     seenMaterials.begin(),
        //     seenMaterials.end(),
        //     [&](uint8_t b1, uint8_t b2) -> bool {
        //         return materialWeightAccumulators[b1].getWeight() > materialWeightAccumulators[b2].getWeight();
        //     }
        // );

        for (size_t i = 0; i < 4; i++)
        {
            if (i < seenMaterials.size())
            {
                const uint8_t &material = seenMaterials[i];
                const float &materialWeight = materialWeightAccumulators[material].getWeight();
                localHeightfield.materials[i] = material;
                localHeightfield.materialsWeights[i] = materialWeight / totalMaterialFactors;
            }
            else
            {
                localHeightfield.materials[i] = 0;
                localHeightfield.materialsWeights[i] = 0;
            }
        }
    }

    // postprocess water
    {
        localHeightfield.waterFactor = sumWaterFactor;
    }


    return localHeightfield;
}
float PGInstance::getHeight(float bx, float bz) {
    const float halfChunkSizeF = (float)chunkSize / 2.f;
    const float maxDistance = std::sqrt(halfChunkSizeF + 1.f);

    std::unordered_map<unsigned char, float> biomeCounts(numBiomes);
    float totalHeightFactors = 0;
    for (float dz = -halfChunkSizeF; dz <= halfChunkSizeF; dz++)
    {
        for (float dx = -halfChunkSizeF; dx <= halfChunkSizeF; dx++)
        {
            float distance = std::sqrt(dx*dx + dz*dz);
            float factor = std::max(1.f - (distance / maxDistance), 0.f);

            if (factor > 0) {
                float ax = bx + dx;
                float az = bz + dz;
                unsigned char b = getBiome(ax, az);

                biomeCounts[b] += factor;
                totalHeightFactors += factor;
            }
        }
    }

    float elevationSum = 0.f;
    vm::vec2 fWorldPosition{bx, bz};
    for (auto const &iter : biomeCounts)
    {
        elevationSum += iter.second * getComputedBiomeHeight(iter.first, fWorldPosition);
    }

    float elevation = elevationSum / totalHeightFactors;
    return elevation;
}

//

/* void PGInstance::getWaterFieldCenter(int bx, int bz, int lod, std::vector<Waterfield> &waterfields) {
    for (int z = 0; z < chunkSize; z++)
    {
        for (int x = 0; x < chunkSize; x++)
        {
            const int index = x + z * chunkSize;
            Waterfield &localWaterfield = waterfields[index];
            localWaterfield = getWaterField(bx + x * lod, bz + z * lod, lod);
        }
    }
}
void PGInstance::getWaterFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, int rowSize, std::vector<Waterfield> &waterfields) {
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = chunkSize * lod / bottomLod;
    const int gridWidthP1 = gridWidth + 1;

    const int gridHeight = chunkSize * lod / rightLod;
    const int gridHeightP1 = gridHeight + 1;

    const int waterfieldsCenterDataOffset = rowSize * rowSize;

    // bottom
    int index = waterfieldsCenterDataOffset;
    {
        for (int dz = 0; dz < 1; dz++) {
            for (int dx = 0; dx < gridWidthP1; dx++) {
                const int z = dz;
                const int x = dx;

                Waterfield &localWaterfieldSeam = waterfields[index];
                
                const int ax = bx + x * bottomLod;
                const int az = bz + chunkSize * lod + z * bottomLod;
                localWaterfieldSeam = getWaterField(ax, az, lod);

                index++;
            }
        }
    }
    // right
    {
        for (int dx = 0; dx < 1; dx++) {
            for (int dz = 0; dz < gridHeightP1; dz++) {
                const int z = dz;
                const int x = dx;
                
                Waterfield &localWaterfieldSeam = waterfields[index];

                const int ax = bx + chunkSize * lod + x * rightLod;
                const int az = bz + z * rightLod;
                localWaterfieldSeam = getWaterField(ax, az, lod);

                index++;
            }
        }
    }
} */

//

/* Waterfield PGInstance::getWaterField(int bx, int bz, int lod) {
    constexpr int waterRange = 4;
    const float maxWaterDistance = (float)std::sqrt((float)waterRange * (float)waterRange);
    constexpr float baseWaterFactor = 0.25;

    float waterFactor = 0;
    float maxWaterFactor = 0;
    for (int dz = -waterRange; dz <= waterRange; dz++) {
        const int az = bz + dz * lod;

        for (int dx = -waterRange; dx <= waterRange; dx++) {
            const int ax = bx + dx * lod;

            const float distance = (float)std::sqrt((float)dx * (float)dx + (float)dz * (float)dz);
            const float localFactor =
              baseWaterFactor +
              (1.f - baseWaterFactor) * (1. - distance / maxWaterDistance);

            unsigned char b = getBiome(ax, az);
            if (isWaterBiome(b)) {
                waterFactor += localFactor;
            }
            maxWaterFactor += localFactor;
        }
    }
    waterFactor /= maxWaterFactor;

    return Waterfield{
        waterFactor
    };
} */

//

/* float PGInstance::getSeed(int x, int z) {
    constexpr int maxLod = 7;
    constexpr int maxLodRange = 1 << (maxLod - 1);
    vm::ivec2 basePosition{
        (x / maxLodRange) * maxLodRange,
        (z / maxLodRange) * maxLodRange
    };
    const SeedNoise &seedNoise = getSeedNoise(
        basePosition.x,
        basePosition.y
    );
    const float seedNoiseFloat = (float)std::abs(seedNoise.seed);
    uint32_t seedNoiseHash;
    MurmurHash3_x86_32(&seedNoiseFloat, sizeof(seedNoiseFloat), 0, &seedNoiseHash);
    const int seedLod = (int)((float)seedNoiseHash / (float)0xFFFFFFFF * (float)maxLod);
    const int seedLodRange = 1 << seedLod;
    return seedNoiseFloat;
} */

//

// 3d caches

/* uint8_t PGInstance::initSkylightField(PGInstance *inst, int x, int y, int z) {
    constexpr uint8_t maxSkyLightu8 = 8;
    constexpr int maxSkyLighti = (int)maxSkyLightu8;
    constexpr float maxSkyLightf = (float)maxSkyLighti;

    skylight = std::min(std::max(skylight, (uint8_t)0), maxSkyLightu8);

    // should flood fill the light

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
    const float &ax = worldPosition.x;
    const float &az = worldPosition.y;

    float biomeHeight = biome.baseHeight +
        noises.elevationNoise1.in2D(ax * biome.amps[0][0], az * biome.amps[0][0]) * biome.amps[0][1] +
        noises.elevationNoise2.in2D(ax * biome.amps[1][0], az * biome.amps[1][0]) * biome.amps[1][1] +
        noises.elevationNoise3.in2D(ax * biome.amps[2][0], az * biome.amps[2][0]) * biome.amps[2][1];
    return biomeHeight;
}

// materials
void PGInstance::getComputedMaterials(Heightfield &localHeightfield, std::vector<MaterialWeightAccumulator> &materialWeightAccumulators, float &totalMaterialFactors, const vm::vec2 &worldPosition)
{
    const std::array<uint8_t, 4> &biomes = localHeightfield.biomesVectorField;
    const std::array<uint8_t, 4> &biomesWeights = localHeightfield.biomesWeightsVectorField;

    const int GRASS = (int)MATERIAL::GRASS;
    const int DIRT = (int)MATERIAL::DIRT;

    for (int i = 0; i < 1; i++)
    {
        const uint8_t b = biomes[i];
        const float bw = (float)biomesWeights[i] / 255.f;
        switch (b)
        {
        // TODO : Define a different set of material rules for each biome, for now we're using these rules as default
        default:
            const float materialNoise = vm::clamp(noises.grassMaterialNoise.in2DWarp(worldPosition), 0.f, 1.f);

            const float grassWeight = (1.f - materialNoise) * bw;
            const float dirtWeight = (materialNoise) * bw;

            MaterialWeightAccumulator &grassWeightAcc = materialWeightAccumulators[GRASS];
            MaterialWeightAccumulator &dirtWeightAcc = materialWeightAccumulators[DIRT];

            grassWeightAcc.addWeight(grassWeight);
            dirtWeightAcc.addWeight(dirtWeight);

            totalMaterialFactors += grassWeight;
            totalMaterialFactors += dirtWeight;

            break;
        }
    }
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
          free(buffer);
        }
    });
    ProcGen::taskQueue.pushTask(trackerUpdateTask);
}