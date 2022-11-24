#include "instance.h"
#include "procgen.h"
#include "octree.h"
#include "biomes.h"
#include "tracker.h"
#include "vector.h"
#include "util.h"
#include "MurmurHash3.h"
#include <emscripten.h>
// #include "peek.h"

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
PGInstance::PGInstance(int seed, int chunkSize) : generator(seed, chunkSize)

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

    const int chunkSize = generator.getChunkSize();

    int minX = worldPositionXZ.x / chunkSize * chunkSize;
    int minZ = worldPositionXZ.y / chunkSize * chunkSize;

    float seed = generator.noises.mobNoise.in2D(minX, minZ);
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

        float noiseValue = generator.noises.mobNoise.in2D(ax, az);

        if (noiseValue < mobRate)
        {
            const float height = generator.getHeight(ax, az);

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

void normalizeNormals(std::vector<vm::vec3> &normals)
{
    for (size_t i = 0, il = normals.size(); i < il; i++)
    {
        Vec vec{
            normals[i].x,
            normals[i].y,
            normals[i].z};
        vec.normalize();
        // normals.setXYZ(i, vec.x, vec.y, vec.z);
        normals[i].x = vec.x;
        normals[i].y = vec.y;
        normals[i].z = vec.z;
    }
}
void computeFaceNormals(std::vector<vm::vec3> &positions, std::vector<vm::vec3> &normals, std::vector<uint32_t> &indices)
{
    // indexed elements
    for (size_t i = 0, il = indices.size(); i < il; i += 3)
    {
        const uint32_t &vA = indices[i];
        const uint32_t &vB = indices[i + 1];
        const uint32_t &vC = indices[i + 2];

        Vec pA{
            positions[vA].x,
            positions[vA].y,
            positions[vA].z};
        Vec pB{
            positions[vB].x,
            positions[vB].y,
            positions[vB].z};
        Vec pC{
            positions[vC].x,
            positions[vC].y,
            positions[vC].z};

        Vec cb = pC - pB;
        Vec ab = pA - pB;
        cb ^= ab;

        Vec nA{
            normals[vA].x,
            normals[vA].y,
            normals[vA].z};
        Vec nB{
            normals[vB].x,
            normals[vB].y,
            normals[vB].z};
        Vec nC{
            normals[vC].x,
            normals[vC].y,
            normals[vC].z};

        nA += cb;
        nB += cb;
        nC += cb;

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

template <typename T>
vm::vec3 calculateCenterPointNormal(const int &x, const int &y, const std::vector<T> &heightfields, const int &rowSize)
{
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

    return vm::normalize(vm::vec3{Lheight - Rheight, 2.0f, Uheight - Dheight});
}

template <typename T>
vm::vec3 calculateBottomPointNormal(const int &x, const int &y, const std::vector<T> &heightfields, const int &heightfieldsCenterDataReadOffset, const int &gridWidth)
{
    const int gridWidthP3 = gridWidth + 3;
    const int Lx = x - 1;
    const int Ly = y;
    const int Lindex =
        heightfieldsCenterDataReadOffset +
        (Ly * gridWidthP3) +
        Lx;
    const float Lheight = heightfields[Lindex].getHeight();

    const int Rx = x + 1;
    const int Ry = y;
    const int Rindex =
        heightfieldsCenterDataReadOffset +
        (Ry * gridWidthP3) +
        Rx;
    const float Rheight = heightfields[Rindex].getHeight();

    const int Ux = x;
    const int Uy = y - 1;
    const int Uindex =
        heightfieldsCenterDataReadOffset +
        (Uy * gridWidthP3) +
        Ux;
    const float Uheight = heightfields[Uindex].getHeight();

    const int Dx = x;
    const int Dy = y + 1;
    const int Dindex =
        heightfieldsCenterDataReadOffset +
        (Dy * gridWidthP3) +
        Dx;
    const float Dheight = heightfields[Dindex].getHeight();

    return vm::normalize(vm::vec3{Lheight - Rheight, 2.0f, Uheight - Dheight});
}

template <typename T>
vm::vec3 calculateRightPointNormal(const int &x, const int &y, const std::vector<T> &heightfields, const int &heightfieldsCenterDataReadOffset, const int &gridWidth, const int &gridHeight)
{
    const int gridWidthP3 = gridWidth + 3;
    const int gridHeightP3 = gridHeight + 3;

    const int Lx = x - 1;
    const int Ly = y;
    const int Lindex =
        heightfieldsCenterDataReadOffset +
        (3 * gridWidthP3) +
        (Lx * gridHeightP3) +
        Ly;
    const float Lheight = heightfields[Lindex].getHeight();

    const int Rx = x + 1;
    const int Ry = y;
    const int Rindex =
        heightfieldsCenterDataReadOffset +
        (3 * gridWidthP3) +
        (Rx * gridHeightP3) +
        Ry;
    const float Rheight = heightfields[Rindex].getHeight();

    const int Ux = x;
    const int Uy = y - 1;
    const int Uindex =
        heightfieldsCenterDataReadOffset +
        (3 * gridWidthP3) +
        (Ux * gridHeightP3) +
        Uy;
    const float Uheight = heightfields[Uindex].getHeight();

    const int Dx = x;
    const int Dy = y + 1;
    const int Dindex =
        heightfieldsCenterDataReadOffset +
        (3 * gridWidthP3) +
        (Dx * gridHeightP3) +
        Dy;
    const float Dheight = heightfields[Dindex].getHeight();

    return vm::normalize(vm::vec3{Lheight - Rheight, 2.0f, Uheight - Dheight});
}

template <typename T>
void calculateCenterNormals(std::vector<T> &heightfields, const int &chunkSize, const int &rowSize)
{
    auto pushCenterPointNormal = [&](const int &x, const int &y) -> void
    {
        const int dx = x + 1;
        const int dy = y + 1;
        const int index = dx + dy * rowSize;
        T &v = heightfields[index];
        v.normal = calculateCenterPointNormal<T>(x, y, heightfields, rowSize);
    };

    for (int y = 0; y < chunkSize; y++)
    {
        for (int x = 0; x < chunkSize; x++)
        {
            pushCenterPointNormal(x, y);
        }
    }
}

template <typename T>
void calculateSeamNormals(std::vector<T> &heightfields,
                          const int &lod,
                          const std::array<int, 2> &lodArray,
                          const int &chunkSize,
                          const int &rowSize)
{
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = chunkSize * lod / bottomLod;
    const int gridWidthP1 = gridWidth + 1;
    const int gridWidthP3 = gridWidth + 3;

    const int gridHeight = chunkSize * lod / rightLod;
    const int gridHeightP1 = gridHeight + 1;
    const int gridHeightP3 = gridHeight + 3;

    const int heightfieldsCenterDataReadOffset = rowSize * rowSize;

    auto pushBottomPointNormal = [&](const int &x, const int &y) -> void
    {
        const int dx = x + 1;
        const int dy = y + 1;
        const int index =
            heightfieldsCenterDataReadOffset +
            (dy * gridWidthP3) +
            dx;
        T &v = heightfields[index];

        v.normal = calculateBottomPointNormal(dx, dy, heightfields, heightfieldsCenterDataReadOffset, gridWidth);
    };

    // bottom
    {
        const int y = 0;
        for (int x = 0; x < gridWidthP1; x++)
        {
            pushBottomPointNormal(x, y);
        }
    }

    auto pushRightPointNormal = [&](const int &x, const int &y) -> void
    {
        const int dx = x + 1;
        const int dy = y + 1;
        const int index =
            heightfieldsCenterDataReadOffset +
            (3 * gridWidthP3) +
            (dx * gridHeightP3) +
            dy;
        T &v = heightfields[index];

        v.normal = calculateRightPointNormal(dx, dy, heightfields, heightfieldsCenterDataReadOffset, gridWidth, gridHeight);
    };

    // right
    {
        const int x = 0;
        for (int y = 0; y < gridHeightP1; y++)
        {
            pushRightPointNormal(x, y);
        }
    }
}

template <typename T>
void calculateSurfaceNormals(std::vector<T> &heightfields,
                             const int &lod,
                             const std::array<int, 2> &lodArray,
                             const int &chunkSize)
{

    const int rowSize = chunkSize + 2;

    calculateCenterNormals<T>(heightfields, chunkSize, rowSize);
    calculateSeamNormals<T>(heightfields, lod, lodArray, chunkSize, rowSize);
}
/* void fillVec3(std::vector<vm::vec3> &array, const vm::vec3 &v) {
    for (size_t i = 0, il = array.size(); i < il; i++) {
        array[i] = v;
    }
} */

//

enum class WindingDirection
{
    CCW,
    CW
};
enum class ComputeNormals
{
    NO,
    YES
};

template <
    typename T,
    typename G,
    WindingDirection windingDirection,
    ComputeNormals computeNormals>
void createPlaneGeometry(
    int width,
    int height,
    int widthSegments,
    int heightSegments,
    int rowSize,
    const std::vector<T> &heightfields,
    G &geometry)
{
    const int &gridX = widthSegments;  // equals chunkSize - 1
    const int &gridY = heightSegments; // equals chunkSize - 1

    const int gridX1 = gridX + 1; // equals chunkSize
    const int gridY1 = gridY + 1; // equals chunkSize

    const int segment_width = width / gridX;   // equals lod
    const int segment_height = height / gridY; // equals lod

    //

    auto pushPoint = [&](int x, int y) -> void
    {
        const int dx = x + 1;
        const int dy = y + 1;
        const int index = dx + dy * rowSize;
        const T &v0 = heightfields[index];

        const int ax = x * segment_width;
        const int ay = y * segment_height;

        // position
        const float height = v0.getHeight();

        geometry.positions.push_back(vm::vec3{
            (float)ax,
            height,
            (float)ay});

        // materials
        const MaterialsArray &materials = v0.materials;
        const MaterialsWeightsArray &materialWeights = v0.materialsWeights;

        geometry.materials.push_back(vm::ivec4{
            materials[0],
            materials[1],
            materials[2],
            materials[3]});
        geometry.materialsWeights.push_back(vm::vec4{
            materialWeights[0],
            materialWeights[1],
            materialWeights[2],
            materialWeights[3]});

        // normal
        const vm::vec3 &normal = v0.normal;
        geometry.normals.push_back(normal);

        // metadata
        geometry.pushPointMetadata(v0);
    };

    // positions
    for (int y = 0; y < gridY1; y++)
    {
        for (int x = 0; x < gridX1; x++)
        {
            pushPoint(x, y);
        }
    }

    auto pushTriangle = [&](int ra, int rb, int rc, int wa, int wb, int wc) -> void
    {
        const T &hfA = heightfields[ra];
        const T &hfB = heightfields[rb];
        const T &hfC = heightfields[rc];
        if (T::acceptIndices(hfA, hfB, hfC))
        {
            if (windingDirection == WindingDirection::CCW)
            {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wb);
                geometry.indices.push_back(wc);
            }
            else
            {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wc);
                geometry.indices.push_back(wb);
            }
        }
    };

    // indices
    for (int y = 0; y < gridY; y++)
    {
        for (int x = 0; x < gridX; x++)
        {
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
template <
    typename T,
    typename G,
    WindingDirection windingDirection,
    ComputeNormals computeNormals>
void createPlaneSeamsGeometry(
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    int rowSize,
    const std::vector<T> &heightfields,
    G &geometry)
{
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

    auto pushBottomPoint = [&](int x, int y) -> void
    {
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
            (float)ay});

        geometry.materials.push_back(vm::ivec4{
            materials[0],
            materials[1],
            materials[2],
            materials[3]});

        geometry.materialsWeights.push_back(vm::vec4{
            materialWeights[0],
            materialWeights[1],
            materialWeights[2],
            materialWeights[3]});

        // normal
        const vm::vec3 &normal = v0.normal;
        geometry.normals.push_back(normal);

        // metadata
        geometry.pushPointMetadata(v0);
    };
    auto pushRightPoint = [&](int x, int y) -> void
    {
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
            (float)ay});

        geometry.materials.push_back(vm::ivec4{
            materials[0],
            materials[1],
            materials[2],
            materials[3]});

        geometry.materialsWeights.push_back(vm::vec4{
            materialWeights[0],
            materialWeights[1],
            materialWeights[2],
            materialWeights[3]});

        // normal
        const vm::vec3 &normal = v0.normal;
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
        for (int x = 0; x < gridWidthP1; x++)
        {
            pushBottomPoint(x, y);
        }
    }
    // right
    {
        const int x = 0;
        for (int y = 0; y < gridHeightP1; y++)
        {
            pushRightPoint(x, y);
        }
    }

    //

    const int chunkSizeM1 = chunkSize - 1;
    const int gridX = chunkSizeM1;
    const int gridY = chunkSizeM1;

    const int gridX1 = gridX + 1; // equals chunkSize
    const int gridY1 = gridY + 1; // equals chunkSize

    auto pushTriangle = [&](int ra, int rb, int rc, int wa, int wb, int wc)
    {
        const T &hfA = heightfields[ra];
        const T &hfB = heightfields[rb];
        const T &hfC = heightfields[rc];
        if (T::acceptIndices(hfA, hfB, hfC))
        {
            if (windingDirection == WindingDirection::CCW)
            {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wb);
                geometry.indices.push_back(wc);
            }
            else
            {
                geometry.indices.push_back(wa);
                geometry.indices.push_back(wc);
                geometry.indices.push_back(wb);
            }
        }
    };

    // indices
    // bottom
    if (bottomLod == lod)
    {
        const int innerPointY = chunkSize - 1;
        for (int innerPointX = 0; innerPointX < chunkSize; innerPointX++)
        {
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
            if (innerPointX != (chunkSize - 1))
            { // only single triangle in corner
                pushTriangle(ra, rc, rd, wa, wc, wd);
            }
        }
    }
    else if (bottomLod > lod)
    {
        const int innerPointY = chunkSize - 1;
        for (int innerPointX = 0; innerPointX < chunkSize; innerPointX++)
        {
            const int outerPointX = innerPointX / 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            const int outerPointDx = outerPointX + 1;

            if (innerPointX % 2 == 0)
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
            else
            {
                if (innerPointX != (chunkSize - 1))
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

                    pushTriangle(ra, rc, rd, wa, wc, wd);
                }
            }
        }
    }
    else if (bottomLod < lod)
    {
        const int innerPointY = chunkSize - 1;
        for (int innerPointX = 0; innerPointX < chunkSize; innerPointX++)
        {
            int outerPointX = innerPointX * 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            int outerPointDx = outerPointX + 1;

            if (innerPointX != (chunkSize - 1))
            {
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
            }
            else
            {
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
    if (rightLod == lod)
    {
        const int innerPointX = chunkSize - 1;
        const int outerPointX = innerPointX;
        for (int innerPointY = 0; innerPointY < chunkSize; innerPointY++)
        {
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
            if (innerPointY != (chunkSize - 1))
            { // only single triangle in corner
                pushTriangle(ra, rb, rc, wa, wb, wc);
            }
        }
    }
    else if (rightLod > lod)
    {
        const int innerPointX = chunkSize - 1;
        const int outerPointX = innerPointX;
        for (int innerPointY = 0; innerPointY < chunkSize; innerPointY++)
        {
            const int outerPointY = innerPointY / 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            const int outerPointDy = outerPointY + 1;

            if (innerPointY % 2 == 0)
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
            else
            {
                if (innerPointY != (chunkSize - 1))
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

                    pushTriangle(ra, rb, rc, wa, wb, wc);
                }
            }
        }
    }
    else if (rightLod < lod)
    {
        const int innerPointX = chunkSize - 1;
        const int outerPointX = innerPointX;
        for (int innerPointY = 0; innerPointY < chunkSize; innerPointY++)
        {
            int outerPointY = innerPointY * 2;

            const int innerPointDx = innerPointX + 1;
            const int innerPointDy = innerPointY + 1;

            int outerPointDy = outerPointY + 1;

            if (innerPointY != (chunkSize - 1))
            {
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
            }
            else
            {
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

template <typename G>
void offsetGeometry(G &geometry, const vm::ivec2 &worldPosition, float height = -(float)WORLD_BASE_HEIGHT)
{
    for (size_t i = 0; i < geometry.positions.size(); i++)
    {
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
    TerrainGeometry &geometry)
{
    const int worldSize = chunkSize * lod;
    const int worldSizeM1 = worldSize - lod;
    const int chunkSizeM1 = chunkSize - 1;
    const int rowSize = chunkSize + 2;
    createPlaneGeometry<
        Heightfield,
        TerrainGeometry,
        WindingDirection::CCW,
        ComputeNormals::YES>(
        worldSizeM1,
        worldSizeM1,
        chunkSizeM1,
        chunkSizeM1,
        rowSize,
        heightfields,
        geometry);
}
void generateHeightfieldSeamsMesh(
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    const std::vector<Heightfield> &heightfields,
    TerrainGeometry &geometry)
{
    const int rowSize = chunkSize + 2;
    createPlaneSeamsGeometry<
        Heightfield,
        TerrainGeometry,
        WindingDirection::CCW,
        ComputeNormals::YES>(
        lod,
        lodArray,
        chunkSize,
        rowSize,
        heightfields,
        geometry);
}

//

void generateWaterfieldCenterMesh(
    int lod,
    int chunkSize,
    const std::vector<Waterfield> &waterfields,
    WaterGeometry &geometry)
{
    const int worldSize = chunkSize * lod;
    const int worldSizeM1 = worldSize - lod;
    const int chunkSizeM1 = chunkSize - 1;
    const int rowSize = chunkSize + 2;
    createPlaneGeometry<
        Waterfield,
        WaterGeometry,
        WindingDirection::CCW,
        ComputeNormals::NO>(
        worldSizeM1,
        worldSizeM1,
        chunkSizeM1,
        chunkSizeM1,
        rowSize,
        waterfields,
        geometry);
}
void generateWaterfieldSeamsMesh(
    int lod,
    const std::array<int, 2> &lodArray,
    int chunkSize,
    const std::vector<Waterfield> &waterfields,
    WaterGeometry &geometry)
{
    const int rowSize = chunkSize + 2;
    createPlaneSeamsGeometry<
        Waterfield,
        WaterGeometry,
        WindingDirection::CCW,
        ComputeNormals::NO>(
        lod,
        lodArray,
        chunkSize,
        rowSize,
        waterfields,
        geometry);
}

//

template <typename G>
void createBoxGeometry(float width, float height, float depth, int widthSegments, int heightSegments, int depthSegments, G &geometry)
{
    int numberOfVertices = 0;

    auto buildPlane = [&](std::function<vm::vec3(float z, float y, float x)> &vectorFn, float udir, float vdir, float width, float height, float depth, int gridX, int gridY) -> void
    {
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

        for (int iy = 0; iy < gridY1; iy++)
        {

            const int y = iy * segmentHeight - heightHalf;

            for (int ix = 0; ix < gridX1; ix++)
            {

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

                vector = vectorFn(0, 0, depth > 0 ? 1 : -1);
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
                    (1.f - (float)iy) * height});

                // counters

                vertexCounter += 1;
            }
        }

        // indices

        // 1. you need three indices to draw a single face
        // 2. a single segment consists of two faces
        // 3. so we need to generate six (2*3) indices per segment

        for (int iy = 0; iy < gridY; iy++)
        {

            for (int ix = 0; ix < gridX; ix++)
            {

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

    std::function<vm::vec3(float z, float y, float x)> pushVertices1 = [&](float z, float y, float x) -> vm::vec3
    {
        return vm::vec3{
            x,
            y,
            z};
    };
    buildPlane(pushVertices1, -1, -1, depth, height, width, depthSegments, heightSegments); // px
    std::function<vm::vec3(float z, float y, float x)> pushVertices2 = [&](float z, float y, float x) -> vm::vec3
    {
        return vm::vec3{
            x,
            y,
            z};
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
    std::function<vm::vec3(float z, float y, float x)> pushVertices5 = [&](float x, float y, float z) -> vm::vec3
    {
        return vm::vec3{
            x,
            y,
            z};
    };
    buildPlane(pushVertices5, 1, -1, width, height, depth, widthSegments, heightSegments); // pz
    std::function<vm::vec3(float z, float y, float x)> pushVertices6 = [&](float x, float y, float z) -> vm::vec3
    {
        return vm::vec3{
            x,
            y,
            z};
    };
    buildPlane(pushVertices6, -1, -1, width, height, -depth, widthSegments, heightSegments); // nz
}

template <typename G>
G &mergeGeometries(G &dst, std::vector<G> &geometries)
{
    size_t numPositions = 0;
    for (size_t i = 0; i < geometries.size(); i++)
    {
        numPositions += geometries[i].positions.size();
    }

    size_t numNormals = 0;
    for (size_t i = 0; i < geometries.size(); i++)
    {
        numNormals += geometries[i].normals.size();
    }

    size_t numIndices = 0;
    for (size_t i = 0; i < geometries.size(); i++)
    {
        numIndices += geometries[i].indices.size();
    }

    size_t numUvs = 0;
    for (size_t i = 0; i < geometries.size(); i++)
    {
        numUvs += geometries[i].uvs.size();
    }

    size_t numPositions2D = 0;
    for (size_t i = 0; i < geometries.size(); i++)
    {
        numPositions2D += geometries[i].positions2D.size();
    }

    dst.positions.reserve(numPositions);
    dst.normals.reserve(numNormals);
    dst.indices.reserve(numIndices);

    for (size_t i = 0; i < geometries.size(); i++)
    {
        G &g = geometries[i];
        for (size_t j = 0; j < g.positions.size(); j++)
        {
            dst.positions.push_back(g.positions[j]);
        }
        for (size_t j = 0; j < g.normals.size(); j++)
        {
            dst.normals.push_back(g.normals[j]);
        }
        size_t positionOffset = dst.positions.size() / 3;
        for (size_t j = 0; j < g.indices.size(); j++)
        {
            dst.indices.push_back(positionOffset + g.indices[j]);
        }
        for (size_t j = 0; j < g.uvs.size(); j++)
        {
            dst.uvs.push_back(g.uvs[j]);
        }
        for (size_t j = 0; j < g.positions2D.size(); j++)
        {
            dst.positions2D.push_back(g.positions2D[j]);
        }
    }

    return dst;
}
void setPositions2D(BarrierGeometry &geometry, const vm::ivec2 position2D)
{
    geometry.positions2D.reserve(geometry.positions.size());
    for (size_t i = 0; i < geometry.positions.size(); i++)
    {
        geometry.positions2D.push_back(position2D);
    }
}

//

void generateBarrierGeometry(
    const vm::ivec2 &worldPosition,
    int chunkSize,
    OctreeContext &octreeContext,
    BarrierGeometry &geometry)
{
    std::vector<OctreeNodePtr> seedLeafNodes = octreeContext.getLeafNodes();

    // geometry
    std::vector<BarrierGeometry> geometries;
    for (size_t i = 0; i < seedLeafNodes.size(); i++)
    {
        OctreeNodePtr node = seedLeafNodes[i];

        BarrierGeometry g;

        const vm::ivec2 &nodePosition = node->min;
        const int &nodeLod = node->lod;

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
            g);
        vm::ivec2 worldOffset{
            width / 2 + worldPosition.x,
            depth / 2 + worldPosition.y};
        offsetGeometry(
            g,
            worldOffset,
            height / 2.f + barrierMinHeight);

        setPositions2D(g, node->min / chunkSize);

        geometries.push_back(std::move(g));
    }
    mergeGeometries(geometry, geometries);

    // leaf nodes
    geometry.leafNodes = seedLeafNodes;

    // leaf node index
    {
        vm::ivec2 leafNodesMin{// in chunks space
                               std::numeric_limits<int>::max(),
                               std::numeric_limits<int>::max()};
        vm::ivec2 leafNodesMax{// in chunks space
                               std::numeric_limits<int>::min(),
                               std::numeric_limits<int>::min()};
        for (size_t i = 0; i < seedLeafNodes.size(); i++)
        {
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

        for (size_t i = 0; i < seedLeafNodes.size(); i++)
        {
            OctreeNodePtr node = seedLeafNodes[i];
            vm::ivec2 nodeChunkPosition = node->min; // in chunks space
            const int &nodeLod = node->lod;

            int minX = nodeChunkPosition.x;
            int minZ = nodeChunkPosition.y;
            int maxX = nodeChunkPosition.x + nodeLod;
            int maxZ = nodeChunkPosition.y + nodeLod;

            for (int z = minZ; z < maxZ; z++)
            {
                for (int x = minX; x < maxX; x++)
                {
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
    TerrainGeometry &geometry)
{
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
    WaterGeometry &geometry)
{
    generateWaterfieldCenterMesh(lod, chunkSize, waterfields, geometry);
    generateWaterfieldSeamsMesh(lod, lodArray, chunkSize, waterfields, geometry);
    offsetGeometry(geometry, worldPosition);
}

//

class HeightfieldSampler
{
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
        const std::vector<Heightfield> &heightfields) : worldPositionXZ(worldPositionXZ),
                                                        lod(lod),
                                                        chunkSize(chunkSize),
                                                        chunkSizeP2(chunkSize + 2),
                                                        heightfields(heightfields)
    {
    }
    float getHeight(float x, float z)
    {
        const vm::vec2 location = getLocalPosition(x, z);

        float result = bilinear<HeightfieldSampler, float>(location, chunkSize, *this);
        return result;
    }
    Heightfield getHeightfield(float x, float z)
    {
        const vm::vec2 location = getLocalPosition(x, z);

        float rx = std::floor(location.x);
        float ry = std::floor(location.y);

        int ix = (int)rx;
        int iy = (int)ry;

        return getHeightfieldByLocalPosition(ix, iy);
    }
    float get(int x, int z)
    {
        const Heightfield heightfield = getHeightfieldByLocalPosition(x, z);
        const float &height = heightfield.height;
        return height;
    }

private:
    vm::vec2 getLocalPosition(float x, float z)
    {
        vm::vec2 location{
            x - worldPositionXZ.x,
            z - worldPositionXZ.y};
        location /= (float)lod;
        return location;
    };
    Heightfield getHeightfieldByLocalPosition(int x, int z)
    {
        const int dx = x + 1;
        const int dz = z + 1;
        const int index = dx + dz * chunkSizeP2;
        const Heightfield &heightfield = heightfields[index];
        return heightfield;
    }
};

//

void pushSplatInstances(const float &ax, const float &az, const float &rot, SplatInstanceGeometry &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler)
{
    auto iterPair = geometry.instances.emplace(std::make_pair(instanceId, SplatInstance{}));
    auto iter = iterPair.first;
    const bool &inserted = iterPair.second;
    SplatInstance &instance = iter->second;
    if (inserted)
    {
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

void pushSubSplatInstances(const float &ax, const float &az, const float &rot, SplatInstanceGeometry &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, const int &numSubObjectsPerObject, std::mt19937 &rng, std::uniform_real_distribution<float> &dis)
{
    for (int i = 0; i < numSubObjectsPerObject; i++)
    {
        const float offsetX = dis(rng) * 2.f - 1.f;
        const float offsetZ = dis(rng) * 2.f - 1.f;

        const float signX = signbit(offsetX) ? -1.f : 1.f;
        const float signZ = signbit(offsetZ) ? -1.f : 1.f;

        const float newX = (BUSH_AROUND_TREE_BASE_OFFSET * signX) + ax + offsetX * BUSH_AROUND_TREE_OFFSET_RANGE;
        const float newZ = (BUSH_AROUND_TREE_BASE_OFFSET * signZ) + az + offsetZ * BUSH_AROUND_TREE_OFFSET_RANGE;

        auto iterPair = geometry.instances.emplace(std::make_pair(instanceId, SplatInstance{}));
        auto iter = iterPair.first;
        const bool &inserted = iterPair.second;
        SplatInstance &instance = iter->second;
        if (inserted)
        {
            instance.instanceId = instanceId;
        }

        const float height = heightfieldSampler.getHeight(newX, newZ) - (float)WORLD_BASE_HEIGHT;

        instance.ps.push_back(newX);
        instance.ps.push_back(height);
        instance.ps.push_back(newZ);

        Quat q = Quat().setFromAxisAngle(Vec{0, 1, 0}, rot);
        instance.qs.push_back(q.x);
        instance.qs.push_back(q.y);
        instance.qs.push_back(q.z);
        instance.qs.push_back(q.w);
    }
}
void generateVegetationInstances(
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int chunkSize,
    const int numVegetationInstances,
    const std::vector<Heightfield> &heightfields,
    Noises &noises,
    VegetationGeometry &treeGeometry,
    VegetationGeometry &bushGeometry)
{
    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y};

    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields);

    for (int dz = 0; dz < lod; dz++)
    {
        for (int dx = 0; dx < lod; dx++)
        {
            int chunkMinX = baseMinX + dx * chunkSize;
            int chunkMinZ = baseMinZ + dz * chunkSize;

            float chunkSeed = noises.uberNoise.hashNoise(chunkMinX, chunkMinZ);
            unsigned int seedInt = *(unsigned int *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < MAX_NUM_VEGGIES_PER_CHUNK; i++)
            {
                const float chunkOffsetX = dis(rng) * (float)chunkSize;
                const float chunkOffsetZ = dis(rng) * (float)chunkSize;
                const float rot = dis(rng) * 2.0f * M_PI;
                const int instanceId = (int)std::round(dis(rng) * (float)(numVegetationInstances - 1));

                const float ax = (float)chunkMinX + chunkOffsetX;
                const float az = (float)chunkMinZ + chunkOffsetZ;

                const Heightfield &heightfield = heightfieldSampler.getHeightfield(ax, az);

                const float slope = heightfield.getSlope();

                if (slope < 0.1f)
                {
                    bool isVisible = noises.uberNoise.treeVisibility(ax, az);
                    if (isVisible)
                    {
                        pushSplatInstances(ax, az, rot, treeGeometry, instanceId, heightfieldSampler);
                        pushSubSplatInstances(ax, az, rot, bushGeometry, instanceId, heightfieldSampler, NUM_BUSHES_AROUND_TREE, rng, dis);
                    }
                }
            }
        }
    }
}

void generateRocksInstances(
    const vm::ivec2 &worldPositionXZ,
    const int lod,
    const int chunkSize,
    const int numRockInstances,
    const std::vector<Heightfield> &heightfields,
    Noises &noises,
    RockGeometry &rockGeometry,
    RockGeometry &stoneGeometry)
{
    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y};
    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields);

    for (int dz = 0; dz < lod; dz++)
    {
        for (int dx = 0; dx < lod; dx++)
        {
            const int chunkMinX = baseMinX + dx * chunkSize;
            const int chunkMinZ = baseMinZ + dz * chunkSize;

            const float chunkSeed = noises.uberNoise.hashNoise(chunkMinX, chunkMinZ);
            unsigned int seedInt = *(unsigned int *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < MAX_NUM_VEGGIES_PER_CHUNK; i++)
            {
                const float chunkOffsetX = dis(rng) * (float)chunkSize;
                const float chunkOffsetZ = dis(rng) * (float)chunkSize;
                const float rot = dis(rng) * 2.0f * M_PI;
                const int instanceId = (int)std::round(dis(rng) * (float)(numRockInstances - 1));

                const float ax = (float)chunkMinX + chunkOffsetX;
                const float az = (float)chunkMinZ + chunkOffsetZ;

                const Heightfield &heightfield = heightfieldSampler.getHeightfield(ax, az);
                const vm::vec3 &normal = heightfield.normal;

                const float slope = heightfield.getSlope();

                if (slope < 0.1f)
                {
                    bool isVisible = noises.uberNoise.rockVisibility(ax, az);
                    if (isVisible)
                    {
                        pushSplatInstances(ax, az, rot, rockGeometry, instanceId, heightfieldSampler);
                        pushSubSplatInstances(ax, az, rot, stoneGeometry, instanceId, heightfieldSampler, NUM_STONES_AROUND_ROCK, rng, dis);
                    }
                }
            }
        }
    }
}
//

template <typename T, typename G>
G &pushMaterialAwareSplatInstances(const float &ax, const float &az, const float &rot, T &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield)
{
    auto iterPair = geometry.instances.emplace(std::make_pair(instanceId, G{}));
    auto iter = iterPair.first;
    const bool &inserted = iterPair.second;
    G &instance = iter->second;
    if (inserted)
    {
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

    const MaterialsArray &heightfieldMaterials = heightfield.materials;
    const MaterialsWeightsArray &heightfieldMaterialsWeights = heightfield.materialsWeights;
    const vm::vec4 materials = vm::vec4{(float)heightfieldMaterials[0], (float)heightfieldMaterials[1], (float)heightfieldMaterials[2], (float)heightfieldMaterials[3]};
    const vm::vec4 materialsWeights = vm::vec4{heightfieldMaterialsWeights[0], heightfieldMaterialsWeights[1], heightfieldMaterialsWeights[2], heightfieldMaterialsWeights[3]};
    instance.materials.push_back(materials);
    instance.materialsWeights.push_back(materialsWeights);

    return instance;
}
void pushGrassInstances(const float &ax, const float &az, const float &rot, GrassGeometry &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield, Noises &noises)
{
    GrassSplatInstance &instance = pushMaterialAwareSplatInstances<GrassGeometry, GrassSplatInstance>(ax, az, rot, geometry, instanceId, heightfieldSampler, heightfield);

    const float heightScale = GRASS_MODEL_BASE_HEIGHT + ((noises.uberNoise.simplexNoise(ax * 7.f, ax * 7.f) * 2.f) - 1.f) * GRASS_HEIGHT_VARIATION_RANGE;

    const float colorVariationNoise = vm::clamp(GRASS_COLOR_VARIATION_BASE + (noises.uberNoise.simplexNoise(ax * 10.f, az * 10.f) * 2.f - 1.f) * GRASS_COLOR_VARIATION_RANGE, 0.0f, 1.f + GRASS_COLOR_VARIATION_RANGE);
    const float randomBladeFactor = noises.uberNoise.hashNoise(ax, az) * 2.f - 1.f;

    const vm::vec3 grassColorMultiplier = (vm::vec3{colorVariationNoise, colorVariationNoise / 1.1f, colorVariationNoise / 1.2f} +
                                           vm::vec3{randomBladeFactor / 10.f, randomBladeFactor / 12.f, randomBladeFactor / 14.f});

    const vm::vec4 grassProps = vm::vec4{grassColorMultiplier.x, grassColorMultiplier.y, grassColorMultiplier.z, heightScale};

    instance.grassProps.push_back(grassProps);
}

// template <typename G, typename I>
// void generateInstances(const vm::ivec2 &worldPositionXZ,
//                        const int lod,
//                        const int chunkSize,
//                        const int numInstances,
//                        const std::vector<Heightfield> &heightfields,
//                        Noises &noises,
//                        G &geometry)
// {
//     int baseMinX = worldPositionXZ.x;
//     int baseMinZ = worldPositionXZ.y;

//     vm::vec2 worldPositionXZf{
//         (float)worldPositionXZ.x,
//         (float)worldPositionXZ.y};

//     HeightfieldSampler heightfieldSampler(
//         worldPositionXZf,
//         lod,
//         chunkSize,
//         heightfields);

//     for (int dz = 0; dz < lod; dz++)
//     {
//         for (int dx = 0; dx < lod; dx++)
//         {
//             const int chunkMinX = baseMinX + dx * chunkSize;
//             const int chunkMinZ = baseMinZ + dz * chunkSize;

//             const float chunkSeed = noises.uberNoise.hashNoise(chunkMinX, chunkMinZ);
//             uint32_t seedInt = *(uint32_t *)&chunkSeed;
//             std::mt19937 rng(seedInt);
//             std::uniform_real_distribution<float> dis(0.f, 1.f);

//             for (int i = 0; i < MAX_NUM_GRASSES_PER_CHUNK; i++)
//             {
//                 const float chunkOffsetX = dis(rng) * (float)chunkSize;
//                 const float chunkOffsetZ = dis(rng) * (float)chunkSize;
//                 const float rot = dis(rng) * 2.0f * M_PI;
//                 const int instanceId = (int)std::round(dis(rng) * (float)(numGrassInstances - 1));

//                 const float ax = (float)chunkMinX + chunkOffsetX;
//                 const float az = (float)chunkMinZ + chunkOffsetZ;

//                 const Heightfield &heightfield = heightfieldSampler.getHeightfield(ax, az);

//                 const float slope = heightfield.getSlope();

//                 if (slope < 0.1f)
//                 {
//                     const bool isVisible = noises.uberNoise.instanceVisibility<I>(ax, az);
//                     if (isVisible)
//                     {
//                         const float crushedGrassNoise = 1.f - noises.uberNoise.rockVisibility(ax, az);
//                         if (crushedGrassNoise > CRUSHED_GRASS_THRESHOLD)
//                         {
//                             pushGrassInstances(ax, az, rot, grassGeometry, instanceId, heightfieldSampler, heightfield, noises, crushedGrassNoise);
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

void generateGrassInstances(const vm::ivec2 &worldPositionXZ,
                            const int lod,
                            const int chunkSize,
                            const int numGrassInstances,
                            const std::vector<Heightfield> &heightfields,
                            Noises &noises,
                            GrassGeometry &grassGeometry)
{
    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y};

    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields);

    for (int dz = 0; dz < lod; dz++)
    {
        for (int dx = 0; dx < lod; dx++)
        {
            const int chunkMinX = baseMinX + dx * chunkSize;
            const int chunkMinZ = baseMinZ + dz * chunkSize;

            const float chunkSeed = noises.uberNoise.hashNoise(chunkMinX, chunkMinZ);
            const uint32_t seedInt = *(uint32_t *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < MAX_NUM_GRASSES_PER_CHUNK; i++)
            {
                const float chunkOffsetX = dis(rng) * (float)chunkSize;
                const float chunkOffsetZ = dis(rng) * (float)chunkSize;
                const float rot = dis(rng) * 2.0f * M_PI;
                const int instanceId = (int)std::round(dis(rng) * (float)(numGrassInstances - 1));

                const float ax = (float)chunkMinX + chunkOffsetX;
                const float az = (float)chunkMinZ + chunkOffsetZ;

                const Heightfield &heightfield = heightfieldSampler.getHeightfield(ax, az);

                const float slope = heightfield.getSlope();

                if (slope < 0.1f)
                {
                    bool isVisible = noises.uberNoise.grassVisibility(ax, az);
                    if (isVisible)
                    {
                        pushGrassInstances(ax, az, rot, grassGeometry, instanceId, heightfieldSampler, heightfield, noises);
                    }
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
    PoiGeometry &poiGeometry)
{
    constexpr int maxNumPoisPerChunk = 16;
    constexpr float poiRate = 0.3;
    // const float poiRate = maxPoiRate / (float)(lod * lod);
    const float poiThrowRate = 1.f / (float)lod;
    // const float poiRate = maxPoiRate;

    int baseMinX = worldPositionXZ.x;
    int baseMinZ = worldPositionXZ.y;

    vm::vec2 worldPositionXZf{
        (float)worldPositionXZ.x,
        (float)worldPositionXZ.y};
    HeightfieldSampler heightfieldSampler(
        worldPositionXZf,
        lod,
        chunkSize,
        heightfields);

    for (int dz = 0; dz < lod; dz++)
    {
        for (int dx = 0; dx < lod; dx++)
        {
            int chunkMinX = baseMinX + dx * chunkSize;
            int chunkMinZ = baseMinZ + dz * chunkSize;

            float chunkSeed = noises.poiSeedNoise.in2D(chunkMinX, chunkMinZ);
            unsigned int seedInt = *(unsigned int *)&chunkSeed;
            std::mt19937 rng(seedInt);
            std::uniform_real_distribution<float> dis(0.f, 1.f);

            for (int i = 0; i < maxNumPoisPerChunk; i++)
            {
                float noiseValue = dis(rng);
                float chunkOffsetX = dis(rng) * (float)chunkSize;
                float chunkOffsetZ = dis(rng) * (float)chunkSize;
                int instanceId = (int)std::round(dis(rng) * (float)(numPoiInstances - 1));

                if (noiseValue < poiRate)
                {
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

        float chunkSeed = generator.noises.numSplitsNoise.in2D(baseNode.x, baseNode.y);
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
    const int chunkSize = generator.getChunkSize();

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
        heightfields = generator.getHeightfields(worldPosition.x, worldPosition.y, lod, lodArray);
        calculateSurfaceNormals(heightfields, lod, lodArray, chunkSize);
        generator.applyMaterials(worldPosition.x, worldPosition.y, lod, lodArray, heightfields);
    }

    // terrain
    if (generateFlags & GF_TERRAIN)
    {
        TerrainGeometry terrainGeometry;

        generateTerrainGeometry(
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

        generateWaterGeometry(
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

        generateVegetationInstances(
            worldPosition,
            lod,
            chunkSize,
            numVegetationInstances,
            heightfields,
            generator.noises,
            treeGeometry,
            bushGeometry);

        result->treeInstancesBuffer = treeGeometry.getBuffer();
        result->bushInstancesBuffer = bushGeometry.getBuffer();
    }
    else
    {
        result->treeInstancesBuffer = nullptr;
        result->bushInstancesBuffer = nullptr;
    }

    // rocks
    if (generateFlags & GF_ROCK)
    {
        RockGeometry rockGeometry;
        RockGeometry stoneGeometry;

        generateRocksInstances(
            worldPosition,
            lod,
            chunkSize,
            numRockInstances,
            heightfields,
            generator.noises,
            rockGeometry,
            stoneGeometry);

        result->rockInstancesBuffer = rockGeometry.getBuffer();
        result->stoneInstancesBuffer = stoneGeometry.getBuffer();
    }
    else
    {
        result->rockInstancesBuffer = nullptr;
        result->stoneInstancesBuffer = nullptr;
    }

    // grass
    if (generateFlags & GF_GRASS)
    {
        GrassGeometry grassGeometry;

        if (lod <= 2)
        {
            generateGrassInstances(
                worldPosition,
                lod,
                chunkSize,
                numGrassInstances,
                heightfields,
                generator.noises,
                grassGeometry);
        }

        result->grassInstancesBuffer = grassGeometry.getBuffer();
    }
    else
    {
        result->grassInstancesBuffer = nullptr;
    }

    // poi
    if (generateFlags & GF_POI)
    {
        PoiGeometry poiGeometry;

        generatePoiInstances(
            worldPosition,
            lod,
            chunkSize,
            numPoiInstances,
            heightfields,
            generator.noises,
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
    const int chunkSize = generator.getChunkSize();

    OctreeContext octreeContext = getChunkSeedOctree(
        worldPosition,
        minLod,
        maxLod,
        chunkSize);

    BarrierGeometry barrierGeometry;

    generateBarrierGeometry(
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
    const int chunkSize = generator.getChunkSize();

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