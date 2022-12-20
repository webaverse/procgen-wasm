#ifndef __POLYGONIZER_H__
#define __POLYGONIZER_H__

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <ctime>
#include <type_traits>

#include "../libs/vectorMath.h"
#include "../libs/vector.h"
#include "../task/octree.h"
#include "../task/tracker.h"
#include "mesh.h"

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

class Polygonizer
{
public:
    Polygonizer();
    ~Polygonizer();

    template <typename T>
    void calculateSurfaceNormals(std::vector<T> &heightfields,
                                 const int &lod,
                                 const std::array<int, 2> &lodArray,
                                 const int &chunkSize)
    {

        const int rowSize = chunkSize + 2;

        // TODO : This is not the best way to calculate normals (it lacks precision)
        calculateCenterNormals<T>(heightfields, chunkSize, rowSize);
        calculateSeamNormals<T>(heightfields, lod, lodArray, chunkSize, rowSize);
    };

    void generateBarrierGeometry(
        const vm::ivec2 &worldPosition,
        int chunkSize,
        OctreeContext &octreeContext,
        BarrierGeometry &geometry);

    void generateTerrainGeometry(
        const vm::ivec2 &worldPosition,
        int lod,
        const std::array<int, 2> &lodArray,
        int chunkSize,
        const std::vector<Heightfield> &heightfields,
        TerrainGeometry &geometry);

    void generateWaterGeometry(
        const vm::ivec2 &worldPosition,
        int lod,
        const std::array<int, 2> &lodArray,
        int chunkSize,
        const std::vector<Waterfield> &waterfields,
        WaterGeometry &geometry);
};

#endif // __POLYGONIZER_H__