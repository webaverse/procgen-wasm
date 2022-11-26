#ifndef __POLYGONIZER_H__
#define __POLYGONIZER_H__

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <ctime>

#include "../libs/vectorMath.h"
#include "../libs/vector.h"
#include "../task/octree.h"
#include "../task/tracker.h"
#include "mesh.h"

class Polygonizer
{
public:
    Polygonizer();
    ~Polygonizer();

    void calculateSurfaceNormals(std::vector<Heightfield> &heightfields,
                                              const int &lod,
                                              const std::array<int, 2> &lodArray,
                                              const int &chunkSize);

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