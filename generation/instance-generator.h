#ifndef __INSTANCED_GENERATOR_H__
#define __INSTANCED_GENERATOR_H__

#include "../polygonization/mesh.h"
#include "noises.h"
#include "../constants.h"
#include "../libs/vectorMath.h"
#include "../libs/vector.h"
#include "../libs/MurmurHash3.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <ctime>
#include <string.h>
#include <memory>

class InstanceGenerator
{
public:
    void generateVegetationInstances(
        const vm::ivec2 &worldPositionXZ,
        const int lod,
        const int chunkSize,
        const int numVegetationInstances,
        const std::vector<Heightfield> &heightfields,
        Noises &noises,
        VegetationGeometry &treeGeometry,
        VegetationGeometry &bushGeometry);

    void generateRocksInstances(
        const vm::ivec2 &worldPositionXZ,
        const int lod,
        const int chunkSize,
        const int numRockInstances,
        const std::vector<Heightfield> &heightfields,
        Noises &noises,
        RockGeometry &rockGeometry,
        RockGeometry &stoneGeometry);

    void generateGrassInstances(const vm::ivec2 &worldPositionXZ,
                                const int lod,
                                const int chunkSize,
                                const int numGrassInstances,
                                const std::vector<Heightfield> &heightfields,
                                Noises &noises,
                                GrassGeometry &grassGeometry);

    void generatePoiInstances(
        const vm::ivec2 &worldPositionXZ,
        const int lod,
        const int chunkSize,
        const int numPoiInstances,
        const std::vector<Heightfield> &heightfields,
        Noises &noises,
        PoiGeometry &poiGeometry);
};

#endif // __INSTANCED-GENERATOR_H__