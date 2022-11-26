#ifndef __HEIGHTFIELD_GENERATOR_H__
#define __HEIGHTFIELD_GENERATOR_H__

#include "noises.h"
#include "heightfield.h"
#include "biomes.h"

const int numBiomes = (int)BIOME::NUM_BIOMES;
const int numLiquids = (int)LIQUID::NUM_LIQUIDS;
const int numMaterials = (int)MATERIAL::NUM_MATERIALS;

enum GenerateFlags
{
    GF_NONE = 0,
    GF_TERRAIN = 1 << 0,
    GF_WATER = 1 << 1,
    GF_VEGETATION = 1 << 2,
    GF_ROCK = 1 << 3,
    GF_GRASS = 1 << 4,
    GF_POI = 1 << 5,
    GF_HEIGHTFIELD = 1 << 6
};

class GenerationSettings
{
public:
    int seed;
    int chunkSize;

    GenerationSettings(int seed, int chunkSize) : seed(seed), chunkSize(chunkSize){};
};

class HeightfieldGenerator
{
public:
    GenerationSettings settings;
    Noises noises;

    HeightfieldGenerator(int seed, int chunkSize) : noises(seed), settings(seed, chunkSize){};

    int getChunkSize();

    std::vector<Heightfield> getHeightfields(
        int x,
        int z,
        int lod,
        const std::array<int, 2> &lodArray);

    BiomeNoiseField getBiomeNoiseField(float bx, float by);
    LiquidNoiseField getLiquidNoiseField(float bx, float by);

    uint8_t getBiome(float bx, float bz);
    uint8_t getLiquid(float bx, float bz, uint8_t biome);

    void getHeightFieldCenter(int bx, int bz, int lod, std::vector<Heightfield> &heightfields);
    void getHeightFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, int rowSize, std::vector<Heightfield> &heightfieldSeams);
    Heightfield getHeightField(float bx, float bz);
    float getHeight(float bx, float bz);

    float getComputedBiomeHeight(unsigned char b, const vm::vec2 &worldPosition);
    float getComputedWaterHeight(const float &height, const float &realHeight, uint8_t liquid);
    float getComputedTerrainHeight(const float &height, const vm::vec2 &worldPosition);

    // materials
    void setHeightfieldMaterial(Heightfield &localHeightfield, const vm::vec2 &position);
    void applyCenterMaterials(const int &bx, const int &bz, const int &lod, std::vector<Heightfield> &heightfields);
    void applySeamMaterials(const int &bx, const int &bz, const int &lod, const std::array<int, 2> &lodArray, const int &rowSize, std::vector<Heightfield> &heightfieldSeams);
    void applyMaterials(const int &x, const int &z, const int &lod, const std::array<int, 2> &lodArray, std::vector<Heightfield> &heightfields);
    void getComputedMaterials(Heightfield &localHeightfield, std::vector<MaterialWeightAccumulator> &materialsCounts, float &totalMaterialFactors, const vm::vec2 &worldPosition);

    float getTemperature(const vm::vec2 &worldPosition, const int &lod);
    float getHumidity(const vm::vec2 &worldPosition, const int &lod);
};

#endif // __HEIGHTFIELD_GENERATOR_H__