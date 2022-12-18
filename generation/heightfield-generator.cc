#include "heightfield-generator.h"

BiomeNoiseField HeightfieldGenerator::getBiomeNoiseField(float bx, float bz)
{
    float tNoise = noises.uberNoise.temperatureNoise(bx, bz);
    float hNoise = noises.uberNoise.humidityNoise(bx, bz);

    return BiomeNoiseField{
        tNoise,
        hNoise};
}

LiquidNoiseField HeightfieldGenerator::getLiquidNoiseField(float bx, float bz)
{
    float oNoise = noises.uberNoise.oceanNoise(bx, bz);
    float rNoise = noises.uberNoise.riverNoise(bx, bz, oNoise);

    return LiquidNoiseField{
        oNoise,
        rNoise};
}

uint8_t HeightfieldGenerator::getBiome(float bx, float bz)
{
    uint8_t biome = 0xFF;

    // const auto &noise = getBiomeNoiseField(bx, bz);
    // const float temperature = noise.temperature;
    // const float coldness = 1.f - temperature;
    // const float humidity = noise.humidity;
    // const float biomeFactor = coldness * humidity;

    // const bool isCold = getNoiseVisibility(biomeFactor, COLD_WARM_BORDER, BIOME_BORDER_MAX);
    // if (isCold)
    // {
    //     biome = (uint8_t)BIOME::TAIGA;
    // }
    // const bool isWarm = getNoiseVisibility(biomeFactor, WARM_HOT_BORDER, COLD_WARM_BORDER);
    // if (isWarm)
    // {
    biome = (uint8_t)BIOME::FOREST;
    // }
    // const bool isHot = getNoiseVisibility(biomeFactor, BIOME_BORDER_MIN, WARM_HOT_BORDER);
    // if (isHot)
    // {
    //     biome = (uint8_t)BIOME::DESERT;
    // }

    return biome;
}

uint8_t HeightfieldGenerator::getLiquid(float bx, float bz, uint8_t biome)
{
    uint8_t liquid = (uint8_t)LIQUID::NULL_LIQUID;

    const auto &noise = getLiquidNoiseField(bx, bz);
    const float ocean = noise.ocean;
    const float river = noise.river;

    const bool isOcean = getNoiseVisibility(ocean, OCEAN_THRESHOLD, 1.f);
    const bool isRiver = getNoiseVisibility(river, RIVER_THRESHOLD, 1.f);

    if (isOcean)
    {
        liquid = (uint8_t)LIQUID::OCEAN;
    }

    if (isRiver)
    {
        if (biome == (uint8_t)BIOME::DESERT)
        {
            return liquid;
        }

        if (liquid == (uint8_t)LIQUID::OCEAN)
        {
            const bool isDeepOcean = getNoiseVisibility(ocean, OCEAN_THRESHOLD + DEEP_OCEAN_AND_WATERFALL_DIFF, 1.f);
            if(isDeepOcean) {
                liquid = (uint8_t)LIQUID::OCEAN;
            }
            else
            {
                liquid = (uint8_t)LIQUID::WATERFALL;
            }
        }
        else
        {
            liquid = (uint8_t)LIQUID::RIVER;
        }
    }

    return liquid;
}

//

void HeightfieldGenerator::getHeightFieldCenter(int bx, int bz, int lod, std::vector<Heightfield> &heightfields)
{
    const int chunkSizeP2 = settings.chunkSize + 2;
    for (int dz = 0; dz < chunkSizeP2; dz++)
    {
        for (int dx = 0; dx < chunkSizeP2; dx++)
        {
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
void HeightfieldGenerator::getHeightFieldSeams(int bx, int bz, int lod, const std::array<int, 2> &lodArray, int rowSize, std::vector<Heightfield> &heightfields)
{
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = settings.chunkSize * lod / bottomLod;
    const int gridWidthP3 = gridWidth + 3;

    const int gridHeight = settings.chunkSize * lod / rightLod;
    const int gridHeightP3 = gridHeight + 3;

    const int heightfieldsCenterDataOffset = rowSize * rowSize;

    // bottom
    int index = heightfieldsCenterDataOffset;
    {
        for (int dz = 0; dz < 3; dz++)
        {
            for (int dx = 0; dx < gridWidthP3; dx++)
            {
                const int x = dx - 1;
                const int z = dz - 1;

                Heightfield &localHeightfieldSeam = heightfields[index];

                const int ax = bx + x * bottomLod;
                const int az = bz + settings.chunkSize * lod + z * bottomLod;
                localHeightfieldSeam = getHeightField(ax, az);

                index++;
            }
        }
    }
    // right
    {
        for (int dx = 0; dx < 3; dx++)
        {
            for (int dz = 0; dz < gridHeightP3; dz++)
            {
                const int x = dx - 1;
                const int z = dz - 1;

                Heightfield &localHeightfieldSeam = heightfields[index];

                const int ax = bx + settings.chunkSize * lod + x * rightLod;
                const int az = bz + z * rightLod;
                localHeightfieldSeam = getHeightField(ax, az);

                index++;
            }
        }
    }
}

std::vector<Heightfield> HeightfieldGenerator::getHeightfields(
    int x,
    int z,
    int lod,
    const std::array<int, 2> &lodArray)
{
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = settings.chunkSize * lod / bottomLod;
    const int gridWidthP1 = gridWidth + 1;
    const int gridWidthP3 = gridWidth + 3;
    const int gridWidthP3T3 = gridWidthP3 * 3;

    const int gridHeight = settings.chunkSize * lod / rightLod;
    const int gridHeightP1 = gridHeight + 1;
    const int gridHeightP3 = gridHeight + 3;
    const int gridHeightP3T3 = gridHeightP3 * 3;

    const int chunkSizeP2 = settings.chunkSize + 2;

    std::vector<Heightfield> heightfields(
        (chunkSizeP2 * chunkSizeP2) +    // center
        (gridWidthP3T3 + gridHeightP3T3) // seams
    );

    getHeightFieldCenter(x, z, lod, heightfields);
    getHeightFieldSeams(x, z, lod, lodArray, chunkSizeP2, heightfields);

    return heightfields;
}

std::vector<uint8_t> sortWeightedTypes(const std::vector<float> &weights)
{
    std::vector<uint8_t> seenTypes;
    for (size_t i = 0; i < weights.size(); i++)
    {
        const uint8_t type = (uint8_t)i;
        const float &typeWeight = weights[type];
        if (typeWeight > 0.f)
        {
            seenTypes.push_back(type);
        }
    }

    // sort by increasing occurrence count of the type
    std::sort(
        seenTypes.begin(),
        seenTypes.end(),
        [&](uint8_t b1, uint8_t b2) -> bool
        {
            return weights[b1] > weights[b2];
        });

    return seenTypes;
}

Heightfield HeightfieldGenerator::getHeightField(float bx, float bz)
{
    Heightfield localHeightfield;

    vm::vec2 fWorldPosition{bx, bz};

    // height
    const float halfChunkSizeF = (float)settings.chunkSize / 2.f;
    const float maxDistance = std::sqrt(halfChunkSizeF);

    // water
    constexpr int waterRange = 4;
    const float maxWaterDistance = (float)std::sqrt((float)waterRange * (float)waterRange);
    constexpr float baseWaterFactor = 0.25;

    // accumulators
    float totalHeightFactors = 0;
    float sumLiquidFactors = 0;

    // weights
    std::vector<float> biomeWeights(numBiomes, 0.f);
    std::vector<float> liquidWeights(numLiquids, 0.f);

    // loop
    for (float dz = -halfChunkSizeF; dz <= halfChunkSizeF; dz++)
    {
        for (float dx = -halfChunkSizeF; dx <= halfChunkSizeF; dx++)
        {
            float distance = std::sqrt(dx * dx + dz * dz);

            float ax = bx + dx;
            float az = bz + dz;

            if (distance < maxDistance)
            {
                uint8_t b = getBiome(ax, az);
                uint8_t l = getLiquid(ax, az, b);

                float heightFactor = 1.f - (distance / maxDistance);

                biomeWeights[b] += heightFactor;
                liquidWeights[l] += heightFactor;
                totalHeightFactors += heightFactor;

                if (l != (uint8_t)LIQUID::NULL_LIQUID)
                {
                    sumLiquidFactors += heightFactor;
                }
            }
        }
    }

    // liquid factor
    sumLiquidFactors /= totalHeightFactors;
    localHeightfield.field.liquidFactor = sumLiquidFactors;

    // postprocess height
    {
        const std::vector<uint8_t> &seenBiomes = sortWeightedTypes(biomeWeights);
        const std::vector<uint8_t> &seenLiquids = sortWeightedTypes(liquidWeights);

        // letting the biome weight fit in an uint8_t
        const float WEIGHT_FITTER = 255.f;

        for (size_t i = 0; i < 4; i++)
        {
            if (i < seenBiomes.size())
            {
                const uint8_t &biome = seenBiomes[i];
                localHeightfield.biomes[i] = biome;
                localHeightfield.biomesWeights[i] = biomeWeights[biome] / totalHeightFactors;
            }
            else
            {
                localHeightfield.biomes[i] = (uint8_t)BIOME::NULL_BIOME;
                localHeightfield.biomesWeights[i] = 0;
            }

            if (i < seenLiquids.size())
            {
                const uint8_t &liquid = seenLiquids[i];
                localHeightfield.liquids[i] = liquid;
                localHeightfield.liquidsWeights[i] = liquidWeights[liquid] / totalHeightFactors;
            }
            else
            {
                localHeightfield.liquids[i] = (uint8_t)LIQUID::NULL_LIQUID;
                localHeightfield.liquidsWeights[i] = 0;
            }
        }

        float elevationSum = 0.f;
        float liquidElevationSum = 0.f;

        vm::vec2 fWorldPosition{bx, bz};

        for (size_t i = 0; i < biomeWeights.size(); i++)
        {
            const uint8_t biome = (uint8_t)i;
            const float &biomeWeight = biomeWeights[biome] / totalHeightFactors;

            if (biomeWeight > 0.f)
            {
                const float computedBiomeHeight = getComputedBiomeHeight(biome, fWorldPosition);
                const float computedTerrainHeight = getComputedTerrainHeight(biome, computedBiomeHeight, fWorldPosition);

                elevationSum += biomeWeight * computedTerrainHeight;

                for (size_t j = 0; j < liquidWeights.size(); j++)
                {
                    const uint8_t liquid = (uint8_t)j;
                    const float &liquidWeight = liquidWeights[liquid] / totalHeightFactors;

                    if (liquidWeight > 0.f)
                    {
                        const float computedLiquidHeight = getComputedWaterHeight(computedBiomeHeight, computedTerrainHeight, liquid);
                        liquidElevationSum += biomeWeight * liquidWeight * computedLiquidHeight;
                    }
                }
            }
        }

        localHeightfield.field.height = elevationSum;
        localHeightfield.field.liquidHeight = liquidElevationSum;
    }

    // caching noises
    {
        localHeightfield.field.hash = noises.uberNoise.hashNoise(fWorldPosition.x, fWorldPosition.y);
        localHeightfield.field.wetness = noises.uberNoise.wetnessNoise(fWorldPosition.x, fWorldPosition.y);
        localHeightfield.field.tree = noises.uberNoise.treeVisibility(fWorldPosition.x, fWorldPosition.y, localHeightfield.field.wetness);
        localHeightfield.field.grass = noises.uberNoise.grassMaterialNoise(fWorldPosition.x, fWorldPosition.y, localHeightfield.field.wetness);
        localHeightfield.field.rock = noises.uberNoise.rockVisibility(fWorldPosition.x, fWorldPosition.y);
        localHeightfield.field.stone = noises.uberNoise.stoneVisibility(fWorldPosition.x, fWorldPosition.y);
        localHeightfield.field.flower = noises.uberNoise.flowerVisibility(fWorldPosition.x, fWorldPosition.y, localHeightfield.field.grass);
    }

    return localHeightfield;
}
float HeightfieldGenerator::getHeight(float bx, float bz)
{
    const float halfChunkSizeF = (float)settings.chunkSize / 2.f;
    const float maxDistance = std::sqrt(halfChunkSizeF + 1.f);

    std::unordered_map<unsigned char, float> biomeCounts(numBiomes);
    float totalHeightFactors = 0;
    for (float dz = -halfChunkSizeF; dz <= halfChunkSizeF; dz++)
    {
        for (float dx = -halfChunkSizeF; dx <= halfChunkSizeF; dx++)
        {
            float distance = std::sqrt(dx * dx + dz * dz);
            float factor = std::max(1.f - (distance / maxDistance), 0.f);

            if (factor > 0)
            {
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

// float randomFromPoint(int x, int y, int z)
// {
//     uint64_t hash = hashOctreeMin(vm::ivec3{x, y, z});
//     uint32_t hash32 = (uint32_t)hash ^ (uint32_t)(hash >> 32);
//     float f = (float)hash32 / (float)UINT32_MAX;
//     return f;
// }

// biomes
float HeightfieldGenerator::getComputedBiomeHeight(const uint8_t &biome, const vm::vec2 &worldPosition)
{
    const float &ax = worldPosition.x;
    const float &az = worldPosition.y;

    switch (biome)
    {
    case (int)BIOME::NULL_BIOME:
        return 0.f;
    case (int)BIOME::DESERT:
        return noises.uberNoise.desertNoise(ax, az);
    case (int)BIOME::FOREST:
        return noises.uberNoise.mountainNoise(ax, az);
    case (int)BIOME::TAIGA:
        return noises.uberNoise.iceMountainNoise(ax, az);
    default:
        return noises.uberNoise.mountainNoise(ax, az);
    }
}

float HeightfieldGenerator::getComputedTerrainHeight(const uint8_t &biome, const float &height, const vm::vec2 &worldPosition)
{
    const float &ax = worldPosition.x;
    const float &az = worldPosition.y;

    float terrainHeight = height;

    switch (biome)
    {
    case (uint8_t)BIOME::DESERT:
        // we don't have water in desert !
        break;
    default:
        terrainHeight -= noises.uberNoise.waterDepthNoise(ax, az);
        break;
    }

    return vm::clamp(terrainHeight, (float)MIN_WORLD_HEIGHT, (float)MAX_WORLD_HEIGHT);
}

float HeightfieldGenerator::getComputedWaterHeight(const float &biomeHeight, const float &terrainHeight, uint8_t liquid)
{

    const float belowBiomeHeight = biomeHeight - WATER_HEIGHT_DIFFERENCE;
    const float belowTerrainHeight = terrainHeight - WATER_HEIGHT_DIFFERENCE;

    float liquidHeight = 0.f;

    // * handle the height of water differently based on different biomes
    switch (liquid)
    {
    case (int)LIQUID::NULL_LIQUID:
        liquidHeight = belowTerrainHeight;
        break;
    case (int)LIQUID::OCEAN:
        liquidHeight = (float)WATER_BASE_HEIGHT;
        break;
    case (int)LIQUID::RIVER:
        liquidHeight = belowBiomeHeight;
        break;
    case (int)LIQUID::LAVA:
        liquidHeight = belowBiomeHeight;
        break;
    case (int)LIQUID::WATERFALL:
        liquidHeight = terrainHeight + WATER_OFFSET;
        break;
    default:
        liquidHeight = belowBiomeHeight;
        break;
    }

    return vm::clamp(liquidHeight, (float)WATER_BASE_HEIGHT, (float)MAX_WORLD_HEIGHT);
}

// materials
void HeightfieldGenerator::setHeightfieldMaterial(Heightfield &localHeightfield, const vm::vec2 &position)
{
    std::vector<MaterialWeightAccumulator> materialWeightAccumulators(numMaterials);
    float totalMaterialFactors = 0;

    getComputedMaterials(localHeightfield, materialWeightAccumulators, totalMaterialFactors, position);

    std::vector<uint8_t> seenMaterials;
    for (size_t i = 0; i < materialWeightAccumulators.size(); i++)
    {
        const uint8_t material = (uint8_t)i;
        const bool &isSeen = materialWeightAccumulators[material].getSeen();
        if (isSeen)
        {
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
void HeightfieldGenerator::applyCenterMaterials(const int &bx, const int &bz, const int &lod, std::vector<Heightfield> &heightfields)
{
    const int chunkSizeP2 = settings.chunkSize + 2;
    for (int dz = 0; dz < chunkSizeP2; dz++)
    {
        for (int dx = 0; dx < chunkSizeP2; dx++)
        {
            const int index = dx + dz * chunkSizeP2;
            Heightfield &localHeightfield = heightfields[index];

            const int x = dx - 1;
            const int z = dz - 1;
            const int ax = bx + x * lod;
            const int az = bz + z * lod;

            const vm::vec2 fWorldPosition{(float)ax, (float)az};

            setHeightfieldMaterial(localHeightfield, fWorldPosition);
        }
    }
}
void HeightfieldGenerator::applySeamMaterials(const int &bx, const int &bz, const int &lod, const std::array<int, 2> &lodArray, const int &rowSize, std::vector<Heightfield> &heightfields)
{
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = settings.chunkSize * lod / bottomLod;
    const int gridWidthP3 = gridWidth + 3;

    const int gridHeight = settings.chunkSize * lod / rightLod;
    const int gridHeightP3 = gridHeight + 3;

    const int heightfieldsCenterDataOffset = rowSize * rowSize;

    // bottom
    int index = heightfieldsCenterDataOffset;
    {
        for (int dz = 0; dz < 3; dz++)
        {
            for (int dx = 0; dx < gridWidthP3; dx++)
            {
                const int x = dx - 1;
                const int z = dz - 1;

                Heightfield &localHeightfieldSeam = heightfields[index];

                const int ax = bx + x * bottomLod;
                const int az = bz + settings.chunkSize * lod + z * bottomLod;

                const vm::vec2 &fWorldPosition{(float)ax, (float)az};
                setHeightfieldMaterial(localHeightfieldSeam, fWorldPosition);

                index++;
            }
        }
    }
    // right
    {
        for (int dx = 0; dx < 3; dx++)
        {
            for (int dz = 0; dz < gridHeightP3; dz++)
            {
                const int x = dx - 1;
                const int z = dz - 1;

                Heightfield &localHeightfieldSeam = heightfields[index];

                const int ax = bx + settings.chunkSize * lod + x * rightLod;
                const int az = bz + z * rightLod;

                const vm::vec2 &fWorldPosition{(float)ax, (float)az};

                setHeightfieldMaterial(localHeightfieldSeam, fWorldPosition);

                index++;
            }
        }
    }
}
void HeightfieldGenerator::applyMaterials(const int &x, const int &z, const int &lod, const std::array<int, 2> &lodArray, std::vector<Heightfield> &heightfields)
{
    const int &bottomLod = lodArray[0];
    const int &rightLod = lodArray[1];

    const int gridWidth = settings.chunkSize * lod / bottomLod;
    const int gridWidthP1 = gridWidth + 1;
    const int gridWidthP3 = gridWidth + 3;
    const int gridWidthP3T3 = gridWidthP3 * 3;

    const int gridHeight = settings.chunkSize * lod / rightLod;
    const int gridHeightP1 = gridHeight + 1;
    const int gridHeightP3 = gridHeight + 3;
    const int gridHeightP3T3 = gridHeightP3 * 3;

    const int chunkSizeP2 = settings.chunkSize + 2;

    applyCenterMaterials(x, z, lod, heightfields);
    applySeamMaterials(x, z, lod, lodArray, chunkSizeP2, heightfields);
}

void HeightfieldGenerator::getComputedMaterials(Heightfield &localHeightfield, std::vector<MaterialWeightAccumulator> &materialWeightAccumulators, float &totalMaterialFactors, const vm::vec2 &worldPosition)
{
    const std::array<uint8_t, 4> &biomes = localHeightfield.biomes;
    const std::array<float, 4> &biomesWeights = localHeightfield.biomesWeights;

    const int GRASS = (int)MATERIAL::GRASS;
    const int DIRT = (int)MATERIAL::DIRT;
    const int ROCK = (int)MATERIAL::ROCK;
    const int STONE = (int)MATERIAL::STONE;

    for (int i = 0; i < 1; i++)
    {
        const uint8_t b = biomes[i];
        const float bw = (float)biomesWeights[i];

        switch (b)
        {
        // TODO : Define a different set of material rules for each biome, for now we're using these rules as default
        default:
            const float grassMaterial = localHeightfield.field.grass;
            const float height = localHeightfield.field.height;

            const float snowFactor = std::max(height - SNOW_START_HEIGHT, 0.f) / SNOW_RANGE;

            const float stiffness = noises.uberNoise.coldnessNoise(worldPosition.x, worldPosition.y) * snowFactor;

            // amplifying the slope of the terrain to blend between ground and cliffs
            const float SLOPE_AMPLIFIER = 2.5f;
            const float mountainAndGroundBlend = vm::clamp(localHeightfield.getSlope() * SLOPE_AMPLIFIER, 0.f, 1.f);

            const float grassWeight = (grassMaterial) * bw * (1.f - mountainAndGroundBlend);
            const float dirtWeight = (1.f - grassMaterial) * bw * (1.f - mountainAndGroundBlend);

            const float stoneWeight = (stiffness) * mountainAndGroundBlend * bw;
            const float rockWeight = (1.f - stiffness) * mountainAndGroundBlend * bw;

            MaterialWeightAccumulator &grassWeightAcc = materialWeightAccumulators[GRASS];
            MaterialWeightAccumulator &dirtWeightAcc = materialWeightAccumulators[DIRT];
            MaterialWeightAccumulator &rockWeightAcc = materialWeightAccumulators[ROCK];
            MaterialWeightAccumulator &stoneWeightAcc = materialWeightAccumulators[STONE];

            grassWeightAcc.addWeight(grassWeight);
            dirtWeightAcc.addWeight(dirtWeight);
            rockWeightAcc.addWeight(rockWeight);
            stoneWeightAcc.addWeight(stoneWeight);

            totalMaterialFactors += grassWeight;
            totalMaterialFactors += dirtWeight;
            totalMaterialFactors += rockWeight;
            totalMaterialFactors += stoneWeight;

            break;
        }
    }
}

int HeightfieldGenerator::getChunkSize()
{
    return settings.chunkSize;
}
