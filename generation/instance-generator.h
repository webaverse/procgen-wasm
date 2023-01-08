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

typedef std::function<void(const float &, const float &, const vm::vec3 &, const float &, const vm::vec3 &, const int &, HeightfieldSampler &, const Heightfield &, std::mt19937 &, std::uniform_real_distribution<float> &)> PushInstancesFunction;
#define INSTANCE_PUSH_FN_PARAMS const float &ax, const float &az, const vm::vec3 &rot, const float &scale, const vm::vec3 &color, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield, std::mt19937 &rng, std::uniform_real_distribution<float> &dis

class InstanceGenerator
{
public:
    template <typename G, typename GG>
    G &getVariationGeometry(const float &ax, const float &az, GG &geometryGroup, Noises &noises)
    {
        // TODO: move this to glsl
        const float simpm5 = noises.uberNoise.simplexNoise(ax * 5.f, az * 5.f);
        const float hashm5d4 = (noises.uberNoise.hashNoise(ax * 5.f, az * 5.f) * 2.f - 1.f) / 4.f;
        const float randomTreePicker = vm::clamp(simpm5 + hashm5d4, 0.f, 1.f);

        const int geometryIndex = std::round((geometryGroup.geometries.size() - 1) * randomTreePicker);
        return geometryGroup.geometries[geometryIndex];
    }
    void pushSplatInstances(const float &ax, const float &az, const vm::vec3 &rot, const float &scale, const vm::vec3 &color, SplatInstanceGeometry &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, std::mt19937 &rng, std::uniform_real_distribution<float> &dis)
    {
        auto iterPair = geometry.instances.emplace(std::make_pair(instanceId, SplatInstance{}));
        auto iter = iterPair.first;
        const bool &inserted = iterPair.second;
        SplatInstance &instance = iter->second;
        if (inserted)
        {
            instance.instanceId = instanceId;
        }

        const float height = heightfieldSampler.getWorldHeight(ax, az);

        instance.set(vm::vec3{ax, height, az}, rot, scale, color);
    }

    void pushSubSplatInstances(const float &ax, const float &az, const vm::vec3 &rot, const float &scale, const vm::vec3 &color, SplatInstanceGeometry &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, const int &numObjects, std::mt19937 &rng, std::uniform_real_distribution<float> &dis, const float &baseOffset, const float &offsetRange)
    {
        for (int i = 0; i < numObjects; i++)
        {
            const float offsetX = dis(rng) * 2.f - 1.f;
            const float offsetZ = dis(rng) * 2.f - 1.f;

            const float signX = signbit(offsetX) ? -1.f : 1.f;
            const float signZ = signbit(offsetZ) ? -1.f : 1.f;

            const float newX = (baseOffset * signX) + ax + offsetX * offsetRange;
            const float newZ = (baseOffset * signZ) + az + offsetZ * offsetRange;

            auto iterPair = geometry.instances.emplace(std::make_pair(instanceId, SplatInstance{}));
            auto iter = iterPair.first;
            const bool &inserted = iterPair.second;
            SplatInstance &instance = iter->second;
            if (inserted)
            {
                instance.instanceId = instanceId;
            }

            const float height = heightfieldSampler.getHeight(newX, newZ) - (float)WORLD_BASE_HEIGHT;

            instance.set(vm::vec3{newX, height, newZ}, rot, scale, color);
        }
    }

    template <typename T, typename G>
    G &pushMaterialAwareSplatInstances(const float &ax, const float &az, const vm::vec3 &rot, const float &scale, const vm::vec3 &color, T &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield)
    {
        auto iterPair = geometry.instances.emplace(std::make_pair(instanceId, G{}));
        auto iter = iterPair.first;
        const bool &inserted = iterPair.second;
        G &instance = iter->second;
        if (inserted)
        {
            instance.instanceId = instanceId;
        }

        const float height = heightfieldSampler.getWorldHeight(ax, az);

        instance.set(vm::vec3{ax, height, az}, rot, scale, color);

        const MaterialsArray &heightfieldMaterials = heightfield.materials;
        const MaterialsWeightsArray &heightfieldMaterialsWeights = heightfield.materialsWeights;
        const vm::vec4 materials = vm::vec4{(float)heightfieldMaterials[0], (float)heightfieldMaterials[1], (float)heightfieldMaterials[2], (float)heightfieldMaterials[3]};
        const vm::vec4 materialsWeights = vm::vec4{heightfieldMaterialsWeights[0], heightfieldMaterialsWeights[1], heightfieldMaterialsWeights[2], heightfieldMaterialsWeights[3]};
        instance.materials.push_back(materials);
        instance.materialsWeights.push_back(materialsWeights);

        return instance;
    }

    virtual bool validateHeightfield(const uint8_t &id, const Heightfield &heightfield)
    {
        // validate heightfield
        return true;
    }

    template <uint8_t I, typename G>
    void generateInstances(
        const vm::ivec2 &worldPositionXZ,
        const int lod,
        const int chunkSize,
        const int numInstances,
        const int maxNumInstancesPerChunk,
        const std::vector<Heightfield> &heightfields,
        Noises &noises,
        const PushInstancesFunction &pushInstancesFunction)
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
                uint32_t seedInt = *(uint32_t *)&chunkSeed;
                std::mt19937 rng(seedInt);
                std::uniform_real_distribution<float> dis(0.f, 1.f);

                for (int i = 0; i < maxNumInstancesPerChunk; i++)
                {
                    const float chunkOffsetX = dis(rng) * (float)chunkSize;
                    const float chunkOffsetZ = dis(rng) * (float)chunkSize;

                    const float ax = (float)chunkMinX + chunkOffsetX;
                    const float az = (float)chunkMinZ + chunkOffsetZ;

                    const Heightfield &heightfield = heightfieldSampler.getHeightfield(ax, az);

                    // TODO :: this shouldn't be a universal rule for instances, this should only apply to instances which don't exist in water
                    const bool hasWater = heightfield.hasWater();
                    if (!hasWater)
                    {
                        if (validateHeightfield(I, heightfield))
                        {
                            if (noises.uberNoise.instanceVisibility<I>(ax, az, heightfield))
                            {
                                // ! Remove numInstances ?
                                const int instanceId = (int)std::round(dis(rng) * (float)(numInstances - 1));
                                const float scale = noises.uberNoise.scaleNoise<I>(ax, az, heightfield);
                                const vm::vec3 rot = noises.uberNoise.rotationNoise<I>(ax, az, heightfield);
                                const vm::vec3 color = noises.uberNoise.colorNoise<I>(ax, az, heightfield);

                                pushInstancesFunction(ax, az, rot, scale, color, instanceId, heightfieldSampler, heightfield, rng, dis);
                            }
                        }
                    }
                }
            }
        }
    }
};

class VegetationGenerator : public InstanceGenerator
{
public:
    bool validateHeightfield(const uint8_t &id, const Heightfield &heightfield) override
    {
        bool valid = true;

        // validate biome
        if (valid)
        {
            const uint8_t &dominantBiome = heightfield.getDominantBiome();

            switch (id)
            {
            case (uint8_t)INSTANCE::GRASS:
                // valid = dominantBiome != (uint8_t)BIOME::DESERT;
                break;
            default:
                break;
            }
        }

        // validate slope
        if (valid)
        {
            const float slope = heightfield.getSlope();
            valid = slope < SLOPE_CUTOFF;
        }

        return valid;
    }

    void pushGrassInstances(const float &ax, const float &az, const vm::vec3 &rot, const float &scale, const vm::vec3 &color, GrassGeometry &geometry, const int &instanceId, HeightfieldSampler &heightfieldSampler, const Heightfield &heightfield, Noises &noises)
    {
        GrassSplatInstance &instance = pushMaterialAwareSplatInstances<GrassGeometry, GrassSplatInstance>(ax, az, rot, scale, color, geometry, instanceId, heightfieldSampler, heightfield);

        const float simplexm10 = noises.uberNoise.simplexNoise(ax * 10.f, az * 10.f) * 2.f - 1.f;
        const float heightScale = GRASS_MODEL_BASE_HEIGHT + simplexm10 * GRASS_HEIGHT_VARIATION_RANGE;

        const vm::vec4 grassProps = vm::vec4{heightfield.field.hash, heightfield.field.hash, heightfield.field.hash, heightScale};
        instance.grassProps.push_back(grassProps);
    }
};

class MineralGenerator : public InstanceGenerator
{
public:
    bool validateHeightfield(const uint8_t &id, const Heightfield &heightfield) override
    {
        bool valid = true;

        // validate slope
        if (valid)
        {
            const float slope = heightfield.getSlope();
            valid = slope < SLOPE_CUTOFF;
        }

        return valid;
    }
};

class PoiGenerator
{
public:
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
                        const float height = heightfieldSampler.getWorldHeight(ax, az);

                        poiGeometry.ps.push_back(ax);
                        poiGeometry.ps.push_back(height);
                        poiGeometry.ps.push_back(az);

                        poiGeometry.instances.push_back(instanceId);
                    }
                }
            }
        }
    }
};

#endif // __INSTANCED-GENERATOR_H__