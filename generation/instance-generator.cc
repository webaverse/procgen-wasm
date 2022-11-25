#include "instance-generator.h"

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

void InstanceGenerator::generateVegetationInstances(
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

void InstanceGenerator::generateRocksInstances(
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
// void InstanceGenerator::generateInstances(const vm::ivec2 &worldPositionXZ,
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

void InstanceGenerator::generateGrassInstances(const vm::ivec2 &worldPositionXZ,
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

void InstanceGenerator::generatePoiInstances(
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
