#include "instance-generator.h"

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
