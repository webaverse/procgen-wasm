#include "glsl.h"
#include "constants.h"
#include "noises.h"

using namespace GLSL;

// ! GLSL INPUTS

const float MAX_TERRAIN_HEIGHT = float(MAX_WORLD_HEIGHT);
const float MIN_TERRAIN_HEIGHT = float(MIN_WORLD_HEIGHT);

// mountains heights
const float TERRAIN_SHORT_HILLS_HEIGHT = SHORT_HILLS_HEIGHT;
const float TERRAIN_SHORT_MOUNTAIN_HEIGHT = SHORT_MOUNTAIN_HEIGHT;
const float TERRAIN_FLAT_SHORT_HILLS_HEIGHT = FLAT_SHORT_HILLS_HEIGHT;
const float TERRAIN_TALL_MOUNTAIN_HEIGHT = TALL_MOUNTAIN_HEIGHT;
const float TERRAIN_ICE_MOUNTAIN_HEIGHT = ICE_MOUNTAIN_HEIGHT;
const float TERRAIN_SAND_MOUNTAIN_HEIGHT = SAND_MOUNTAIN_HEIGHT;
// biomes borders
const float TERRAIN_BIOME_BORDER_MIN = BIOME_BORDER_MIN;
const float TERRAIN_BIOME_BORDER_MAX = BIOME_BORDER_MAX;
const float TERRAIN_COLD_WARM_BORDER = COLD_WARM_BORDER;
const float TERRAIN_WARM_HOT_BORDER = WARM_HOT_BORDER;
const float TERRAIN_BIOME_BORDERS_FLATTENER_DELTA = BIOME_BORDERS_FLATTENER_DELTA;
// water
const float WATER_SURROUNDING_HEIGHT = float(WORLD_BASE_HEIGHT / 4);
const float TERRAIN_OCEAN_DEPTH = float(OCEAN_DEPTH);
const float TERRAIN_MAX_WATER_DEPTH = TERRAIN_OCEAN_DEPTH;
const float TERRAIN_RIVER_DEPTH = float(RIVER_DEPTH);
// ? making sure the terrain surface is above water before water depth subtraction
const float TERRAIN_BASE_HEIGHT = float(WATER_SURROUNDING_HEIGHT + WATER_BASE_HEIGHT);
const float TERRAIN_OCEAN_THRESHOLD = OCEAN_THRESHOLD;
const float TERRAIN_RIVER_THRESHOLD = RIVER_BASE;
const float TERRAIN_WATER_THRESHOLD = WATER_THRESHOLD;
const float TERRAIN_RIVER_EDGES_DELTA = RIVER_EDGES_DELTA;
// grass
const float TERRAIN_GRASS_THRESHOLD = GRASS_THRESHOLD;
// flowers
const float TERRAIN_FLOWER_THRESHOLD = FLOWER_THRESHOLD;
// trees
const float TERRAIN_TREE_THRESHOLD = VEGGIE_THRESHOLD;
// rocks
const float TERRAIN_ROCK_THRESHOLD = ROCK_THRESHOLD;
const float TERRAIN_STONE_THRESHOLD = STONE_THRESHOLD;
const float TERRAIN_ROCK_NOISE_SHARPNESS = ROCK_NOISE_SHARPNESS;
// flattener
const float TERRAIN_FLATTENER_DEPTH = float(WATER_BASE_HEIGHT * 2);

// ----------------------------------

// * START GLSL ---------------------

// constants
const float NOISE_SCALE = 512.f;

// util functions

// ? smooth and soft clamp implementation from : https://www.shadertoy.com/view/Ws3Xzr

// reaches the limit inside the domain bounds, but has slope 2 at the midpoint
float smoothClamp(float x, float a, float b)
{
    float t = clamp(x, a, b);
    return t != x ? t : b + (a - b) / (1. + exp((b - a) * (2. * x - a - b) / ((x - a) * (b - x))));
}

// approaches the limit more slowly with slope <= 1 everywhere
float softClamp(float x, float a, float b)
{
    float mid = (a + b) * 0.5;
    return mid + smoothClamp((x - mid) * 0.5, a - mid, b - mid);
}

float smoothEdge(float value, float edge)
{
    float clampedValue = smoothClamp(value, edge, 1.f);
    // making the value fit in the 0-1 range
    return ((clampedValue - edge) * (1.f / (1.f - edge)));
}

float getSmoothBorder(float biomeFactor, float border, float delta)
{
    float multiplier = 1.f / delta;
    return clamp(1.f - (smoothClamp(biomeFactor, border, border + delta) - smoothClamp(biomeFactor, border - delta, border)) * multiplier, 0.f, 1.f);
}

// noise functions

// ? Simplex noise implementation from : https://gist.github.com/patriciogonzalezvivo/670c22f3966e662d2f83

// (sqrt(3)-1)/2, (3-sqrt(3))/6
const vec4 C = vec4(0.211324865405187f, 0.366025403784439f, -0.577350269189626f, 0.024390243902439f);

// ! The simplex noise function below is slower than the one being used (Based on profiling)

// #define ffSNOISE 1.79284291400159f
// #define fdSNOISE 0.85373472095314f

// vec3 permute(vec3 x) { return mod(((x * 34.f) + 1.f) * x, 289.f); }

// float snoise(vec2 v)
// {
//     vec2 i = floor(v + dot(v, vec2(C.y, C.y)));
//     vec2 x0 = v - i + dot(i, vec2(C.x, C.x));
//     vec2 i1 = (x0.x > x0.y) ? vec2(1.f, 0.f) : vec2(0.f, 1.f);
//     vec4 x12 = vec4(x0.x, x0.y, x0.x, x0.y) + vec4(C.x, C.x, C.z, C.z) - vec4(i1.x, i1.y, 0.f, 0.f);
//     i = mod(i, 289.f);
//     vec3 p = permute(permute(vec3(0.f, i1.y, 1.f) + i.y) + i.x + vec3(0.f, i1.x, 1.f));
//     vec3 m = max(vec3(0.5f, 0.5f, 0.5f) -
//                      vec3(dot(x0, x0), dot(vec2(x12.x, x12.y), vec2(x12.x, x12.y)), dot(vec2(x12.z, x12.w), vec2(x12.z, x12.w))),
//                  0.f);
//     m = m * m;
//     m = m * m;
//     vec3 x = fract(p * vec3(C.w, C.w, C.w)) * 2.f - 1.f;
//     p = abs(x) - 0.5f;
//     vec3 ox = floor(x + 0.5f);
//     vec3 a0 = x - ox;
//     m *= vec3(ffSNOISE, ffSNOISE, ffSNOISE) - (a0 * a0 + p * p) * fdSNOISE;
//     i1 = vec2(a0.y, a0.z) * vec2(x12.x, x12.z) + vec2(p.y, p.z) * vec2(x12.y, x12.w);
//     p = vec3(a0.x * x0.x + p.x * x0.y, i1.x, i1.y);
//     return 135.f * dot(m, p);
// }

float hash2D_1(vec2 p)
{
    float h = dot(p, vec2(127.1, 311.7));
    return fract(sin(h) * 43758.5453123);
}

vec2 hash2D_2(vec2 p)
{
    p = vec2(dot(p, vec2(127.1, 311.7)), dot(p, vec2(269.5, 183.3)));
    return vec2(-1.f, -1.f) + fract(sin(p) * 43758.5453123) * 2.f;
}
float snoise(vec2 p)
{
    vec2 i = floor(p + (p.x + p.y) * C.y);
    vec2 a = p - i + (i.x + i.y) * C.x;
    float m = step(a.y, a.x);
    vec2 o = vec2(m, 1.f - m);
    vec2 b = a - o + C.x;
    vec2 c = a - 1.f + 2.f * C.x;
    vec3 h = max(vec3(-dot(a, a), -dot(b, b), -dot(c, c)) + 0.5f, 0.f);
    vec3 n = h * h * h * h * vec3(dot(a, hash2D_2(i)), dot(b, hash2D_2(i + o)), dot(c, hash2D_2(i + 1.f)));

    // ? Original Implementation
    // return dot(n, vec3(70.f, 70.f, 70.f));

    // ? Modified Implementation for better results
    return clamp((dot(n, vec3(140.f, 140.f, 140.f)) + 0.45f) / 2.f, 0.f, 1.f);
}
// float lerp(float a, float b, float t)
// {
// 	return a + t * (b - a);
// }
// float snoise(vec2 p)
// {
// 	vec2 i = floor(p);
// 	vec2 f = fract(p);

//     //grid points
//     vec2 p0 = vec2(0.0, 0.0);
//     vec2 p1 = vec2(1.0, 0.0);
//     vec2 p2 = vec2(0.0, 1.0);
//     vec2 p3 = vec2(1.0, 1.0);

//     //distance vectors to each grid point
//     vec2 s0 = f - p0;
//     vec2 s1 = f - p1;
//     vec2 s2 = f - p2;
//     vec2 s3 = f - p3;

//     //random gradient vectors on each grid point
//     vec2 g0 = hash2D_2(i + p0);
//     vec2 g1 = hash2D_2(i + p1);
//     vec2 g2 = hash2D_2(i + p2);
//     vec2 g3 = hash2D_2(i + p3);

//     //gradient values
//     float q0 = dot(s0, g0);
//     float q1 = dot(s1, g1);
//     float q2 = dot(s2, g2);
//     float q3 = dot(s3, g3);

//     //interpolant weights
//     vec2 u = f * f * (vec2(-f.x, -f.y) * 2.0 + 3.0);

//     //bilinear interpolation
//     float l0 = lerp(q0, q1, u.x);
//     float l1 = lerp(q2, q3, u.x);
//     float l2 = lerp(l0, l1, u.y);

//     return clamp(l2 * 1.5f, 0.f, 1.f);
// }

float simplex(vec2 position)
{
    const float frequency = 0.74f;
    const float SCALE = 1.f / NOISE_SCALE;
    vec2 p = position * SCALE * frequency;
    return snoise(p);
}

float voronoi(vec2 x)
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    float res = 8.0;
    for (int j = -1; j <= 1; j++)
    {
        for (int i = -1; i <= 1; i++)
        {
            vec2 b = vec2(float(i), float(j));
            float hash = hash2D_1(p + (b));
            vec2 r = b - f + vec2(hash, hash);
            float d = (r.x * r.x + r.y * r.y);
            res = min(res, d);
        }
    }
    return (res);
}

float voronoiNoise(vec2 position)
{
    const float SCALE = 1.f / NOISE_SCALE;
    vec2 p = position * SCALE;
    return voronoi(p);
}

float FBM_8(vec2 position)
{
    const int octaves = 8;
    const float frequency = 0.6f;
    const float lacunarity = 3.f;
    const float persistence = 0.3f;

    const float SCALE = 1.f / NOISE_SCALE;

    vec2 p = position * SCALE;
    float noise = 0.f;

    float amplitude = 1.f;
    p *= frequency;

    for (int i = 0; i < octaves; i++)
    {
        noise += snoise(p) * amplitude;
        p *= lacunarity;
        amplitude *= persistence;
    }

    return noise;
}

float FBM_4(vec2 position)
{
    const int octaves = 4;
    const float frequency = 0.8f;
    const float lacunarity = 3.f;
    const float persistence = 0.3f;

    const float SCALE = 1.f / NOISE_SCALE;

    vec2 p = position * SCALE;
    float noise = 0.f;

    float amplitude = 1.f;
    p *= frequency;

    for (int i = 0; i < octaves; i++)
    {
        noise += snoise(p) * amplitude;
        p *= lacunarity;
        amplitude *= persistence;
    }

    return noise;
}

float FBM_2(vec2 position)
{
    const int octaves = 2;
    const float frequency = 0.74f;
    const float lacunarity = 3.f;
    const float persistence = 0.5f;

    const float SCALE = 1.f / NOISE_SCALE;

    vec2 p = position * SCALE;
    float noise = 0.f;

    float amplitude = 1.f;
    p *= frequency;

    for (int i = 0; i < octaves; i++)
    {
        noise += snoise(p) * amplitude;
        p *= lacunarity;
        amplitude *= persistence;
    }

    return noise;
}

// custom noise functions

// ? Domain Warping : https://www.shadertoy.com/view/4s23zz
float warpNoise1Layer_1(vec2 position)
{
    vec2 q = vec2(simplex(position + vec2(0.0f, 0.0f)),
                  simplex(position + vec2(7.4f, 30.2f)));

    return simplex(position + q * NOISE_SCALE * 2.f);
}

float warpNoise1Layer_2(vec2 position)
{
    vec2 q = vec2(FBM_2(position + vec2(0.0f, 0.0f)),
                  FBM_2(position + vec2(7.4f, 30.2f)));

    return FBM_2(position + q * NOISE_SCALE * 2.f);
}

float warpNoise1Layer_4(vec2 position)
{
    vec2 q = vec2(FBM_4(position + vec2(0.0f, 0.0f)),
                  FBM_4(position + vec2(7.4f, 30.2f)));

    return FBM_4(position + q * NOISE_SCALE * 2.f);
}

float warpNoise2Layer_2(vec2 position)
{
    vec2 q = vec2(FBM_2(position + vec2(0.0f, 0.0f)),
                  FBM_2(position + vec2(7.4f, 30.2f)));

    vec2 r = vec2(FBM_2(position + q * NOISE_SCALE * 2.f + vec2(1.7, 9.2)),
                  FBM_2(position + q * NOISE_SCALE * 2.f + vec2(8.3, 2.8)));

    return FBM_2(position + r * NOISE_SCALE * 2.f);
}

float warpNoise2Layer_4(vec2 position)
{
    vec2 q = vec2(FBM_4(position + vec2(0.0f, 0.0f)),
                  FBM_4(position + vec2(7.4f, 30.2f)));

    vec2 r = vec2(FBM_4(position + q * NOISE_SCALE * 2.f + vec2(1.7, 9.2)),
                  FBM_4(position + q * NOISE_SCALE * 2.f + vec2(8.3, 2.8)));

    return FBM_4(position + r * NOISE_SCALE * 2.f);
}

// terrain noises

vec3 getWaterFlow(vec2 position)
{
    return vec3(1.f, 0.f, 0.0f);
}

float getOceanNoise(vec2 position)
{
    return clamp(simplex(position / 8.f) * 2.f, 0.f, 1.f);
}

float getRiverNoise(vec2 position, float ocean)
{
    float connectOceans = (1. - ocean);
    float riverNoise = simplex(position / 2.);
    float river = 1.0 - abs(riverNoise * 8.f * connectOceans - 0.5);
    const float RIVER_EDGE = 0.6;
    float smoothRiver = smoothEdge(river, RIVER_EDGE);
    return clamp(smoothRiver, 0.f, 1.f);
}

float getSmoothRiverBorder(float biomeFactor, float border)
{
    const float delta = TERRAIN_RIVER_EDGES_DELTA;
    return getSmoothBorder(biomeFactor, border, delta);
}

float getRiverEdges(vec2 position)
{
    float ocean = getOceanNoise(position);

    float smoothRiverEdges = getSmoothRiverBorder(getRiverNoise(position, ocean), TERRAIN_RIVER_THRESHOLD);

    // we only want river edges so we need to make sure where there's an ocean there's no river edges
    float isOceanVisible = step(TERRAIN_WATER_THRESHOLD, ocean);

    float breakingPattern = clamp(simplex(position * 10.0) + 0.2f, 0.f, 1.f);
    return smoothRiverEdges * breakingPattern - isOceanVisible;
}

bool getOceanVisibility(vec2 position)
{
    const float waterEdge = TERRAIN_WATER_THRESHOLD;
    float ocean = getOceanNoise(position);
    return ocean >= waterEdge;
}

bool getRiverVisibility(vec2 position, float ocean)
{
    const float waterEdge = TERRAIN_WATER_THRESHOLD;
    float river = getRiverNoise(position, ocean);
    return river >= waterEdge;
}

bool getWaterVisibility(vec2 position)
{
    const float waterEdge = TERRAIN_WATER_THRESHOLD;
    float ocean = getOceanNoise(position);
    float river = getRiverNoise(position, ocean);
    float water = ocean + river;
    return water >= waterEdge;
}

float getWaterDepth(vec2 position)
{
    const float waterEdge = TERRAIN_OCEAN_THRESHOLD;
    float ocean = getOceanNoise(position);
    float river = getRiverNoise(position, ocean);
    float softOcean = smoothEdge(ocean, waterEdge);
    float softRiver = smoothEdge(river, waterEdge);
    float water = softOcean * TERRAIN_OCEAN_DEPTH + softRiver * TERRAIN_RIVER_DEPTH;
    return clamp(water, 0.f, TERRAIN_MAX_WATER_DEPTH);
}

float getWetness(vec2 position)
{
    return clamp(simplex(position / 2.f) + 0.2f, 0.f, 1.f);
}

float getGrassMaterial(vec2 position, float wetness)
{
    float noise = warpNoise1Layer_2(position * 10.f) * wetness;
    return clamp(noise, 0.f, 1.f);
}

float getGrassObject(vec2 position, float wetness)
{
    return getGrassMaterial(position, wetness);
}

bool getGrassVisibility(vec2 position, float wetness)
{
    float grass = getGrassObject(position, wetness);
    return grass >= TERRAIN_GRASS_THRESHOLD;
}

float getFlowerObject(vec2 position, float grass)
{
    vec2 q = vec2(grass*grass, grass*grass*grass);
    float flower = warpNoise1Layer_2(position + q * NOISE_SCALE / 5.f);
    return clamp(flower, 0.f, 1.f);
}

bool getFlowerVisibility(vec2 position, float grass)
{
    float flower = getFlowerObject(position, grass);
    return flower >= TERRAIN_FLOWER_THRESHOLD;
}

float getTreeObject(vec2 position, float wetness)
{
    return clamp(simplex(position * 20.f) * wetness, 0.f, 1.f);
}

bool getTreeVisibility(vec2 position, float wetness)
{
    float tree = getTreeObject(position, wetness);
    return tree >= TERRAIN_TREE_THRESHOLD;
}

float getRockObject(vec2 position)
{
    return clamp(simplex(position * 15.f) * TERRAIN_ROCK_NOISE_SHARPNESS - (TERRAIN_ROCK_NOISE_SHARPNESS - 1.f), 0.f, 1.f);
}
bool getRockVisibility(vec2 position)
{
    float rock = getRockObject(position);
    return rock >= TERRAIN_ROCK_THRESHOLD;
}

float getStoneObject(vec2 position)
{
    return getRiverEdges(position);
}

bool getStoneVisibility(vec2 position)
{
    float stone = getStoneObject(position);
    return stone >= TERRAIN_STONE_THRESHOLD;
}

float getColdness(vec2 position)
{
    return clamp(warpNoise1Layer_1(position * 5.f) + 0.4f, 0.f, 1.f);
}

float getHumidity(vec2 position)
{
    float noise = simplex(position / 25.f);
    return noise;
}

float getTemperature(vec2 position)
{
    return simplex(position / 65.f);
}

float getSmoothBiomeBorder(float biomeFactor, float border)
{
    const float delta = TERRAIN_BIOME_BORDERS_FLATTENER_DELTA;
    return getSmoothBorder(biomeFactor, border, delta);
}

float getBiomeFactor(vec2 position)
{
    float cold = 1.f - getTemperature(position);
    float humidity = getHumidity(position);
    return cold * humidity;
}

float getBiomeBorders(vec2 position)
{
    float biomeFactor = getBiomeFactor(position);
    float coldWarmBorder = getSmoothBiomeBorder(biomeFactor, TERRAIN_COLD_WARM_BORDER);
    float warmHotBorder = getSmoothBiomeBorder(biomeFactor, TERRAIN_WARM_HOT_BORDER);
    return coldWarmBorder + warmHotBorder;
}

float getFlattenerNoise(vec2 position)
{
    float simd5clamped = clamp(simplex(position / 5.f) - 0.2f, 0.f, 1.f);
    // making sure the transitions between biomes are smooth
    float biomeBorders = getBiomeBorders(position);
    float flattener = simd5clamped + biomeBorders;
    float flatAreas = flattener * TERRAIN_FLATTENER_DEPTH;
    return clamp(flatAreas, 0.f, TERRAIN_FLATTENER_DEPTH);
}

// terrain height noises

float terrainHeightWrapper(vec2 position, float height)
{
    float flatAreas = getFlattenerNoise(position);
    float flattenedHeight = TERRAIN_BASE_HEIGHT + height - flatAreas;
    return clamp(flattenedHeight, TERRAIN_BASE_HEIGHT, MAX_TERRAIN_HEIGHT);
}

float getSandMountainHeight(vec2 position)
{
    // calculating noises
    float warp1l1 = warpNoise1Layer_1(position);
    float fbm4d3 = FBM_4(position / 3.f);
    float subtract = (fbm4d3 + warp1l1) / 5.f;
    float vorphalf = voronoiNoise(position) + 0.5f;

    // layering noises
    float sandNoise = clamp(vorphalf - subtract, 0.f, 1.f);
    float sandMountain = sandNoise * TERRAIN_SAND_MOUNTAIN_HEIGHT;

    return terrainHeightWrapper(position, sandMountain);
}

float getIceMountainHeight(vec2 position)
{
    // calculating noises
    float fbm4d3 = FBM_4(position / 3.f);
    float fbm8d3 = FBM_8(position / 3.f);
    float fbm2d6 = FBM_2(position / 6.f);
    float vord3 = voronoiNoise(position / 3.f);

    // layering noises
    float snowMountains = clamp(((2.f - (fbm4d3 + fbm8d3 + fbm2d6) * 2.f) + vord3) / 4.f, 0.f, 1.f);
    float snowMountainsLayer = snowMountains * TERRAIN_ICE_MOUNTAIN_HEIGHT;
    return terrainHeightWrapper(position, snowMountainsLayer);
}

float getMountainHillsHeight(vec2 position)
{
    // calculating noises
    float fbm2d2 = FBM_2(position / 2.f);
    float fbm2d2clamped = clamp(fbm2d2, 0.f, 1.f);

    float fbm2d3 = FBM_2(position / 3.f);

    float fbm2d4 = FBM_2(position / 4.f);
    float fbm2d4clamped = clamp(fbm2d4, 0.f, 1.f);
    float fbm2d4clampedm2 = clamp(fbm2d4 - 0.15f, 0.f, 0.5f) * 2.f;

    float warp2d4 = warpNoise1Layer_1(position / 4.f);

    // layering noises

    float shortHills = clamp(fbm2d2 * 2.f - fbm2d2clamped / 5.f + 0.7f, 0.f, 1.f);

    float shortHillsLayer = shortHills * TERRAIN_SHORT_HILLS_HEIGHT;

    float shortMountains = (fbm2d4 + warp2d4 / 5.f) * TERRAIN_SHORT_MOUNTAIN_HEIGHT;

    float shortHillsMountainsBlender = clamp(fbm2d3 * 5.f, 0.f, 1.f);
    float shortMountainsLayer = mix(shortHillsLayer, shortMountains, shortHillsMountainsBlender);

    float flatAreaShortHills = clamp(fbm2d2, 0.f, 0.5f);
    float flatAreaLayer = flatAreaShortHills * TERRAIN_FLAT_SHORT_HILLS_HEIGHT;

    float tallMountains = (fbm2d3 - warp2d4 / 10.f);
    float tallMountainsLayer = tallMountains * TERRAIN_TALL_MOUNTAIN_HEIGHT;

    float shortAndTallMountainsBlender = fbm2d2clamped;
    float tallAreaLayer = mix(tallMountainsLayer, shortMountainsLayer, shortAndTallMountainsBlender);

    float highAndFlatAreaBlender = fbm2d4clampedm2;
    float mountainHeight = mix(tallAreaLayer, flatAreaLayer, highAndFlatAreaBlender);

    return terrainHeightWrapper(position, mountainHeight);
}

// * END GLSL ---------------------

// --------------------------------

// ! The name of the c++ functions should be different than the name of the glsl functions

float GLSL::humidityNoise(const vec2 &position)
{
    return getHumidity(position);
}

float GLSL::grassMaterialNoise(const vec2 &position, const float &wetness)
{
    return getGrassMaterial(position, wetness);
}

float GLSL::wetnessNoise(const vec2 &position)
{
    return getWetness(position);
}

bool GLSL::grassVisibility(const vec2 &position, const float &wetness)
{
    return getGrassVisibility(position, wetness);
}

bool GLSL::treeVisibility(const vec2 &position, const float &wetness)
{
    return getTreeVisibility(position, wetness);
}

float GLSL::coldnessNoise(const vec2 &position)
{
    return getColdness(position);
}

float GLSL::temperatureNoise(const vec2 &position)
{
    return getTemperature(position);
}

float GLSL::desertNoise(const vec2 &position)
{
    return getSandMountainHeight(position);
}

float GLSL::mountainNoise(const vec2 &position)
{
    return getMountainHillsHeight(position);
}

float GLSL::iceMountainNoise(const vec2 &position)
{
    return getIceMountainHeight(position);
}

float GLSL::hashNoise(const vec2 &position)
{
    return hash2D_1(position);
}

float GLSL::simplexNoise(const vec2 &position)
{
    return simplex(position);
}
float GLSL::waterDepthNoise(const vec2 &position)
{
    return getWaterDepth(position);
}
float GLSL::oceanNoise(const vec2 &position)
{
    return getOceanNoise(position);
}
float GLSL::riverNoise(const vec2 &position, const float &ocean)
{
    return getRiverNoise(position, ocean);
}

vm::vec3 GLSL::flowNoise(const vec2 &position)
{
    return getWaterFlow(position);
}

bool GLSL::flowerVisibility(const vec2 &position, const float &grass)
{
    return getFlowerVisibility(position, grass);
}
bool GLSL::waterVisibilityNoise(const vec2 &position)
{
    return getWaterVisibility(position);
}
bool GLSL::rockVisibility(const vec2 &position)
{
    return getRockVisibility(position);
}
bool GLSL::stoneVisibility(const vec2 &position)
{
    return getStoneVisibility(position);
}