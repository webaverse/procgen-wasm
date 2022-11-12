#include "glsl.h"
#include "constants.h"
#include "noises.h"

using namespace GLSL;

// ! GLSL INPUTS

const float MAX_TERRAIN_HEIGHT = float(MAX_WORLD_HEIGHT);
const float MIN_TERRAIN_HEIGHT = float(MIN_WORLD_HEIGHT);
const float TERRAIN_WATER_HEIGHT = float(WATER_BASE_HEIGHT);
const float WATER_SURROUNDING_HEIGHT = float(WORLD_BASE_HEIGHT / 4);
const float TERRAIN_WATER_DEPTH = float(WATER_BASE_HEIGHT * 8);
// ? making sure the terrain surface is above water before water depth subtraction
const float TERRAIN_BASE_HEIGHT = float(WATER_SURROUNDING_HEIGHT + WATER_BASE_HEIGHT); 
const float TERRAIN_FLATTENER_DEPTH = float(WATER_BASE_HEIGHT * 2);
const float TERRAIN_OCEAN_THRESHOLD = OCEAN_THRESHOLD;
const float TERRAIN_STONE_THRESHOLD = ROCK_THRESHOLD;

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
    return t != x ? t : b + (a - b)/(1. + exp((b - a)*(2.*x - a - b)/((x - a)*(b - x))));
}

// approaches the limit more slowly with slope <= 1 everywhere
float softClamp(float x, float a, float b)
{
 	float mid = (a + b)*0.5;
    return mid + smoothClamp((x - mid)*0.5, a - mid, b - mid);
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

float hash2D_1(vec2 p) {
  float h = dot(p, vec2(127.1, 311.7));
  return fract(sin(h) * 43758.5453123);
}

vec2 hash2D_2(vec2 p)
{
	p = vec2(dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)));
	return vec2(-1.f, -1.f) + fract(sin(p)*43758.5453123) * 2.f;
}

float snoise(vec2 p)
{
	vec2  i = floor(p + (p.x+p.y)*C.y);
    vec2  a = p - i + (i.x+i.y)*C.x;
    float m = step(a.y,a.x); 
    vec2  o = vec2(m,1.f - m);
    vec2  b = a - o + C.x;
	vec2  c = a - 1.0 + 2.0*C.x;
    vec3  h = max(vec3(-dot(a,a), -dot(b,b), -dot(c,c)) + 0.5f, 0.f);
	vec3  n = h*h*h*h*vec3(dot(a,hash2D_2(i+0.f)), dot(b,hash2D_2(i+o)), dot(c,hash2D_2(i+1.f)));
    return clamp(dot(n, vec3(70.f, 70.f, 70.f)), 0.f, 1.f);
}

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
            vec2 r = b - f + hash2D_2(p + (b));
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

    return (noise + 0.5f) * 0.5f;
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

    return (noise + 0.5f) * 0.5f;
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

    return (noise + 0.5f) * 0.5f;
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

    vec2 r = vec2(FBM_2(position + q * NOISE_SCALE * 2.f + vec2(1.7,9.2)),
                  FBM_2(position + q * NOISE_SCALE * 2.f + vec2(8.3,2.8)));

    return FBM_2(position + r * NOISE_SCALE * 2.f);
}

float warpNoise2Layer_4(vec2 position)
{
    vec2 q = vec2(FBM_4(position + vec2(0.0f, 0.0f)), 
                  FBM_4(position + vec2(7.4f, 30.2f)));

    vec2 r = vec2(FBM_4(position + q * NOISE_SCALE * 2.f + vec2(1.7,9.2)),
                  FBM_4(position + q * NOISE_SCALE * 2.f + vec2(8.3,2.8)));

    return FBM_4(position + r * NOISE_SCALE * 2.f);
}

// terrain noises
float flattenerNoise(vec2 position) 
{
    float simd5 = simplex(position/5.f);
    float flatAreas = (clamp(simd5, 0.f, 0.5f) * 2.f) * TERRAIN_FLATTENER_DEPTH;
    return flatAreas;
}
float waterNoise(vec2 position)
{
    float noise = clamp(simplex(position / 8.f) * 2.f, 0.f, 1.f);
    return noise;
}

float softWater(vec2 position)
{
    const float edge = TERRAIN_OCEAN_THRESHOLD;
    float oceanNoise = waterNoise(position);
    float softOcean = smoothClamp(oceanNoise, edge, 1.f);
    // making the value fit in the 0-1 range
    float ocean = ((softOcean - edge) * (1.f / (1.f - edge)));
    return ocean;
}

bool getWaterVisibility(vec2 position)
{
    const float edge = TERRAIN_OCEAN_THRESHOLD;
    float ocean = waterNoise(position);
    bool visibility = bool(step(edge, ocean));
    return visibility;
}

float waterPitNoise(vec2 position){
    float sWater = softWater(position);
    float waterPit = clamp(sWater, 0.f, 1.f);
    return waterPit * TERRAIN_WATER_DEPTH;
}

float wetNoise(vec2 position)
{
    return clamp(simplex(position / 2.f) * 2.f + 0.75f, 0.f, 1.f);
}

float grassObjNoise(vec2 position)
{
    if(getWaterVisibility(position)) {
        return 0.f;
    }
    float wetness = wetNoise(position);
    return clamp(simplex(position * 40.f) * wetness, 0.f, 1.f);
}

float treeObjNoise(vec2 position)
{
    if(getWaterVisibility(position)) {
        return 0.f;
    }
    float wetness = wetNoise(position);
    return clamp(simplex(position * 20.f) * wetness, 0.f, 1.f);
}

float grassMatNoise(vec2 position)
{
    float wetness = clamp(wetNoise(position) + 0.25f, 0.f, 1.f);
    float noise = warpNoise1Layer_2(position * 10.f) * wetness;
    return clamp(noise * 2.f - 0.5f, 0.f, 1.f);
}

bool getStoneVisibility(vec2 position)
{
    if(getWaterVisibility(position)) {
        return 0.f;
    }
    const float edge = TERRAIN_STONE_THRESHOLD;
    float dryness = 1.f - (wetNoise(position) - 0.5f);
    float stone = clamp(simplex(position * 6.f) * dryness, 0.f, 1.f);
    bool visibility = bool(step(edge, stone));
    return visibility;
}

float stiffNoise(vec2 position)
{
    return warpNoise1Layer_1(position * 5.f);
}

float humNoise(vec2 position)
{
    float noise = simplex(position/18.f);
    return clamp(noise + 0.15f, 0.f, 1.f);
}

float heatNoise(vec2 position)
{
    return simplex(position/30.f);
}

// terrain height noises

float terrainHeightWrapper(vec2 position, float height)
{
    float flatAreas = flattenerNoise(position);
    float flattenedHeight = TERRAIN_BASE_HEIGHT + height - flatAreas;
    float flattenedHeightClamped = clamp(flattenedHeight, TERRAIN_BASE_HEIGHT, MAX_TERRAIN_HEIGHT);
    float waterPit = waterPitNoise(position);
    float finalHeight = flattenedHeightClamped - waterPit;
    return clamp(finalHeight, MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);
}

float sandMountainNoise(vec2 position)
{
    const float desertMountainHeight = MAX_TERRAIN_HEIGHT / 10.f;

    // calculating noises
    float fbm2d2 = FBM_2(position / 2.f);

    float fbm4d2 = FBM_4(position / 2.f);
    float fbm4clamped = clamp(fbm4d2 , 0.f, 1.f);

    float fbm8d2 = FBM_8(position / 2.f);

    // layering noises

    float detailedNoise = clamp(fbm2d2 * 2.f - fbm4clamped / 4.f + 0.6f, 0.f, 1.f);
    float sandMountain = (detailedNoise - fbm8d2) * desertMountainHeight;
    
    return terrainHeightWrapper(position, sandMountain);
}

float snowNoise(vec2 position)
{
    const float iceMountainHeight = MAX_TERRAIN_HEIGHT / 4.f;

    // calculating noises
    float fbm4 = FBM_4(position);
    float fbm8d2 = FBM_8(position / 2.f);
    float fbm2d4 = FBM_2(position / 4.f);

    float vor = voronoiNoise(position);

    // layering noises
    float snowMountains = clamp(((2.f - (fbm4 + fbm8d2 + fbm2d4)*2.f) + vor) / 4.f, 0.f, 1.f) * iceMountainHeight;
    return terrainHeightWrapper(position, snowMountains);
}

// float valleyNoise(vec2 position)
// {
//     // calculating noises
//     float fbm2d2 = FBM_4(position);
//     float fbm8d2 = FBM_8(position / 2.f);
//     float fbm8d4 = FBM_8(position / 4.f);

//     float vord300 = voronoi(position / 150.f);

//     // layering noises
//     float snowMountains = (3.f - (fbm2d2 + fbm8d2 + fbm8d4)) + vord300/2.f;

//     float finalHeight = snowMountains * (MAX_TERRAIN_HEIGHT/4.f);

//     return clamp(finalHeight, MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);
// }

float mountainHillsNoise(vec2 position)
{
    // defining terrain height parameters
    const float highMountainsHeight = MAX_TERRAIN_HEIGHT;
    const float lowMountainsHeight = MAX_TERRAIN_HEIGHT / 2.f;
    const float smallHillsHeight = MAX_TERRAIN_HEIGHT / 3.f;

    // calculating noises
    float fbm2d2 = FBM_2(position/2.f);
    float fbm2d2clamped = clamp(fbm2d2 , 0.f, 1.f);

    float fbm2d4 = FBM_2(position/4.f);

    float fbm2d6 = FBM_2(position/6.f);
    float fbm2d6clamped = clamp(fbm2d6, 0.f, 1.f);
    float fbm2d6clampedm2  = clamp(fbm2d6, 0.f, 0.5f) * 2.f;

    float warp2d4 = warpNoise1Layer_1(position/4.f);

    // layering noises

    float smallHills = clamp(fbm2d2 * 2.f - fbm2d2clamped  / 5.f + 0.7f, 0.f, 1.f);
    float smallHillsLayer = smallHills * smallHillsHeight;
    float lowMountains = (fbm2d6 + warp2d4 / 5.f) * lowMountainsHeight;

    float smallHillsMountainsBlender = clamp(fbm2d4 * 5.f, 0.f, 1.f);
    float lowMountainsLayer = mix(smallHillsLayer, lowMountains, smallHillsMountainsBlender);

    float flatAreaSmallHills = clamp(fbm2d2, 0.f, 0.4f);
    float flatAreaLayer = flatAreaSmallHills * (lowMountainsHeight / 2.f);

    float lowAndHighMountainsBlender = fbm2d6clamped ;
    float highMountainsLayer = (fbm2d4 - warp2d4 / 7.f) * highMountainsHeight;
    float highAreasHeight = mix(highMountainsLayer, lowMountainsLayer, lowAndHighMountainsBlender);
    float highAreaLayer = highAreasHeight;

    float highAndFlatAreaBlender = fbm2d6clampedm2;
    float mountainHeight = mix(highAreaLayer, flatAreaLayer, highAndFlatAreaBlender);
    
    return terrainHeightWrapper(position, mountainHeight);
}

// * END GLSL ---------------------

// --------------------------------

// ! The name of the c++ functions should be different than the name of the glsl functions

float GLSL::humidityNoise(const vec2 &position)
{
    return humNoise(position);
}

float GLSL::grassMaterialNoise(const vec2 &position)
{
    return grassMatNoise(position);
}

float GLSL::grassObjectNoise(const vec2 &position)
{
    return grassObjNoise(position);
}

float GLSL::treeObjectNoise(const vec2 &position)
{
    return treeObjNoise(position);
}

float GLSL::stiffnessNoise(const vec2 &position)
{
    return stiffNoise(position);
}

float GLSL::temperatureNoise(const vec2 &position)
{
    return heatNoise(position);
}

float GLSL::desertNoise(const vec2 &position)
{
    return sandMountainNoise(position);
}

float GLSL::mountainNoise(const vec2 &position)
{
    return mountainHillsNoise(position);
}

float GLSL::iceMountainNoise(const vec2 &position)
{
    return snowNoise(position);
}

float GLSL::hashNoise(const vec2 &position)
{
    return hash2D_1(position);
}

float GLSL::simplexNoise(const vec2 &position)
{
    return simplex(position);
}

bool GLSL::oceanNoise(const vec2 &position)
{
    return getWaterVisibility(position);
}

bool GLSL::stoneNoise(const vec2 &position)
{
    return getStoneVisibility(position);
}