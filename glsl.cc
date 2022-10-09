#include "glsl.h"
#include "constants.h"
#include "noises.h"

using namespace GLSL;

// ! GLSL INPUTS

const float MAX_TERRAIN_HEIGHT = float(MAX_WORLD_HEIGHT);
const float MIN_TERRAIN_HEIGHT = float(MIN_WORLD_HEIGHT);
const float TERRAIN_BASE_HEIGHT = float(WORLD_BASE_HEIGHT + 5); // moving the terrain 5 meters up ( to make it go above the water )
const float TERRAIN_OCEAN_THRESHOLD = OCEAN_THRESHOLD;

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
vec3 permute(vec3 x) { return mod(((x * 34.0) + 1.0) * x, 289.0); }
float snoise(vec2 v)
{
    const vec4 C = vec4(0.211324865405187, 0.366025403784439, -0.577350269189626, 0.024390243902439);
    vec2 i = floor(v + dot(v, vec2(C.y, C.y)));
    vec2 x0 = v - i + dot(i, vec2(C.x, C.x));
    vec2 i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
    vec4 x12 = vec4(x0.x, x0.y, x0.x, x0.y) + vec4(C.x, C.x, C.z, C.z) - vec4(i1.x, i1.y, 0.f, 0.f);
    i = mod(i, 289.0);
    vec3 p = permute(permute(vec3(0.0, i1.y, 1.0) + i.y) + i.x + vec3(0.0, i1.x, 1.0));
    vec3 m = max(vec3(0.5f, 0.5f, 0.5f) -
                     vec3(dot(x0, x0), dot(vec2(x12.x, x12.y), vec2(x12.x, x12.y)), dot(vec2(x12.z, x12.w), vec2(x12.z, x12.w))),
                 0.0);
    m = m * m;
    m = m * m;
    vec3 x = fract(p * vec3(C.w, C.w, C.w)) * 2.f - 1.0;
    vec3 h = abs(x) - 0.5;
    vec3 ox = floor(x + 0.5);
    vec3 a0 = x - ox;
    const float ff = 1.79284291400159;
    const float fd = 0.85373472095314;
    m *= vec3(ff, ff, ff) - (a0 * a0 + h * h) * fd;
    vec2 gf = vec2(a0.y, a0.z) * vec2(x12.x, x12.z) + vec2(h.y, h.z) * vec2(x12.y, x12.w);
    vec3 g = vec3(a0.x * x0.x + h.x * x0.y, gf.x, gf.y);
    return 135.f * dot(m, g);
}

vec2 rand2(vec2 p)
{
    return fract(vec2(sin(p.x * 591.32f + p.y * 154.077f), cos(p.x * 391.32f + p.y * 49.077f)));
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
            vec2 r = b - f + rand2(p + (b));
            float d = (r.x * r.x + r.y * r.y);
            res = min(res, d);
        }
    }
    return (res);
}

float FBM_8(vec2 position)
{
    int octaves = 8;
    float frequency = 0.6f;
    float lacunarity = 3.f;
    float persistence = 0.3f;

    float SCALE = 1.f / NOISE_SCALE;
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
    int octaves = 4;
    float frequency = 0.8f;
    float lacunarity = 3.f;
    float persistence = 0.3f;

    float SCALE = 1.f / NOISE_SCALE;
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
    int octaves = 2;
    float frequency = 0.74f;
    float lacunarity = 3.f;
    float persistence = 0.5f;

    float SCALE = 1.f / NOISE_SCALE;
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

float wetnessNoise(vec2 position)
{
    float noise = warpNoise1Layer_4(position);
    return clamp(noise + 0.65f, 0.f, 1.f);
}

float waterNoise(vec2 position)
{
    float noise = FBM_2(position/40.f);
    return clamp(1.f - noise, 0.f, 1.f);
}

float sandMountainNoise(vec2 position)
{
    const float desertMountainHeight = MAX_TERRAIN_HEIGHT / 20.f;
    // calculating noises
    float fbm2 = FBM_2(position / 2.f);

    float fbm4 = FBM_4(position / 2.f);
    float fbm4clamped = clamp(fbm4 , 0.f, 1.f);

    float fbm8 = FBM_8(position / 2.f);

    // layering noises
    float detailedNoise = clamp(fbm2 * 2.f - fbm4clamped / 4.f + 0.6f, 0.f, 1.f);
    float finalHeight = (detailedNoise - fbm8) * desertMountainHeight;
    
    return clamp(finalHeight, MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);
}

float snowNoise(vec2 position)
{
    const float iceMountainHeight = MAX_TERRAIN_HEIGHT / 20.f;
    // calculating noises
    float fbm2d2 = FBM_4(position);
    float fbm8d2 = FBM_8(position / 2.f);
    float fbm8d4 = FBM_8(position / 4.f);

    float vord300 = voronoi(position / 150.f);

    // layering noises
    float snowMountains = (2.f - (fbm2d2 + fbm8d2 + fbm8d4)*2.f) + vord300;

    float finalHeight = snowMountains * iceMountainHeight;

    return clamp(finalHeight, MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);
}

float valleyNoise(vec2 position)
{
    // calculating noises
    float fbm2d2 = FBM_4(position);
    float fbm8d2 = FBM_8(position / 2.f);
    float fbm8d4 = FBM_8(position / 4.f);

    float vord300 = voronoi(position / 150.f);

    // layering noises
    float snowMountains = (3.f - (fbm2d2 + fbm8d2 + fbm8d4)) + vord300/2.f;

    float finalHeight = snowMountains * (MAX_TERRAIN_HEIGHT/4.f);

    return clamp(finalHeight, MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);
}

float flattenerNoise(vec2 position){
    float fbm2d5 = FBM_2(position/5.f);
    float terrainFlattener = clamp(fbm2d5, 0.f, 0.5f) * 2.f;
    return terrainFlattener;
}

float terrainNoise(vec2 position)
{
    // defining terrain height parameters
    const float highMountainsHeight = MAX_TERRAIN_HEIGHT;
    const float lowMountainsHeight = MAX_TERRAIN_HEIGHT / 2.f;
    const float smallHillsHeight = MAX_TERRAIN_HEIGHT / 10.f;
    const float oceanDepth = MAX_TERRAIN_HEIGHT / 15.f;

    // calculating noises
    float fbm2 = FBM_2(position/2.f);

    float fbm2d2 = FBM_2(position/4.f);

    float fbm2d3 = FBM_2(position/6.f);
    float fbm2d3clamped = clamp(fbm2d3, 0.f, 1.f);


    float fbm4 = FBM_4(position/2.f);
    float fbm4clamped = clamp(fbm4 , 0.f, 1.f);

    float fbm8 = FBM_8(position/2.f);

    float warp2d2 = warpNoise1Layer_2(position/4.f);

    // layering noises

    float smallHills = clamp(fbm2 * 2.f - fbm4clamped / 5.f + 0.7f, 0.f, 1.f);
    float smallHillsLayer = smallHills * smallHillsHeight;
    float lowMountains = (fbm2d3 + warp2d2 / 8.f) * lowMountainsHeight;

    float smallHillsMountainsBlender = clamp(fbm2d2 * 5.f, 0.f, 1.f);
    float lowMountainsLayer = mix(smallHillsLayer, lowMountains, smallHillsMountainsBlender);

    float flatAreaSmallHills = clamp(fbm8, 0.f, 0.4f);
    float flatAreaLayer = flatAreaSmallHills * (lowMountainsHeight / 2.f);

    float lowAndHighMountainsBlender = fbm2d3clamped ;
    float highMountainsLayer = fbm2d2 * highMountainsHeight;
    float highAreasHeight = mix(highMountainsLayer, lowMountainsLayer, lowAndHighMountainsBlender);
    float highAreaLayer = highAreasHeight;

    float terrainFlattener = flattenerNoise(position);
    float terrainHeight = clamp(mix(flatAreaLayer, highAreaLayer, terrainFlattener), MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);

    float oceanNoise = waterNoise(position) * (1.f - terrainFlattener);
    float ocean = (softClamp(oceanNoise, TERRAIN_OCEAN_THRESHOLD - 0.1f, 1.f) - TERRAIN_OCEAN_THRESHOLD) * (1.f / (1.f - TERRAIN_OCEAN_THRESHOLD));
    float oceanLayer = ocean * oceanDepth;

    float finalHeight = TERRAIN_BASE_HEIGHT + terrainHeight - oceanLayer;

    return clamp(finalHeight, MIN_TERRAIN_HEIGHT, MAX_TERRAIN_HEIGHT);
}

// * END GLSL ---------------------

// --------------------------------

float GLSL::simplex2D(const vec2 &position)
{
    return glslNoise.in2D((double)position.x, (double)position.y);
}

// ! The name of the c++ functions should be different than the name of the glsl functions

float GLSL::humidityNoise(const vec2 &position)
{
    return wetnessNoise(position);
}

float GLSL::desertNoise(const vec2 &position)
{
    return sandMountainNoise(position);
}

float GLSL::mountainNoise(const vec2 &position)
{
    return terrainNoise(position);
}

float GLSL::iceMountainNoise(const vec2 &position)
{
    return snowNoise(position);
}

float GLSL::oceanNoise(const vec2 &position)
{
    float terrainFlattener = flattenerNoise(position);
    float oceanNoise = waterNoise(position) * (1.f - terrainFlattener);
    return oceanNoise;
}