#include "glsl.h"

using namespace GLSL;

//

// * START GLSL ---------------------

const float NOISE_SCALE = 512.f;

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

// ? Domain Warping : https://www.shadertoy.com/view/4s23zz
float warpNoise2D(vec2 position)
{
    vec2 q = vec2(FBM_4(position + vec2(0.0f, 0.0f)), FBM_4(position + vec2(7.4f, 30.2f)));
    return FBM_4(position + q * NOISE_SCALE * 2.f);
}

float terrainNoise(vec2 position)
{
    float fbm2 = FBM_2(position);
    float fbm2d2 = FBM_2(position/2.f);
    float fbm2d3 = FBM_2(position/3.f);
    float fbm2d5 = FBM_2(position/5.f);

    float fbm4 = FBM_4(position);

    float noiseArea = clamp(fbm2d2 * 5.f, 0.f, 1.f);
    float noise2D = clamp(fbm2d3, 0.f, 1.f);
    float bigMountains = clamp(fbm4 , 0.f, 1.f);
    float detailedNoise = clamp(fbm2 *2.f - bigMountains/5.f + 0.7f, 0.f, 1.f);
    float cliffNoise = clamp(fbm2 *3.f - bigMountains/ 2.f + 0.6f, 0.f, 1.f);
    float smallHills = detailedNoise * (1.f - noiseArea) * 50.f;
    float mountains = (fbm2d3 * 300.f) * noiseArea ;

    float terrainBlend = clamp(fbm2d3 * 2.f, 0.f, 1.f);
    float terrainFlatter = clamp(fbm2d5 * 5.f, 0.2f + clamp(fbm2/5.f, 0.f, 1.f), 1.f);

    float terrain1 = (smallHills + mountains) * terrainBlend;
    // terrain1 = 0.f;
    float terrain2 = (cliffNoise + fbm2d2 * 500.f) * (1.f - terrainBlend);
    // terrain2 = 0.f;
    return (terrain1 + terrain2) * terrainFlatter;
}

// * END GLSL ---------------------

//

float GLSL::simplex2D(const vec2 &position)
{
    return glslNoise.in2D((double)position.x, (double)position.y);
}

float GLSL::in2DWarp(const vec2 &position)
{
    return warpNoise2D(position);
}

float GLSL::in2DTerrain(const vec2 &position)
{
    return terrainNoise(position);
}