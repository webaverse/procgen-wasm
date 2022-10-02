#include "glsl.h"

using namespace GLSL;

// 

// * START GLSL ---------------------

float FBM(vec2 position)
{
    int octaves = 4;
    float frequency = 0.54f;
    float lacunarity = 3.f;
    float persistence = 0.7f;

    float SCALE = 1.f / 128.f;
    vec2 p = position * SCALE;
    double noise = 0.f;

    float amplitude = 1.f;
    p *= frequency;

    for (int i = 0; i < octaves; i++)
    {
        noise += simplex2D(position) * amplitude;
        p *= lacunarity;
        amplitude *= persistence;
    }

    return noise;
}

// ? Domain Warping : https://www.shadertoy.com/view/4s23zz
float warpNoise2D(vec2 position)
{
    vec2 q = vec2(FBM(position + vec2(0.0, 0.0)), FBM(position + vec2(7.2, 2.3)));
    return FBM(position + q * 80.f);
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