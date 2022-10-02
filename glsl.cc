#include "glsl.h"

using namespace GLSL;

//

// * START GLSL ---------------------

const float NOISE_SCALE = 256.f;

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

float FBM(vec2 position)
{
    int octaves = 4;
    float frequency = 0.54f;
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

    return (noise * 0.5f) + 0.5f;
}

// ? Domain Warping : https://www.shadertoy.com/view/4s23zz
float warpNoise2D(vec2 position)
{
    vec2 q = vec2(FBM(position + vec2(0.0, 0.0)), FBM(position + vec2(27.2, 52.3)));
    return FBM(position + q * NOISE_SCALE * 2.f);
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