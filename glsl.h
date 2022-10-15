#ifndef GLSL_H
#define GLSL_H

#include "vectorMath.h"
#include "noise.h"

using namespace vm;

namespace GLSL
{
    #define vec2(a, b) vec2{a, b}
    #define vec3(a, b, c) vec3{a, b, c}
    #define vec4(a, b, c, d) vec4{a, b, c, d}

    float humidityNoise(const vec2 &position);
    float temperatureNoise(const vec2 &position);
    float wetnessNoise(const vec2 &position);
    float stiffnessNoise(const vec2 &position);
    float desertNoise(const vec2 &position);
    float mountainNoise(const vec2 &position);
    float iceMountainNoise(const vec2 &position);
    float oceanNoise(const vec2 &position);
};

#endif // GLSL_H