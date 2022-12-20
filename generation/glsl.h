#ifndef GLSL_H
#define GLSL_H

#include "../libs/vectorMath.h"
#include "noise.h"

using namespace vm;

namespace GLSL
{
    // static Noise glslNoise(0, 0.01, 1);

    #define vec2(a, b) vec2{a, b}
    #define vec3(a, b, c) vec3{a, b, c}
    #define vec4(a, b, c, d) vec4{a, b, c, d}

    // float snoise(vec2 position);

    float humidityNoise(const vec2 &position);
    float temperatureNoise(const vec2 &position);
    float grassMaterialNoise(const vec2 &position, const float &wetness);
    float wetnessNoise(const vec2 &position);
    float coldnessNoise(const vec2 &position);
    float desertNoise(const vec2 &position);
    float mountainNoise(const vec2 &position);
    float iceMountainNoise(const vec2 &position);
    float hashNoise(const vec2 &position);
    float simplexNoise(const vec2 &position);
    float waterDepthNoise(const vec2 &position);
    float oceanNoise(const vec2 &position);
    float riverNoise(const vec2 &position, const float &ocean);

    vm::vec3 flowNoise(const vec2 &position);

    bool flowerVisibility(const vec2 &position, const float &grass);
    bool grassVisibility(const vec2 &position, const float &wetness);
    bool treeVisibility(const vec2 &position, const float &wetness);
    bool waterVisibilityNoise(const vec2 &position);
    bool rockVisibility(const vec2 &position);
    bool stoneVisibility(const vec2 &position);
};

#endif // GLSL_H