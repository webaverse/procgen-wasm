#ifndef GLSL_H
#define GLSL_H

#include "vectorMath.h"
#include "noise.h"

using namespace vm;

namespace GLSL
{
    static Noise glslNoise(0, 0.02, 1);

    #define gl_FragColor(a) return a
    #define vec2(a, b) vec2{a, b}
    #define vec3(a, b, c) vec3{a, b, c}
    #define vec4(a, b, c, d) vec4{a, b, c, d}
    #define simplex2D(a, b) simplex2D(a, b)

    float simplex2D(float a, float b);
    float in2DWarp(const vec2 &position);
};

#endif // GLSL_H