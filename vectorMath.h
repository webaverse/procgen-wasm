#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#include <cmath>
#include <algorithm>
#include <cstdio>

namespace vm
{
    struct vec2
    {
        // vec2()
        // {
        //     x = 0;
        //     y = 0;
        // }
        // vec2(const float &_x, const float &_y) : x(_x), y(_y) {}
        // vec2(const float &_v) : x(_v), y(_v) {}
        float x;
        float y;
    };
    struct vec3
    {
        // vec3()
        // {
        //     x = 0;
        //     y = 0;
        //     z = 0;
        // }
        // vec3(const float &_x, const float &_y, const float &_z) : x(_x), y(_y), z(_z) {}
        // vec3(const float &_v) : x(_v), y(_v), z(_v) {}
        float x;
        float y;
        float z;
    };
    struct vec4
    {
        // vec4()
        // {
        //     x = 0;
        //     y = 0;
        //     z = 0;
        //     w = 0;
        // }
        // vec4(const float &_x, const float &_y, const float &_z, const float &_w) : x(_x), y(_y), z(_z), w(_w) {}
        // vec4(const float &_v) : x(_v), y(_v), z(_v), w(_v) {}
        float x;
        float y;
        float z;
        float w;
    };
    struct ivec2
    {
        // ivec2()
        // {
        //     x = 0;
        //     y = 0;
        // }
        // ivec2(const int &_x, const int &_y) : x(_x), y(_y) {}
        // ivec2(const int &_v) : x(_v), y(_v) {}
        int x;
        int y;
    };
    struct ivec3
    {
        // ivec3()
        // {
        //     x = 0;
        //     y = 0;
        //     z = 0;
        // }
        // ivec3(const int &_x, const int &_y, const int &_z) : x(_x), y(_y), z(_z) {}
        // ivec3(const int _v) : x(_v), y(_v), z(_v) {}
        int x;
        int y;
        int z;
    };
    struct ivec4
    {
        // ivec4()
        // {
        //     x = 0;
        //     y = 0;
        //     z = 0;
        //     w = 0;
        // }
        // ivec4(const int &_x, const int &_y, const int &_z, const int &_w) : x(_x), y(_y), z(_z), w(_w) {}
        // ivec4(const int &_v) : x(_v), y(_v), z(_v), w(_v) {}
        int x;
        int y;
        int z;
        int w;
    };
    struct ibox3
    {
        // ibox3()
        // {
        //     min = ivec3(0);
        //     max = ivec3(0);
        // }
        // ibox3(const ivec3 &_min, const ivec3 &_max) : min(_min), max(_max) {}
        ivec3 min;
        ivec3 max;
    };
    struct box3
    {
        // box3()
        // {
        //     min = vec3(0);
        //     max = vec3(0);
        // }
        // box3(const vec3 &_min, const vec3 &_max) : min(_min), max(_max) {}
        vec3 min;
        vec3 max;
    };

    float min(const float &v1, const float &v2);
    vec2 min(const vec2 &v1, const vec2 &v2);
    vec3 min(const vec3 &v1, const vec3 &v2);
    vec4 min(const vec4 &v1, const vec4 &v2);
    ivec2 min(const ivec2 &v1, const ivec2 &v2);
    ivec3 min(const ivec3 &v1, const ivec3 &v2);
    ivec4 min(const ivec4 &v1, const ivec4 &v2);

    vec2 min(const vec2 &v1, const float _m);
    vec3 min(const vec3 &v1, const float _m);
    vec4 min(const vec4 &v1, const float _m);
    ivec2 min(const ivec2 &v1, const int _m);
    ivec3 min(const ivec3 &v1, const int _m);
    ivec4 min(const ivec4 &v1, const int _m);

    vec2 max(const vec2 &v1, const vec2 &v2);
    vec3 max(const vec3 &v1, const vec3 &v2);
    vec4 max(const vec4 &v1, const vec4 &v2);
    ivec2 max(const ivec2 &v1, const ivec2 &v2);
    ivec3 max(const ivec3 &v1, const ivec3 &v2);
    ivec4 max(const ivec4 &v1, const ivec4 &v2);

    vec2 max(const vec2 &v1, const float _m);
    vec3 max(const vec3 &v1, const float _m);
    vec4 max(const vec4 &v1, const float _m);
    ivec2 max(const ivec2 &v1, const int _m);
    ivec3 max(const ivec3 &v1, const int _m);
    ivec4 max(const ivec4 &v1, const int _m);

    vec2 abs(const vec2 &v);
    vec3 abs(const vec3 &v);
    vec4 abs(const vec4 &v);
    ivec2 abs(const ivec2 &v);
    ivec3 abs(const ivec3 &v);
    ivec4 abs(const ivec4 &v);

    vec2 normalize(const vec2 &v);
    vec3 normalize(const vec3 &v);
    vec4 normalize(const vec4 &v);
    ivec2 normalize(const ivec2 &v);
    ivec3 normalize(const ivec3 &v);
    ivec4 normalize(const ivec4 &v);

    float step(float edge, float x);
    vec2 step(vec2 edge, vec2 x);
    vec3 step(vec3 edge, vec3 x);
    vec4 step(vec4 edge, vec4 x);

    vec2 step(float edge, vec2 x);
    vec3 step(float edge, vec3 x);
    vec4 step(float edge, vec4 x);

    float length(const vec2 v);
    float length(const vec3 v);
    float length(const vec4 v);
    float length(const ivec2 v);
    float length(const ivec3 v);
    float length(const ivec4 v);
    float length(const float v);
    float lengthSq(const ivec3 &v);
    float lengthSq(const vec3 &v);

    float distance(const vm::vec2 &v, const vm::vec2 &o);
    float distance(const vm::vec2 &v, const vm::ivec2 &o);
    float distance(const vm::ivec2 &v, const vm::ivec2 &o);

    float dot(const vm::vec2 &v, const vm::vec2 &o);
    float dot(const vm::vec3 &v, const vm::vec3 &o);

    vec2 mod(const vec2 v, const float m);
    vec3 mod(const vec3 v, const float m);
    vec4 mod(const vec4 v, const float m);
    ivec2 mod(const ivec2 v, const float m);
    ivec3 mod(const ivec3 v, const float m);
    ivec4 mod(const ivec4 v, const float m);

    float mix(const float v1, const float v2, const float m);
    vec2 mix(const vec2 v1, const vec2 v2, const float m);
    vec3 mix(const vec3 v1, const vec3 v2, const float m);
    vec4 mix(const vec4 v1, const vec4 v2, const float m);
    vec2 mix(const ivec2 v1, const ivec2 v2, const float m);
    vec3 mix(const ivec3 v1, const ivec3 v2, const float m);
    vec4 mix(const ivec4 v1, const ivec4 v2, const float m);

    vec2 floor(const vec2 v);
    vec3 floor(const vec3 v);
    vec4 floor(const vec4 v);

    vec2 fract(const vec2 v);
    vec3 fract(const vec3 v);
    vec4 fract(const vec4 v);

    template <typename T>
    T clamp(const T v, const T min, const T max) {
       return std::max(min, std::min(v, max));
    }
    // template<>
    // double clamp(const double v, const double min, const double max) = delete;

    // + operator
    vec2 operator+(const vec2 &v1, const vec2 &v2);
    vec3 operator+(const vec3 &v1, const vec3 &v2);
    vec4 operator+(const vec4 &v1, const vec4 &v2);
    ivec2 operator+(const ivec2 &v1, const ivec2 &v2);
    ivec3 operator+(const ivec3 &v1, const ivec3 &v2);
    ivec4 operator+(const ivec4 &v1, const ivec4 &v2);

    vec2 operator+(const vec2 &v1, const float _m);
    vec3 operator+(const vec3 &v1, const float _m);
    vec4 operator+(const vec4 &v1, const float _m);
    ivec2 operator+(const ivec2 &v1, const int _m);
    ivec3 operator+(const ivec3 &v1, const int _m);
    ivec4 operator+(const ivec4 &v1, const int _m);
    // - operator
    vec2 operator-(const vec2 &v1, const vec2 &v2);
    vec3 operator-(const vec3 &v1, const vec3 &v2);
    vec4 operator-(const vec4 &v1, const vec4 &v2);
    ivec2 operator-(const ivec2 &v1, const ivec2 &v2);
    ivec3 operator-(const ivec3 &v1, const ivec3 &v2);
    ivec4 operator-(const ivec4 &v1, const ivec4 &v2);

    vec2 operator-(const vec2 &v1, const float _m);
    vec3 operator-(const vec3 &v1, const float _m);

    vec4 operator-(const vec4 &v1, const float _m);
    ivec2 operator-(const ivec2 &v1, const int _m);
    ivec3 operator-(const ivec3 &v1, const int _m);
    ivec4 operator-(const ivec4 &v1, const int _m);

    // * operator
    vec2 operator*(const vec2 &v1, const vec2 &v2);
    vec3 operator*(const vec3 &v1, const vec3 &v2);
    vec4 operator*(const vec4 &v1, const vec4 &v2);
    ivec2 operator*(const ivec2 &v1, const ivec2 &v2);
    ivec3 operator*(const ivec3 &v1, const ivec3 &v2);
    ivec4 operator*(const ivec4 &v1, const ivec4 &v2);
    vec2 operator*(const vec2 &v1, const float _m);
    vec3 operator*(const vec3 &v1, const float _m);
    vec4 operator*(const vec4 &v1, const float _m);
    ivec2 operator*(const ivec2 &v1, const int _m);
    ivec3 operator*(const ivec3 &v1, const int _m);
    ivec4 operator*(const ivec4 &v1, const int _m);

    // / operator
    vec2 operator/(const vec2 &v1, const vec2 &v2);
    vec3 operator/(const vec3 &v1, const vec3 &v2);
    vec4 operator/(const vec4 &v1, const vec4 &v2);
    ivec2 operator/(const ivec2 &v1, const ivec2 &v2);
    ivec3 operator/(const ivec3 &v1, const ivec3 &v2);
    ivec4 operator/(const ivec4 &v1, const ivec4 &v2);
    vec2 operator/(const vec2 &v1, const float _m);
    vec3 operator/(const vec3 &v1, const float _m);
    vec4 operator/(const vec4 &v1, const float _m);
    ivec2 operator/(const ivec2 &v1, const int _m);
    ivec3 operator/(const ivec3 &v1, const int _m);
    ivec4 operator/(const ivec4 &v1, const int _m);

    // += operator
    vec2 operator+=(vec2 &v1, const vec2 &v2);
    vec3 operator+=(vec3 &v1, const vec3 &v2);
    vec4 operator+=(vec4 &v1, const vec4 &v2);
    ivec2 operator+=(ivec2 &v1, const ivec2 &v2);
    ivec3 operator+=(ivec3 &v1, const ivec3 &v2);
    ivec4 operator+=(ivec4 &v1, const float _m);
    vec2 operator+=(vec2 &v1, const float _m);
    vec3 operator+=(vec3 &v1, const float _m);
    vec4 operator+=(vec4 &v1, const float _m);
    ivec2 operator+=(ivec2 &v1, const int _m);
    ivec3 operator+=(ivec3 &v1, const int _m);
    ivec4 operator+=(ivec4 &v1, const int _m);
    // -= operator
    vec2 operator-=(vec2 &v1, const vec2 &v2);
    vec3 operator-=(vec3 &v1, const vec3 &v2);
    vec4 operator-=(vec4 &v1, const vec4 &v2);
    ivec2 operator-=(ivec2 &v1, const ivec2 &v2);
    ivec3 operator-=(ivec3 &v1, const ivec3 &v2);
    ivec4 operator-=(ivec4 &v1, const float _m);
    vec2 operator-=(vec2 &v1, const float _m);
    vec3 operator-=(vec3 &v1, const float _m);
    vec4 operator-=(vec4 &v1, const float _m);
    ivec2 operator-=(ivec2 &v1, const int _m);
    ivec3 operator-=(ivec3 &v1, const int _m);
    ivec4 operator-=(ivec4 &v1, const int _m);
    // *= operator
    vec2 operator*=(vec2 &v1, const vec2 &v2);
    vec3 operator*=(vec3 &v1, const vec3 &v2);
    vec4 operator*=(vec4 &v1, const vec4 &v2);
    ivec2 operator*=(ivec2 &v1, const ivec2 &v2);
    ivec3 operator*=(ivec3 &v1, const ivec3 &v2);
    ivec4 operator*=(ivec4 &v1, const float _m);
    vec2 operator*=(vec2 &v1, const float _m);
    vec3 operator*=(vec3 &v1, const float _m);
    vec4 operator*=(vec4 &v1, const float _m);
    ivec2 operator*=(ivec2 &v1, const int _m);
    ivec3 operator*=(ivec3 &v1, const int _m);
    ivec4 operator*=(ivec4 &v1, const int _m);
    // /= operator
    vec2 operator/=(vec2 &v1, const vec2 &v2);
    vec3 operator/=(vec3 &v1, const vec3 &v2);
    vec4 operator/=(vec4 &v1, const vec4 &v2);
    ivec2 operator/=(ivec2 &v1, const ivec2 &v2);
    ivec3 operator/=(ivec3 &v1, const ivec3 &v2);
    ivec4 operator/=(ivec4 &v1, const float _m);
    vec2 operator/=(vec2 &v1, const float _m);
    vec3 operator/=(vec3 &v1, const float _m);
    vec4 operator/=(vec4 &v1, const float _m);
    ivec2 operator/=(ivec2 &v1, const int _m);
    ivec3 operator/=(ivec3 &v1, const int _m);
    ivec4 operator/=(ivec4 &v1, const int _m);
    // == operator
    bool operator==(const vec2 &v1, const vec2 &v2);
    bool operator==(const vec3 &v1, const vec3 &v2);
    bool operator==(const vec4 &v1, const vec4 &v2);
    bool operator==(const ivec2 &v1, const ivec2 &v2);
    bool operator==(const ivec3 &v1, const ivec3 &v2);
    bool operator==(const ivec4 &v1, const ivec4 &v2);

    bool operator!=(const vec2 &v1, const vec2 &v2);
    bool operator!=(const vec3 &v1, const vec3 &v2);
    bool operator!=(const vec4 &v1, const vec4 &v2);
    bool operator!=(const ivec2 &v1, const ivec2 &v2);
    bool operator!=(const ivec3 &v1, const ivec3 &v2);
    bool operator!=(const ivec4 &v1, const ivec4 &v2);

    // % operator
    ivec2 operator%(const ivec2 &v1, const ivec2 &v2);
    ivec3 operator%(const ivec3 &v1, const ivec3 &v2);
    ivec4 operator%(const ivec4 &v1, const ivec4 &v2);
    ivec2 operator%(const ivec2 &v1, const int _m);
    ivec3 operator%(const ivec3 &v1, const int _m);
    ivec4 operator%(const ivec4 &v1, const int _m);

} // namespace vm

#endif // VECTOR_MATH_H