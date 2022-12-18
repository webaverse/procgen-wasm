#include "vectorMath.h"

float vm::min(const float &v1, const float &v2)
{
    return std::min(v1, v2);
};

vm::vec2 vm::min(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return vm::vec2{std::min(v1.x, v2.x), std::min(v1.y, v2.y)};
};
vm::vec3 vm::min(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return vm::vec3{std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z)};
};
vm::vec4 vm::min(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return vm::vec4{std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z), std::min(v1.w, v2.w)};
};
vm::ivec2 vm::min(const vm::ivec2 &v1, const vm::ivec2 &v2)
{
    return vm::ivec2{std::min(v1.x, v2.x), std::min(v1.y, v2.y)};
};
vm::ivec3 vm::min(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z)};
};
vm::ivec4 vm::min(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z), std::min(v1.w, v2.w)};
};

vm::vec2 vm::min(const vm::vec2 &v1, const float _m)
{
    return vm::vec2{std::min(v1.x, _m), std::min(v1.y, _m)};
};
vm::vec3 vm::min(const vm::vec3 &v1, const float _m)
{
    return vm::vec3{std::min(v1.x, _m), std::min(v1.y, _m), std::min(v1.z, _m)};
};
vm::vec4 vm::min(const vm::vec4 &v1, const float _m)
{
    return vm::vec4{std::min(v1.x, _m), std::min(v1.y, _m), std::min(v1.z, _m), std::min(v1.w, _m)};
};
vm::ivec2 vm::min(const vm::ivec2 &v1, const int _m)
{
    return vm::ivec2{std::min(v1.x, _m), std::min(v1.y, _m)};
};
vm::ivec3 vm::min(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{std::min(v1.x, _m), std::min(v1.y, _m), std::min(v1.z, _m)};
};
vm::ivec4 vm::min(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{std::min(v1.x, _m), std::min(v1.y, _m), std::min(v1.z, _m), std::min(v1.w, _m)};
};

vm::vec2 vm::max(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return vm::vec2{std::max(v1.x, v2.x), std::max(v1.y, v2.y)};
};
vm::vec3 vm::max(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return vm::vec3{std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z)};
};
vm::vec4 vm::max(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return vm::vec4{std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z), std::max(v1.w, v2.w)};
};
vm::ivec2 vm::max(const vm::ivec2 &v1, const vm::ivec2 &v2)
{
    return vm::ivec2{std::max(v1.x, v2.x), std::max(v1.y, v2.y)};
};
vm::ivec3 vm::max(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z)};
};
vm::ivec4 vm::max(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z), std::max(v1.w, v2.w)};
};

vm::vec2 vm::max(const vm::vec2 &v1, const float _m)
{
    return vm::vec2{std::max(v1.x, _m), std::max(v1.y, _m)};
};
vm::vec3 vm::max(const vm::vec3 &v1, const float _m)
{
    return vm::vec3{std::max(v1.x, _m), std::max(v1.y, _m), std::max(v1.z, _m)};
};
vm::vec4 vm::max(const vm::vec4 &v1, const float _m)
{
    return vm::vec4{std::max(v1.x, _m), std::max(v1.y, _m), std::max(v1.z, _m), std::max(v1.w, _m)};
};
vm::ivec2 vm::max(const vm::ivec2 &v1, const int _m)
{
    return vm::ivec2{std::max(v1.x, _m), std::max(v1.y, _m)};
};
vm::ivec3 vm::max(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{std::max(v1.x, _m), std::max(v1.y, _m), std::max(v1.z, _m)};
};
vm::ivec4 vm::max(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{std::max(v1.x, _m), std::max(v1.y, _m), std::max(v1.z, _m), std::max(v1.w, _m)};
};

vm::vec2 vm::abs(const vm::vec2 &v)
{
    return vm::vec2{std::abs(v.x), std::abs(v.y)};
}
vm::vec3 vm::abs(const vm::vec3 &v)
{
    return vm::vec3{std::abs(v.x), std::abs(v.y), std::abs(v.z)};
}
vm::vec4 vm::abs(const vm::vec4 &v)
{
    return vm::vec4{std::abs(v.x), std::abs(v.y), std::abs(v.z), std::abs(v.w)};
}
vm::ivec2 vm::abs(const vm::ivec2 &v)
{
    return vm::ivec2{std::abs(v.x), std::abs(v.y)};
}
vm::ivec3 vm::abs(const vm::ivec3 &v)
{
    return vm::ivec3{std::abs(v.x), std::abs(v.y), std::abs(v.z)};
}
vm::ivec4 vm::abs(const vm::ivec4 &v)
{
    return vm::ivec4{std::abs(v.x), std::abs(v.y), std::abs(v.z), std::abs(v.w)};
}

vm::vec2 vm::normalize(const vm::vec2 &v)
{
    float length = vm::length(v);

    if (!length) {
        return vm::vec2{0.f, 0.f};
    }

    return vm::vec2{v.x / length, v.y / length};
}
vm::vec3 vm::normalize(const vm::vec3 &v)
{
    float length = vm::length(v);

    if (!length) {
        return vm::vec3{0.f, 0.f, 0.f};
    }

    return vm::vec3{v.x / length, v.y / length, v.z / length};
}
vm::vec4 vm::normalize(const vm::vec4 &v)
{
    float length = vm::length(v);

    if (!length) {
        return vec4{0.f, 0.f, 0.f, 0.f};
    }

    return vm::vec4{v.x / length, v.y / length, v.z / length, v.w / length};
}
vm::ivec2 vm::normalize(const vm::ivec2 &v)
{
    float length = vm::length(v);

    if (!length) {
        return vm::ivec2{(int)0, (int)0};
    }

    return vm::ivec2{(int)(v.x / length), (int)(v.y / length)};
}
vm::ivec3 vm::normalize(const vm::ivec3 &v)
{
    float length = vm::length(v);

    if (!length) {
        return vm::ivec3{(int)0, (int)0, (int)0};
    }

    return vm::ivec3{(int)(v.x / length), (int)(v.y / length), (int)(v.z / length)};
}
vm::ivec4 vm::normalize(const vm::ivec4 &v)
{
    float length = vm::length(v);

    if (!length) {
        return vm::ivec4{(int)0, (int)0, (int)0, (int)0};
    }

    return vm::ivec4{(int)(v.x / length), (int)(v.y / length), (int)(v.z / length), (int)(v.w / length)};
}

float vm::step(float edge, float x)
{
    return x < edge ? 0.0f : 1.0f;
};
vm::vec2 step(vm::vec2 edge, vm::vec2 x)
{
    return vm::vec2{vm::step(edge.x, x.x), vm::step(edge.y, x.y)};
};
vm::vec3 step(vm::vec3 edge, vm::vec3 x)
{
    return vm::vec3{vm::step(edge.x, x.x), vm::step(edge.y, x.y), vm::step(edge.z, x.z)};
};
vm::vec4 step(vm::vec4 edge, vm::vec4 x)
{
    return vm::vec4{vm::step(edge.x, x.x), vm::step(edge.y, x.y), vm::step(edge.z, x.z), vm::step(edge.w, x.w)};
};

vm::vec2 step(float edge, vm::vec2 x)
{
    return vm::vec2{vm::step(edge, x.x), vm::step(edge, x.y)};
};
vm::vec3 step(float edge, vm::vec3 x)
{
    return vm::vec3{vm::step(edge, x.x), vm::step(edge, x.y), vm::step(edge, x.z)};
};
vm::vec4 step(float edge, vm::vec4 x)
{
    return vm::vec4{vm::step(edge, x.x), vm::step(edge, x.y), vm::step(edge, x.z), vm::step(edge, x.w)};
};

float vm::length(const vm::vec2 v)
{
    return std::sqrt(v.x * v.x + v.y * v.y);
};
float vm::length(const vm::vec3 v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
};
float vm::length(const vm::vec4 v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
};
float vm::length(const vm::ivec2 v)
{
    return std::sqrt(v.x * v.x + v.y * v.y);
};
float vm::length(const vm::ivec3 v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
};
float vm::length(const vm::ivec4 v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
};
float vm::length(const float v)
{
    return std::sqrt(v * v);
};
float vm::lengthSq(const vm::ivec3 &v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
float vm::lengthSq(const vm::vec3 &v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

float vm::distance(const vm::vec2 &v, const vm::vec2 &o)
{
    const float dx = v.x - o.x;
    const float dy = v.y - o.y;
    return std::sqrt(dx * dx + dy * dy);
}
float vm::distance(const vm::vec2 &v, const vm::ivec2 &o)
{
    const float dx = v.x - o.x;
    const float dy = v.y - o.y;
    return std::sqrt(dx * dx + dy * dy);
}
float vm::distance(const vm::ivec2 &v, const vm::ivec2 &o)
{
    const int dx = v.x - o.x;
    const int dy = v.y - o.y;
    return std::sqrt(dx * dx + dy * dy);
}

float vm::dot(const vm::vec2 &v, const vm::vec2 &o)
{
    return v.x * o.x + v.y * o.y;
}
float vm::dot(const vm::vec3 &v, const vm::vec3 &o)
{
    return v.x * o.x + v.y * o.y + v.z * o.z;
}

vm::vec3 vm::cross(const vm::vec3 &v, const vm::vec3 &o)
{
    return vm::vec3{v.y * o.z - v.z * o.y, v.z * o.x - v.x * o.z, v.x * o.y - v.y * o.x};
}

vm::vec2 vm::sin(const vm::vec2 &v)
{
    return vm::vec2{std::sin(v.x), std::sin(v.y)};
};

vm::vec3 vm::sin(const vm::vec3 &v)
{
    return vm::vec3{std::sin(v.x), std::sin(v.y), std::sin(v.z)};
};

vm::vec4 vm::sin(const vm::vec4 &v)
{
    return vm::vec4{std::sin(v.x), std::sin(v.y), std::sin(v.z), std::sin(v.w)};
};

vm::vec2 vm::sin(const vm::ivec2 &v)
{
    return vm::vec2{std::sin((float)v.x), std::sin((float)v.y)};
};

vm::vec3 vm::sin(const vm::ivec3 &v)
{
    return vm::vec3{std::sin((float)v.x), std::sin((float)v.y), std::sin((float)v.z)};
};

vm::vec4 vm::sin(const vm::ivec4 &v)
{
    return vm::vec4{std::sin((float)v.x), std::sin((float)v.y), std::sin((float)v.z), std::sin((float)v.w)};
};

vm::vec2 vm::mod(const vec2 v, const float m)
{
    return vm::vec2{(float)((int)v.x % (int)m), (float)((int)v.y % (int)m)};
};
vm::vec3 vm::mod(const vec3 v, const float m)
{
    return vm::vec3{(float)((int)v.x % (int)m), (float)((int)v.y % (int)m), (float)((int)v.z % (int)m)};
};
vm::vec4 vm::mod(const vec4 v, const float m)
{
    return vm::vec4{(float)((int)v.x % (int)m), (float)((int)v.y % (int)m), (float)((int)v.z % (int)m), (float)((int)v.w % (int)m)};
};
vm::ivec2 vm::mod(const ivec2 v, const float m)
{
    return vm::ivec2{(int)v.x % (int)m, (int)v.y % (int)m};
};
vm::ivec3 vm::mod(const ivec3 v, const float m)
{
    return vm::ivec3{(int)v.x % (int)m, (int)v.y % (int)m, (int)v.z % (int)m};
};
vm::ivec4 vm::mod(const ivec4 v, const float m)
{
    return vm::ivec4{(int)v.x % (int)m, (int)v.y % (int)m, (int)v.z % (int)m, (int)v.w % (int)m};
};

float vm::mix(const float v1, const float v2, const float m)
{
    return v1 * (1.f - m) + v2 * m;
};
vm::vec2 vm::mix(const vm::vec2 v1, const vm::vec2 v2, const float m)
{
    return vm::vec2{v1.x * (1.f - m) + v2.x * m, v1.y * (1.f - m) + v2.y * m};
};
vm::vec3 vm::mix(const vm::vec3 v1, const vm::vec3 v2, const float m)
{
    return vm::vec3{v1.x * (1.f - m) + v2.x * m, v1.y * (1.f - m) + v2.y * m, v1.y * (1.f - m) + v2.y * m};
};
vm::vec4 vm::mix(const vm::vec4 v1, const vm::vec4 v2, const float m)
{
    return vm::vec4{v1.x * (1.f - m) + v2.x * m, v1.y * (1.f - m) + v2.y * m, v1.y * (1.f - m) + v2.y * m, v1.w * (1.f - m) + v2.w * m};
};
vm::vec2 vm::mix(const vm::ivec2 v1, const vm::ivec2 v2, const float m)
{
    return vm::vec2{v1.x * (1.f - m) + v2.x * m, v1.y * (1.f - m) + v2.y * m};
};
vm::vec3 vm::mix(const vm::ivec3 v1, const vm::ivec3 v2, const float m)
{
    return vm::vec3{v1.x * (1.f - m) + v2.x * m, v1.y * (1.f - m) + v2.y * m, v1.y * (1.f - m) + v2.y * m};
};
vm::vec4 vm::mix(const vm::ivec4 v1, const vm::ivec4 v2, const float m)
{
    return vm::vec4{v1.x * (1.f - m) + v2.x * m, v1.y * (1.f - m) + v2.y * m, v1.y * (1.f - m) + v2.y * m, v1.w * (1.f - m) + v2.w * m};
};

vm::vec2 vm::floor(const vm::vec2 v)
{
    return vm::vec2{std::floor(v.x), std::floor(v.y)};
};
vm::vec3 vm::floor(const vm::vec3 v)
{
    return vm::vec3{std::floor(v.x), std::floor(v.y), std::floor(v.z)};
};
vm::vec4 vm::floor(const vm::vec4 v)
{
    return vm::vec4{std::floor(v.x), std::floor(v.y), std::floor(v.z), std::floor(v.w)};
};

float vm::fract(const float v)
{
    return v - std::floor(v);
};
vm::vec2 vm::fract(const vm::vec2 v)
{
    return vm::vec2{v.x - std::floor(v.x), v.y - std::floor(v.y)};
};
vm::vec3 vm::fract(const vm::vec3 v)
{
    return vm::vec3{v.x - std::floor(v.x), v.y - std::floor(v.y), v.z - std::floor(v.z)};
};
vm::vec4 vm::fract(const vm::vec4 v)
{
    return vm::vec4{v.x - std::floor(v.x), v.y - std::floor(v.y), v.z - std::floor(v.z), v.w - std::floor(v.w)};
};

// + operator
vm::vec2 vm::operator+(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return vm::vec2{v1.x + v2.x, v1.y + v2.y};
}
vm::vec3 vm::operator+(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return vm::vec3{v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}
vm::vec4 vm::operator+(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return vm::vec4{v1.x + v2.x,
                    v1.y + v2.y,
                    v1.z + v2.z,
                    v1.w + v2.w};
}
vm::ivec2 vm::operator+(const vm::ivec2 &v1, const vm::ivec2 &v2)
{
    return vm::ivec2{v1.x + v2.x, v1.y + v2.y};
}
vm::ivec3 vm::operator+(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}
vm::ivec4 vm::operator+(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w};
}

vm::vec2 vm::operator+(const vm::vec2 &v1, const float _m)
{
    return vm::vec2{v1.x + _m, v1.y + _m};
}
vm::vec3 vm::operator+(const vm::vec3 &v1, const float _m)
{
    return vm::vec3{v1.x + _m, v1.y + _m, v1.z + _m};
}
vm::vec4 vm::operator+(const vm::vec4 &v1, const float _m)
{
    return vm::vec4{v1.x + _m,
                    v1.y + _m,
                    v1.z + _m,
                    v1.w + _m};
}
vm::ivec2 vm::operator+(const vm::ivec2 &v1, const int _m)
{
    return vm::ivec2{v1.x + _m, v1.y + _m};
}
vm::ivec3 vm::operator+(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{v1.x + _m, v1.y + _m, v1.z + _m};
}
vm::ivec4 vm::operator+(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{v1.x + _m, v1.y + _m, v1.z + _m, v1.w + _m};
}
// - operator
vm::vec2 vm::operator-(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return vm::vec2{v1.x - v2.x, v1.y - v2.y};
}
vm::vec3 vm::operator-(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return vm::vec3{v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}
vm::vec4 vm::operator-(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return vm::vec4{v1.x - v2.x,
                    v1.y - v2.y,
                    v1.z - v2.z,
                    v1.w - v2.w};
}
vm::ivec2 vm::operator-(const vm::ivec2 &v1, const vm::ivec2 &v2)
{

    return vm::ivec2{v1.x - v2.x,
                     v1.y - v2.y};
}
vm::ivec3 vm::operator-(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{v1.x - v2.x,
                     v1.y - v2.y,
                     v1.z - v2.z};
}
vm::ivec4 vm::operator-(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{v1.x - v2.x,
                     v1.y - v2.y,
                     v1.z - v2.z,
                     v1.w - v2.w};
}

vm::vec2 vm::operator-(const vm::vec2 &v1, const float _m)
{
    return vm::vec2{v1.x - _m, v1.y - _m};
}
vm::vec3 vm::operator-(const vm::vec3 &v1, const float _m)
{
    return vm::vec3{v1.x - _m, v1.y - _m, v1.z - _m};
}
vm::vec4 vm::operator-(const vm::vec4 &v1, const float _m)
{
    return vm::vec4{v1.x - _m,
                    v1.y - _m,
                    v1.z - _m,
                    v1.w - _m};
}
vm::ivec2 vm::operator-(const vm::ivec2 &v1, const int _m)
{

    return vm::ivec2{v1.x - _m,
                     v1.y - _m};
}
vm::ivec3 vm::operator-(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{v1.x - _m,
                     v1.y - _m,
                     v1.z - _m};
}
vm::ivec4 vm::operator-(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{v1.x - _m,
                     v1.y - _m,
                     v1.z - _m,
                     v1.w - _m};
}
// * operator
vm::vec2 vm::operator*(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return vm::vec2{v1.x * v2.x, v1.y * v2.y};
}
vm::vec3 vm::operator*(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return vm::vec3{v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
}
vm::vec4 vm::operator*(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return vm::vec4{v1.x * v2.x,
                    v1.y * v2.y,
                    v1.z * v2.z,
                    v1.w * v2.w};
}
vm::ivec2 vm::operator*(const vm::ivec2 &v1, const vm::ivec2 &v2)
{

    return vm::ivec2{v1.x * v2.x,
                     v1.y * v2.y};
}
vm::ivec3 vm::operator*(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{v1.x * v2.x,
                     v1.y * v2.y,
                     v1.z * v2.z};
}
vm::ivec4 vm::operator*(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{v1.x * v2.x,
                     v1.y * v2.y,
                     v1.z * v2.z,
                     v1.w * v2.w};
}

vm::vec2 vm::operator*(const vm::vec2 &v1, const float _m)
{
    return vm::vec2{v1.x * _m, v1.y * _m};
}
vm::vec3 vm::operator*(const vm::vec3 &v1, const float _m)
{
    return vm::vec3{v1.x * _m, v1.y * _m, v1.z * _m};
}
vm::vec4 vm::operator*(const vm::vec4 &v1, const float _m)
{
    return vm::vec4{v1.x * _m,
                    v1.y * _m,
                    v1.z * _m,
                    v1.w * _m};
}
vm::ivec2 vm::operator*(const vm::ivec2 &v1, const int _m)
{

    return vm::ivec2{v1.x * _m,
                     v1.y * _m};
}
vm::ivec3 vm::operator*(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{v1.x * _m,
                     v1.y * _m,
                     v1.z * _m};
}
vm::ivec4 vm::operator*(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{v1.x * _m,
                     v1.y * _m,
                     v1.z * _m,
                     v1.w * _m};
}
// / operator
vm::vec2 vm::operator/(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return vm::vec2{v1.x / v2.x, v1.y / v2.y};
}
vm::vec3 vm::operator/(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return vm::vec3{v1.x / v2.x, v1.y / v2.y, v1.z / v2.z};
}
vm::vec4 vm::operator/(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return vm::vec4{v1.x / v2.x,
                    v1.y / v2.y,
                    v1.z / v2.z,
                    v1.w / v2.w};
}
vm::ivec2 vm::operator/(const vm::ivec2 &v1, const vm::ivec2 &v2)
{

    return vm::ivec2{v1.x / v2.x,
                     v1.y / v2.y};
}
vm::ivec3 vm::operator/(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{v1.x / v2.x,
                     v1.y / v2.y,
                     v1.z / v2.z};
}
vm::ivec4 vm::operator/(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{v1.x / v2.x,
                     v1.y / v2.y,
                     v1.z / v2.z,
                     v1.w / v2.w};
}

vm::vec2 vm::operator/(const vm::vec2 &v1, const float _m)
{
    return vm::vec2{v1.x / _m, v1.y / _m};
}
vm::vec3 vm::operator/(const vm::vec3 &v1, const float _m)
{
    return vm::vec3{v1.x / _m, v1.y / _m, v1.z / _m};
}
vm::vec4 vm::operator/(const vm::vec4 &v1, const float _m)
{
    return vm::vec4{v1.x / _m,
                    v1.y / _m,
                    v1.z / _m,
                    v1.w / _m};
}
vm::ivec2 vm::operator/(const vm::ivec2 &v1, const int _m)
{
    return vm::ivec2{v1.x / _m,
                     v1.y / _m};
}
vm::ivec3 vm::operator/(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{v1.x / _m,
                     v1.y / _m,
                     v1.z / _m};
}
vm::ivec4 vm::operator/(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{v1.x / _m,
                     v1.y / _m,
                     v1.z / _m,
                     v1.w / _m};
}
// += operator
vm::vec2 &vm::operator+=(vm::vec2 &v1, const vm::vec2 &v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    return v1;
}
vm::vec3 &vm::operator+=(vm::vec3 &v1, const vm::vec3 &v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}
vm::vec4 &vm::operator+=(vm::vec4 &v1, const vm::vec4 &v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    v1.w += v2.w;
    return v1;
}
vm::ivec2 &vm::operator+=(vm::ivec2 &v1, const vm::ivec2 &v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    return v1;
}
vm::ivec3 &vm::operator+=(vm::ivec3 &v1, const vm::ivec3 &v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

vm::ivec4 &vm::operator+=(vm::ivec4 &v1, const float _m)
{
    v1.x += _m;
    v1.y += _m;
    v1.z += _m;
    v1.w += _m;
    return v1;
}
vm::vec2 &vm::operator+=(vm::vec2 &v1, const float _m)
{
    v1.x += _m;
    v1.y += _m;
    return v1;
}
vm::vec3 &vm::operator+=(vm::vec3 &v1, const float _m)
{
    v1.x += _m;
    v1.y += _m;
    v1.z += _m;
    return v1;
}
vm::vec4 &vm::operator+=(vm::vec4 &v1, const float _m)
{
    v1.x += _m;
    v1.y += _m;
    v1.z += _m;
    v1.w += _m;
    return v1;
}
vm::ivec2 &vm::operator+=(vm::ivec2 &v1, const int _m)
{
    v1.x += _m;
    v1.y += _m;
    return v1;
}
vm::ivec3 &vm::operator+=(vm::ivec3 &v1, const int _m)
{
    v1.x += _m;
    v1.y += _m;
    v1.z += _m;
    return v1;
}
vm::ivec4 &vm::operator+=(vm::ivec4 &v1, const int _m)
{
    v1.x += _m;
    v1.y += _m;
    v1.z += _m;
    v1.w += _m;
    return v1;
}
// -= operator
vm::vec2 &vm::operator-=(vm::vec2 &v1, const vm::vec2 &v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    return v1;
}
vm::vec3 &vm::operator-=(vm::vec3 &v1, const vm::vec3 &v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}
vm::vec4 &vm::operator-=(vm::vec4 &v1, const vm::vec4 &v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    v1.w -= v2.w;
    return v1;
}
vm::ivec2 &vm::operator-=(vm::ivec2 &v1, const vm::ivec2 &v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    return v1;
}
vm::ivec3 &vm::operator-=(vm::ivec3 &v1, const vm::ivec3 &v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}

vm::ivec4 &vm::operator-=(vm::ivec4 &v1, const float _m)
{
    v1.x -= _m;
    v1.y -= _m;
    v1.z -= _m;
    v1.w -= _m;
    return v1;
}
vm::vec2 &vm::operator-=(vm::vec2 &v1, const float _m)
{
    v1.x -= _m;
    v1.y -= _m;
    return v1;
}
vm::vec3 &vm::operator-=(vm::vec3 &v1, const float _m)
{
    v1.x -= _m;
    v1.y -= _m;
    v1.z -= _m;
    return v1;
}
vm::vec4 &vm::operator-=(vm::vec4 &v1, const float _m)
{
    v1.x -= _m;
    v1.y -= _m;
    v1.z -= _m;
    v1.w -= _m;
    return v1;
}
vm::ivec2 &vm::operator-=(vm::ivec2 &v1, const int _m)
{
    v1.x -= _m;
    v1.y -= _m;
    return v1;
}
vm::ivec3 &vm::operator-=(vm::ivec3 &v1, const int _m)
{
    v1.x -= _m;
    v1.y -= _m;
    v1.z -= _m;
    return v1;
}
vm::ivec4 &vm::operator-=(vm::ivec4 &v1, const int _m)
{
    v1.x -= _m;
    v1.y -= _m;
    v1.z -= _m;
    v1.w -= _m;
    return v1;
}
// *= operator
vm::vec2 &vm::operator*=(vm::vec2 &v1, const vm::vec2 &v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    return v1;
}
vm::vec3 &vm::operator*=(vm::vec3 &v1, const vm::vec3 &v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}
vm::vec4 &vm::operator*=(vm::vec4 &v1, const vm::vec4 &v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    v1.w *= v2.w;
    return v1;
}
vm::ivec2 &vm::operator*=(vm::ivec2 &v1, const vm::ivec2 &v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    return v1;
}
vm::ivec3 &vm::operator*=(vm::ivec3 &v1, const vm::ivec3 &v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}

vm::ivec4 &vm::operator*=(vm::ivec4 &v1, const float _m)
{
    v1.x *= _m;
    v1.y *= _m;
    v1.z *= _m;
    v1.w *= _m;
    return v1;
}
vm::vec2 &vm::operator*=(vm::vec2 &v1, const float _m)
{
    v1.x *= _m;
    v1.y *= _m;
    return v1;
}
vm::vec3 &vm::operator*=(vm::vec3 &v1, const float _m)
{
    v1.x *= _m;
    v1.y *= _m;
    v1.z *= _m;
    return v1;
}
vm::vec4 &vm::operator*=(vm::vec4 &v1, const float _m)
{
    v1.x *= _m;
    v1.y *= _m;
    v1.z *= _m;
    v1.w *= _m;
    return v1;
}
vm::ivec2 &vm::operator*=(vm::ivec2 &v1, const int _m)
{
    v1.x *= _m;
    v1.y *= _m;
    return v1;
}
vm::ivec3 &vm::operator*=(vm::ivec3 &v1, const int _m)
{
    v1.x *= _m;
    v1.y *= _m;
    v1.z *= _m;
    return v1;
}
vm::ivec4 &vm::operator*=(vm::ivec4 &v1, const int _m)
{
    v1.x *= _m;
    v1.y *= _m;
    v1.z *= _m;
    v1.w *= _m;
    return v1;
}
// /= operator
vm::vec2 &vm::operator/=(vm::vec2 &v1, const vm::vec2 &v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    return v1;
}
vm::vec3 &vm::operator/=(vm::vec3 &v1, const vm::vec3 &v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}
vm::vec4 &vm::operator/=(vm::vec4 &v1, const vm::vec4 &v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    v1.w /= v2.w;
    return v1;
}
vm::ivec2 &vm::operator/=(vm::ivec2 &v1, const vm::ivec2 &v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    return v1;
}
vm::ivec3 &vm::operator/=(vm::ivec3 &v1, const vm::ivec3 &v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}

vm::ivec4 &vm::operator/=(vm::ivec4 &v1, const float _m)
{
    v1.x /= _m;
    v1.y /= _m;
    v1.z /= _m;
    v1.w /= _m;
    return v1;
}
vm::vec2 &vm::operator/=(vm::vec2 &v1, const float _m)
{
    v1.x /= _m;
    v1.y /= _m;
    return v1;
}
vm::vec3 &vm::operator/=(vm::vec3 &v1, const float _m)
{
    v1.x /= _m;
    v1.y /= _m;
    v1.z /= _m;
    return v1;
}
vm::vec4 &vm::operator/=(vm::vec4 &v1, const float _m)
{
    v1.x /= _m;
    v1.y /= _m;
    v1.z /= _m;
    v1.w /= _m;
    return v1;
}
vm::ivec2 &vm::operator/=(vm::ivec2 &v1, const int _m)
{
    v1.x /= _m;
    v1.y /= _m;
    return v1;
}
vm::ivec3 &vm::operator/=(vm::ivec3 &v1, const int _m)
{
    v1.x /= _m;
    v1.y /= _m;
    v1.z /= _m;
    return v1;
}
vm::ivec4 &vm::operator/=(vm::ivec4 &v1, const int _m)
{
    v1.x /= _m;
    v1.y /= _m;
    v1.z /= _m;
    v1.w /= _m;
    return v1;
}

// ^= operator
vm::vec2 &vm::operator^=(vm::vec2 &v1, const vm::vec2 &v2)
{
    v1.x = v1.x * v2.x - v1.y * v2.y;
    v1.y = v1.x * v2.y + v1.y * v2.x;
    return v1;
}
vm::vec3 &vm::operator^=(vm::vec3 &v1, const vm::vec3 &v2)
{
    v1.x = v1.x * v2.x - v1.y * v2.y - v1.z * v2.z;
    v1.y = v1.x * v2.y + v1.y * v2.x + v1.z * v2.z;
    v1.z = v1.x * v2.z + v1.y * v2.z + v1.z * v2.x;
    return v1;
}
vm::vec4 &vm::operator^=(vm::vec4 &v1, const vm::vec4 &v2)
{
    v1.x = v1.x * v2.x - v1.y * v2.y - v1.z * v2.z - v1.w * v2.w;
    v1.y = v1.x * v2.y + v1.y * v2.x + v1.z * v2.z + v1.w * v2.w;
    v1.z = v1.x * v2.z + v1.y * v2.z + v1.z * v2.x + v1.w * v2.w;
    v1.w = v1.x * v2.w + v1.y * v2.w + v1.z * v2.w + v1.w * v2.x;
    return v1;
}
vm::ivec2 &vm::operator^=(vm::ivec2 &v1, const vm::ivec2 &v2)
{
    v1.x = v1.x * v2.x - v1.y * v2.y;
    v1.y = v1.x * v2.y + v1.y * v2.x;
    return v1;
}

vm::ivec3 &vm::operator^=(vm::ivec3 &v1, const vm::ivec3 &v2)
{
    v1.x = v1.x * v2.x - v1.y * v2.y - v1.z * v2.z;
    v1.y = v1.x * v2.y + v1.y * v2.x + v1.z * v2.z;
    v1.z = v1.x * v2.z + v1.y * v2.z + v1.z * v2.x;
    return v1;
}
vm::ivec4 &vm::operator^=(vm::ivec4 &v1, const vm::ivec4 &v2)
{
    v1.x = v1.x * v2.x - v1.y * v2.y - v1.z * v2.z - v1.w * v2.w;
    v1.y = v1.x * v2.y + v1.y * v2.x + v1.z * v2.z + v1.w * v2.w;
    v1.z = v1.x * v2.z + v1.y * v2.z + v1.z * v2.x + v1.w * v2.w;
    v1.w = v1.x * v2.w + v1.y * v2.w + v1.z * v2.w + v1.w * v2.x;
    return v1;
}

vm::vec2 &vm::operator^=(vm::vec2 &v1, const float _m)
{
    v1.x = v1.x * _m - v1.y * _m;
    v1.y = v1.x * _m + v1.y * _m;
    return v1;
}
vm::vec3 &vm::operator^=(vm::vec3 &v1, const float _m)
{
    v1.x = v1.x * _m - v1.y * _m - v1.z * _m;
    v1.y = v1.x * _m + v1.y * _m + v1.z * _m;
    v1.z = v1.x * _m + v1.y * _m + v1.z * _m;
    return v1;
}
vm::vec4 &vm::operator^=(vm::vec4 &v1, const float _m)
{
    v1.x = v1.x * _m - v1.y * _m - v1.z * _m - v1.w * _m;
    v1.y = v1.x * _m + v1.y * _m + v1.z * _m + v1.w * _m;
    v1.z = v1.x * _m + v1.y * _m + v1.z * _m + v1.w * _m;
    v1.w = v1.x * _m + v1.y * _m + v1.z * _m + v1.w * _m;
    return v1;
}
vm::ivec2 &vm::operator^=(vm::ivec2 &v1, const int _m)
{
    v1.x = v1.x * _m - v1.y * _m;
    v1.y = v1.x * _m + v1.y * _m;
    return v1;
}
vm::ivec3 &vm::operator^=(vm::ivec3 &v1, const int _m)
{
    v1.x = v1.x * _m - v1.y * _m - v1.z * _m;
    v1.y = v1.x * _m + v1.y * _m + v1.z * _m;
    v1.z = v1.x * _m + v1.y * _m + v1.z * _m;
    return v1;
}
vm::ivec4 &vm::operator^=(vm::ivec4 &v1, const int _m)
{
    v1.x = v1.x * _m - v1.y * _m - v1.z * _m - v1.w * _m;
    v1.y = v1.x * _m + v1.y * _m + v1.z * _m + v1.w * _m;
    v1.z = v1.x * _m + v1.y * _m + v1.z * _m + v1.w * _m;
    v1.w = v1.x * _m + v1.y * _m + v1.z * _m + v1.w * _m;
    return v1;
}

// == operator
bool vm::operator==(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y);
}
bool vm::operator==(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}
bool vm::operator==(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z) && (v1.w == v2.w);
}
bool vm::operator==(const vm::ivec2 &v1, const vm::ivec2 &v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y);
}
bool vm::operator==(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}
bool vm::operator==(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z) && (v1.w == v2.w);
}

// != operator
bool vm::operator!=(const vm::vec2 &v1, const vm::vec2 &v2)
{
    return !(v1 == v2);
}
bool vm::operator!=(const vm::vec3 &v1, const vm::vec3 &v2)
{
    return !(v1 == v2);
}
bool vm::operator!=(const vm::vec4 &v1, const vm::vec4 &v2)
{
    return !(v1 == v2);
}
bool vm::operator!=(const vm::ivec2 &v1, const vm::ivec2 &v2)
{
    return !(v1 == v2);
}
bool vm::operator!=(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return !(v1 == v2);
}
bool vm::operator!=(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return !(v1 == v2);
}

// % operator
vm::ivec2 vm::operator%(const vm::ivec2 &v1, const vm::ivec2 &v2)
{
    return vm::ivec2{v1.x % v2.x,
                     v1.y % v2.y};
}
vm::ivec3 vm::operator%(const vm::ivec3 &v1, const vm::ivec3 &v2)
{
    return vm::ivec3{v1.x % v2.x,
                     v1.y % v2.y,
                     v1.z % v2.z};
}
vm::ivec4 vm::operator%(const vm::ivec4 &v1, const vm::ivec4 &v2)
{
    return vm::ivec4{v1.x % v2.x,
                     v1.y % v2.y,
                     v1.z % v2.z,
                     v1.w % v2.w};
}
vm::ivec2 vm::operator%(const vm::ivec2 &v1, const int _m)
{
    return vm::ivec2{v1.x % _m,
                     v1.y % _m};
}
vm::ivec3 vm::operator%(const vm::ivec3 &v1, const int _m)
{
    return vm::ivec3{v1.x % _m,
                     v1.y % _m, v1.z % _m};
}
vm::ivec4 vm::operator%(const vm::ivec4 &v1, const int _m)
{
    return vm::ivec4{v1.x % _m,
                     v1.y % _m, v1.z % _m, v1.w % _m};
}
