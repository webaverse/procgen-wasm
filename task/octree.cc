#include "octree.h"
#include <cmath>

vm::ivec2 chunkMinForPosition(const vm::ivec2 &p, const int &lod)
{
    return vm::ivec2{
        (int)std::floor((float)p.x / (float)lod) * lod,
        (int)std::floor((float)p.y / (float)lod) * lod
    };
}

uint64_t hashOctreeMin(const vm::ivec2 &min)
{
    uint64_t result = uint16_t(min.x);
    result = (result << 16) | uint16_t(min.y);
    return result;
}
uint64_t hashOctreeMin(const vm::ivec3 &min)
{
    uint64_t result = uint16_t(min.x);
    result = (result << 16) | uint16_t(min.y);
    result = (result << 16) | uint16_t(min.z);
    return result;
}

uint64_t hashOctreeMinLod(const vm::ivec2 &min, int lod)
{
    uint64_t result = uint16_t(min.x);
    result = (result << 16) | uint16_t(min.y);
    result = (result << 16) | uint16_t(lod);
    return result;
}
uint64_t hashOctreeMinLod(const vm::ivec3 &min, int lod)
{
    uint64_t result = uint16_t(min.x);
    result = (result << 16) | uint16_t(min.y);
    result = (result << 16) | uint16_t(min.z);
    result = (result << 16) | uint16_t(lod);
    return result;
}
