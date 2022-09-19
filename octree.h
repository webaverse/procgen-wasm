#ifndef OCTREE_H
#define OCTREE_H

#include "vectorMath.h"
#define _USE_MATH_DEFINES

#include "procgen.h"
// #include "qef.h"
#include "mesh.h"
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <stdint.h>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <set>
// #include "density.h"
#include <emscripten.h>

#define M_PI 3.14159265358979323846

//


const vm::ivec3 chunkMinForPosition(const vm::ivec3 &p, const int &lod);
const vm::ivec3 chunkMinForPosition(const vm::vec3 &p, const int &lod);

uint64_t hashOctreeMin(const vm::ivec2 &min);
uint64_t hashOctreeMin(const vm::ivec3 &min);
uint64_t hashOctreeMinLod(const vm::ivec2 &min, int lod);
uint64_t hashOctreeMinLod(const vm::ivec3 &min, int lod);

//

/* enum OctreeNodeType
{
    Node_Internal,
    Node_Leaf
}; */

//

class OctreeNodeSpec {
public:
    vm::ivec2 min;
    int lod;
    int lodArray[2];
};

class OctreeNode;
typedef std::shared_ptr<OctreeNode> OctreeNodePtr;

class OctreeNode : public OctreeNodeSpec
{
public:
    OctreeNode(const vm::ivec2 &min, int lod) :
      OctreeNodeSpec{min, lod, {0, 0}}
      {}

    /* bool isLeaf() const {
        for (const auto &child : children) {
            if (child) {
                return false;
            }
        }
        return true;
    } */
    int getPriority() const {
        return 0;
    }
    Box3 getBox() const {
        /* float halfSize = (float)lod / 2.0f;
        float halfHeight = (float)(MIN_WORLD_HEIGHT + MAX_WORLD_HEIGHT) / 2.f;
        Vec center{
            (float)min.x + halfSize,
            halfHeight,
            (float)min.y + halfSize,
        };
        return Box3(
            center - Vec{halfSize, halfHeight, halfSize},
            center + Vec{halfSize, halfHeight, halfSize}
        ); */
        float height = (float)(MAX_WORLD_HEIGHT - (float)MIN_WORLD_HEIGHT);
        Vec boxMin{(float)min.x, (float)MIN_WORLD_HEIGHT, (float)min.y};
        return Box3(
            boxMin,
            boxMin + Vec{(float)lod, height, (float)lod}
        );
    }
    uint8_t *getBuffer() const {
        constexpr size_t size = sizeof(min) + sizeof(lod);
        uint8_t *buffer = (uint8_t *)malloc(size);
        int index = 0;
        memcpy(buffer + index, &min, sizeof(min));
        index += sizeof(min);
        memcpy(buffer + index, &lod, sizeof(lod));
        index += sizeof(lod);

        // std::cout << "output node " << min.x << " " << min.y << " : " << lod << std::endl;

        return buffer;
    }
};

//

template <typename DCContextType>
vm::vec3 calculateSurfaceNormal(const vm::vec3 &p, const int lod, DCInstance *inst)
{
    // finding the surface normal with the derivative
    constexpr float H = 0.001f;
    const float dx = DCContextType::densityFn(p + vm::vec3{H, 0.f, 0.f}, lod, inst) -
                     DCContextType::densityFn(p - vm::vec3{H, 0.f, 0.f}, lod, inst);
    const float dy = DCContextType::densityFn(p + vm::vec3{0.f, H, 0.f}, lod, inst) -
                     DCContextType::densityFn(p - vm::vec3{0.f, H, 0.f}, lod, inst);
    const float dz = DCContextType::densityFn(p + vm::vec3{0.f, 0.f, H}, lod, inst) -
                     DCContextType::densityFn(p - vm::vec3{0.f, 0.f, H}, lod, inst);
    return vm::normalize(vm::vec3{dx, dy, dz});
}

#endif // OCTREE_H