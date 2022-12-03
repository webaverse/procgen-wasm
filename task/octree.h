#ifndef OCTREE_H
#define OCTREE_H

#define _USE_MATH_DEFINES

#include "../libs/vectorMath.h"
#include "../libs/vector.h"
#include "../constants.h"
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <stdint.h>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <set>
#include <emscripten.h>

#define M_PI 3.14159265358979323846

//


const vm::ivec3 chunkMinForPosition(const vm::ivec3 &p, const int &lod);
const vm::ivec3 chunkMinForPosition(const vm::vec3 &p, const int &lod);

uint64_t hashOctreeMin(const vm::ivec2 &min);
uint64_t hashOctreeMin(const vm::ivec3 &min);
uint64_t hashOctreeMinLod(const vm::ivec2 &min, int lod);
uint64_t hashOctreeMinLod(const vm::ivec3 &min, int lod);

class OctreeNodeSpec {
public:
    vm::ivec2 min;
    int lod;
    int lodArray[2];
};

class OctreeNode : public OctreeNodeSpec {
public:
    OctreeNode(const vm::ivec2 &min, int lod) :
      OctreeNodeSpec{min, lod, {0, 0}}
      {}
    int getPriority() const {
        return 0;
    }
    Box3 getBox() const {
        float height = (float)(MAX_WORLD_HEIGHT - (float)MIN_WORLD_HEIGHT);
        Vec boxMin{(float)min.x, (float)MIN_WORLD_HEIGHT, (float)min.y};
        return Box3(
            boxMin,
            boxMin + Vec{(float)lod, height, (float)lod}
        );
    }
};
typedef std::shared_ptr<OctreeNode> OctreeNodePtr;

//

#endif // OCTREE_H