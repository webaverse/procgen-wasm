#ifndef _TRACKER_H_
#define _TRACKER_H_

#include "octree.h"
#include "vectorMath.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <functional>
#include <emscripten/atomic.h>

//

class OctreeNode;

//

constexpr size_t numCachedOctreeNodes = 64 * 1024;
typedef std::shared_ptr<OctreeNode> OctreeNodePtr;

class OctreeNodeAllocator {
public: 
  static OctreeNodePtr alloc(const vm::ivec2 &min, int lod);
};

class OctreeContext {
public:
  std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>> nodeMap;

  OctreeNodePtr alloc(const vm::ivec2 &min, int lod);
};

//

class DataRequest {
public:
  OctreeNodePtr node;
  
  std::vector<uint8_t> getBuffer() const;
};

typedef std::shared_ptr<DataRequest> DataRequestPtr;

class DataRequestUpdate {
public:
  std::vector<DataRequestPtr> newDataRequests;
  std::vector<DataRequestPtr> keepDataRequests;
  std::vector<DataRequestPtr> cancelDataRequests;

  std::vector<uint8_t> getBuffer() const;
};

//

class Dominator {
public:
  OctreeNodePtr node;
  std::vector<OctreeNodePtr> newChunks;
  std::vector<OctreeNodePtr> oldChunks;

  Dominator(OctreeNodePtr node);

  std::vector<uint8_t> getBuffer() const;
};

//

extern std::atomic<int> nextTrackerId;

/* enum TrackerTaskType {
  ADD = 1,
  REMOVE = 2,
  OUTRANGE = 3
};

class TrackerTask {
public:
  int id;
  TrackerTaskType type;
  OctreeNodePtr maxLodNode;
  std::vector<OctreeNodePtr> oldNodes;
  std::vector<OctreeNodePtr> newNodes;

  bool isNop() const;

  std::vector<uint8_t> getBuffer() const;
};
typedef std::shared_ptr<TrackerTask> TrackerTaskPtr; */

class TrackerUpdate {
public:
  std::vector<OctreeNodePtr> leafNodes;

  std::vector<DataRequestPtr> newDataRequests;
  std::vector<DataRequestPtr> keepDataRequests;
  std::vector<DataRequestPtr> cancelDataRequests;

  // std::unordered_map<uint64_t, Dominator> dominators;

  uint8_t *getBuffer() const;
};

//

bool containsPoint(const OctreeNode &node, const vm::ivec3 &p);
bool containsNode(const OctreeNode &node, const OctreeNode &other);
bool equalsNode(const OctreeNode &node, const OctreeNode &other);
bool equalsNodeLod(const OctreeNode &node, const OctreeNode &other);
// bool intersectsNode(const OctreeNode &node, const OctreeNode &other);

//

/* extern std::array<vm::ivec3, 8> lodOffsets;
extern OctreeNodeAllocator octreeNodeAllocator;

//

enum ConstructTreeFlags {
  ConstructTreeFlags_None = 0,
  ConstructTreeFlags_Children = 1,
  ConstructTreeFlags_Peers = 2,
}; */

// OctreeNodePtr getLeafNodeFromPoint(const std::vector<OctreeNode *> &leafNodes, const vm::ivec2 &p);
// OctreeNodePtr getNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod);
// OctreeNodePtr createNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod, bool isLeaf);
// OctreeNodePtr getOrCreateNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod, bool isLeaf);
// void ensureChildren(OctreeContext &octreeContext, OctreeNode *parentNode);
// void ensurePeers(OctreeContext &octreeContext, OctreeNode *node, int maxLod);
// void constructTreeUpwards(OctreeContext &octreeContext, const vm::ivec2 &leafPosition, int minLod, int maxLod, int flags);

std::vector<OctreeNodePtr> constructOctreeForLeaf(const vm::ivec2 &position, int lod1Range, int maxLod);
OctreeNodePtr getMaxLodNode(const std::vector<OctreeNodePtr> &newLeafNodes, const std::vector<OctreeNodePtr> &oldLeafNodes, const vm::ivec3 &min);
// std::vector<TrackerTaskPtr> diffLeafNodes(const std::vector<OctreeNodePtr> &newLeafNodes, const std::vector<OctreeNodePtr> &oldLeafNodes);
// std::vector<TrackerTaskPtr> sortTasks(const std::vector<TrackerTaskPtr> &tasks, const vm::vec3 &worldPosition);
// std::pair<std::vector<OctreeNodePtr>, std::vector<TrackerTaskPtr>> updateChunks(const std::vector<OctreeNodePtr> &oldChunks, const std::vector<TrackerTaskPtr> &tasks);

//

class Tracker {
public:
  PGInstance *inst;
  int lods;
  int lod1Range;
  
  vm::ivec2 lastCoord;
  std::vector<OctreeNodePtr> leafNodes;
  std::unordered_map<uint64_t, DataRequestPtr> dataRequests;

  //

  Tracker(PGInstance *inst, int lods, int lod1Range);

  // static methods

  vm::ivec2 getCurrentCoord(const vm::vec3 &position);

  // dynamic methods

  void sortNodes(std::vector<OctreeNodePtr> &nodes);
  DataRequestUpdate updateDataRequests(
    std::unordered_map<uint64_t, DataRequestPtr> &dataRequests,
    const std::vector<OctreeNodePtr> &leafNodes
  );
  /* std::unordered_map<uint64_t, Dominator> updateDominators(
    const std::vector<OctreeNodePtr> &oldRenderedChunks,
    std::unordered_map<uint64_t, DataRequestPtr> &dataRequests
  ); */
  TrackerUpdate updateCoord(const vm::ivec2 &currentCoord);
  TrackerUpdate update(const vm::vec3 &position);
};

#endif // _TRACKER_H_