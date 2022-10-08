#ifndef _TRACKER_H_
#define _TRACKER_H_

#include "vectorMath.h"
#include "octree.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <functional>
#include <emscripten/atomic.h>

//

class PGInstance;

//

constexpr size_t numCachedOctreeNodes = 64 * 1024;

class OctreeNodeAllocator {
public: 
  static OctreeNodePtr alloc(const vm::ivec2 &min, int lod);
};

class OctreeContext {
public:
  std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>> nodeMap;

  OctreeNodePtr alloc(const vm::ivec2 &min, int lod);
  std::vector<OctreeNodePtr> getLeafNodes();
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
  std::unordered_map<uint64_t, DataRequestPtr> dataRequests;
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
  // vm::ivec2 chunkPosition;

  // std::unordered_map<uint64_t, Dominator> dominators;

  uint8_t *getBuffer() const;
};

//

bool containsPoint(const vm::ivec2 &min, const int lod, const vm::ivec2 &p);
bool containsPoint(const vm::vec2 &min, const int lod, const vm::vec2 &p);
bool containsPoint(const OctreeNode &node, const vm::ivec2 &p);
bool containsNode(const OctreeNode &node, const OctreeNode &other);
bool equalsNode(const OctreeNode &node, const OctreeNode &other);
bool equalsNodeLod(const OctreeNode &node, const OctreeNode &other);
// bool intersectsNode(const OctreeNode &node, const OctreeNode &other);
bool chunksIntersect(const vm::ivec2 &min1, const int lod1, const vm::ivec2 &min2, const int lod2);

void constructSeedTree(
  OctreeContext &octreeContext,
  const std::vector<vm::ivec2> &maxLodChunkPositions,
  const int maxLod,
  const std::vector<std::pair<vm::ivec2, int>> &lodSplits
);
// std::vector<OctreeNodePtr> constructOctreeForSeed(const vm::ivec2 &maxLodCenter, const int maxLod, const std::vector<std::pair<vm::ivec2, int>> &lodSplits);
std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>>::iterator findNodeIterAtPoint(OctreeContext &octreeContext, const vm::ivec2 &position);

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
  
  // vm::ivec2 lastCoord;
  std::unordered_map<uint64_t, DataRequestPtr> dataRequests;
  // int lastEpoch;
  std::mutex mutex;

  //

  Tracker(PGInstance *inst);

  // dynamic methods

  void sortNodes(std::vector<OctreeNodePtr> &nodes);
  DataRequestUpdate updateDataRequests(
    const std::unordered_map<uint64_t, DataRequestPtr> &dataRequests,
    const std::vector<OctreeNodePtr> &leafNodes
  );
  /* std::unordered_map<uint64_t, Dominator> updateDominators(
    const std::vector<OctreeNodePtr> &oldRenderedChunks,
    std::unordered_map<uint64_t, DataRequestPtr> &dataRequests
  ); */
  TrackerUpdate update(const vm::vec3 &position, int lods, int lod1Range);
};

#endif // _TRACKER_H_