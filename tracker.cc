#include "tracker.h"
#include "octree.h"
#include "sort.h"
#include <iostream>

/*
tracker = await dcWorkerManager.createTracker(3, 2, true); trackerUpdate = await dcWorkerManager.trackerUpdate(tracker, new THREE.Vector3(1, 2, 3));
*/

//

// std::atomic<int> nextTrackerId(0);

//

/* bool TrackerTask::isNop() const {
  auto task = this;
  return task->newNodes.size() == task->oldNodes.size() &&
    std::all_of(task->newNodes.begin(), task->newNodes.end(), [&](OctreeNodePtr newNode) -> bool {
      return std::any_of(task->oldNodes.begin(), task->oldNodes.end(), [&](OctreeNodePtr oldNode) -> bool {
        return equalsNodeLod(*oldNode, *newNode);
      });
    });
}
std::vector<uint8_t> TrackerTask::getBuffer() const {
  size_t size = 0;
  // header
  size += sizeof(int); // id
  size += sizeof(int); // type
  // max lod node
  size += sizeof(vm::ivec3); // min
  size += sizeof(int); // lod
  size += sizeof(int); // isLeaf
  size += sizeof(int[8]); // lodArray
  // old nodes
  size += sizeof(uint32_t); // numOldNodes
  size += sizeof(vm::ivec3) * oldNodes.size(); // min
  size += sizeof(int) * oldNodes.size(); // lod
  size += sizeof(int) * oldNodes.size(); // isLeaf
  size += sizeof(int[8]) * oldNodes.size(); // lodArray
  // new nodes
  size += sizeof(uint32_t); // numNewNodes
  size += sizeof(vm::ivec3) * newNodes.size(); // min
  size += sizeof(int) * newNodes.size(); // lod
  size += sizeof(int) * newNodes.size(); // isLeaf
  size += sizeof(int[8]) * newNodes.size(); // lodArray

  std::vector<uint8_t> result(size);
  int index = 0;
  // id
  *((int *)(result.data() + index)) = id;
  index += sizeof(int);
  // type
  *((int *)(result.data() + index)) = type;
  index += sizeof(int);
  // max lod node
  std::memcpy(result.data() + index, &maxLodNode->min, sizeof(vm::ivec3));
  index += sizeof(vm::ivec3);
  *((int *)(result.data() + index)) = maxLodNode->size;
  index += sizeof(int);
  *((int *)(result.data() + index)) = (maxLodNode->type == Node_Leaf) ? 1 : 0;
  index += sizeof(int);
  std::memcpy(result.data() + index, &maxLodNode->lodArray[0], sizeof(int[8]));
  index += sizeof(int[8]);
  // old nodes
  *((uint32_t *)(result.data() + index)) = oldNodes.size();
  index += sizeof(uint32_t);
  for (auto oldNode : oldNodes) {
    std::memcpy(result.data() + index, &oldNode->min, sizeof(vm::ivec3));
    index += sizeof(vm::ivec3);
    *((int *)(result.data() + index)) = oldNode->size;
    index += sizeof(int);
    *((int *)(result.data() + index)) = (oldNode->type == Node_Leaf) ? 1 : 0;
    index += sizeof(int);
    std::memcpy(result.data() + index, &oldNode->lodArray[0], sizeof(int[8]));
    index += sizeof(int[8]);
  }
  // new nodes
  *((uint32_t *)(result.data() + index)) = newNodes.size();
  index += sizeof(uint32_t);
  for (auto newNode : newNodes) {
    std::memcpy(result.data() + index, &newNode->min, sizeof(vm::ivec3));
    index += sizeof(vm::ivec3);
    *((int *)(result.data() + index)) = newNode->size;
    index += sizeof(int);
    *((int *)(result.data() + index)) = (newNode->type == Node_Leaf) ? 1 : 0;
    index += sizeof(int);
    std::memcpy(result.data() + index, &newNode->lodArray[0], sizeof(int[8]));
    index += sizeof(int[8]);
  }

  return result;
} */
void serializeOctreeNodes(const std::vector<OctreeNodePtr> &datas, uint8_t *ptr, int &index) {
  // std::cout << "serialize octree nodes" << datas.size() << std::endl;
  
  *((int32_t *)(ptr + index)) = datas.size();  
  index += sizeof(int32_t);

  for (const auto &data : datas) {
    std::memcpy(ptr + index, &data->min, sizeof(vm::ivec2));
    index += sizeof(vm::ivec2);

    *((int *)(ptr + index)) = data->lod;
    index += sizeof(int);

    std::memcpy(ptr + index, data->lodArray, sizeof(int[2]));
    index += sizeof(int[2]);
  }
}
void serializeDataRequests(const std::vector<DataRequestPtr> &datas, uint8_t *ptr, int &index) {
  // std::cout << "serialize data requests" << datas.size() << std::endl;

  *((int32_t *)(ptr + index)) = datas.size();  
  index += sizeof(int32_t);

  for (const auto &data : datas) {
    std::memcpy(ptr + index, &data->node->min, sizeof(vm::ivec2));
    index += sizeof(vm::ivec2);

    *((int *)(ptr + index)) = data->node->lod;
    index += sizeof(int);
    
    std::memcpy(ptr + index, data->node->lodArray, sizeof(int[2]));
    index += sizeof(int[2]);
  }
}
/* void serializeChunkPosition(const vm::ivec2 &chunkPosition, uint8_t *ptr, int &index) {
  memcpy(ptr + index, &chunkPosition, sizeof(chunkPosition));
  index += sizeof(vm::ivec2);
} */

uint8_t *TrackerUpdate::getBuffer() const {
  // compute size
  size_t size = 0;

  constexpr size_t octreeNodeSize = sizeof(vm::ivec2) + sizeof(int) + sizeof(int[2]);

  size += sizeof(int32_t); // numLeafNodes
  size += octreeNodeSize * leafNodes.size();

  size += sizeof(int32_t); // numNewDataRequests
  size += octreeNodeSize * newDataRequests.size();

  size += sizeof(int32_t); // numKeepDataRequests
  size += octreeNodeSize * keepDataRequests.size();

  size += sizeof(int32_t); // numCancelDataRequests
  size += octreeNodeSize * cancelDataRequests.size();

  // size += sizeof(vm::ivec2); // chunkPosition

  // serialize
  uint8_t *ptr = (uint8_t *)malloc(size);
  int index = 0;

  serializeOctreeNodes(leafNodes, ptr, index);
  serializeDataRequests(newDataRequests, ptr, index);
  serializeDataRequests(keepDataRequests, ptr, index);
  serializeDataRequests(cancelDataRequests, ptr, index);
  // serializeChunkPosition(chunkPosition, ptr, index);

  return ptr;
}

//

/* std::vector<uint8_t> DataRequest::getBuffer() const {
  size_t size = 0;
  size += sizeof(vm::ivec3); // min
  size += sizeof(int); // lod

  std::vector<uint8_t> result(size);
  int index = 0;
  
  // min
  std::memcpy(result.data() + index, &node->min, sizeof(vm::ivec3));
  index += sizeof(vm::ivec3);
  
  // lod
  *((int *)(result.data() + index)) = node->lod;
  index += sizeof(int);
  
  return result;
}

std::vector<uint8_t> DataRequestUpdate::getBuffer() const {
  return std::vector<uint8_t>(); // XXX
} */

//

Dominator::Dominator(OctreeNodePtr node) :
  node(node)
  {}

//

std::vector<uint8_t> Dominator::getBuffer() const {
  size_t size = 0;
  size += sizeof(uint32_t); // numNewChunks
  size += sizeof(uint32_t); // numOldChunks
  size += sizeof(uint32_t) * newChunks.size(); // new chunk hashes
  size += sizeof(uint32_t) * oldChunks.size(); // old chunk hashes

  std::vector<uint8_t> result(size);
  int index = 0;
  // numNewChunks
  *((uint32_t *)(result.data() + index)) = newChunks.size();
  index += sizeof(uint32_t);
  // numOldChunks
  *((uint32_t *)(result.data() + index)) = oldChunks.size();
  index += sizeof(uint32_t);
  // new chunk hashes
  for (auto chunk : newChunks) {
    uint64_t hash = hashOctreeMinLod(chunk->min, chunk->lod);
    *((uint32_t *)(result.data() + index)) = hash;
    index += sizeof(uint32_t);
  }
  // old chunk hashes
  for (auto chunk : oldChunks) {
    uint64_t hash = hashOctreeMinLod(chunk->min, chunk->lod);
    *((uint32_t *)(result.data() + index)) = hash;
    index += sizeof(uint32_t);
  }
  
  return result;
}

//

bool containsPoint(const vm::ivec2 &min, const int lod, const vm::ivec2 &p) {
  return p.x >= min.x && p.x < min.x + lod &&
    p.y >= min.y && p.y < min.y + lod;
}
bool containsPoint(const vm::vec2 &min, const int lod, const vm::vec2 &p) {
  return p.x >= min.x && p.x < min.x + (float)lod &&
    p.y >= min.y && p.y < min.y + (float)lod;
}
bool containsPoint(const OctreeNode &node, const vm::ivec2 &p) {
  return containsPoint(node.min, node.lod, p);
}
bool containsNode(const OctreeNode &node, const OctreeNode &other) {
  return containsPoint(node, other.min);
}
bool equalsNode(const OctreeNode &node, const OctreeNode &other) {
  if (node.min == other.min && node.lod == other.lod) {
    return true;
  } else {
    return false;
  }
}
bool equalsNodeLod(const OctreeNode &node, const OctreeNode &other) {
  if (equalsNode(node, other)) {
    if (node.lod != other.lod) {
      return false;
    }
    /* for (size_t i = 0; i < node.lodArray.size(); i++) {
      if (node.lodArray[i] != other.lodArray[i]) {
        return false;
      }
    } */
    return true;
  } else {
    return false;
  }
}
bool equalsNodeLodArray(const OctreeNode &node, const OctreeNode &other) {
  if (equalsNode(node, other)) {
    if (
      (node.lod != other.lod) ||
      (node.lodArray[0] != other.lodArray[0] || node.lodArray[1] != other.lodArray[1])
    ) {
      return false;
    }
    /* for (size_t i = 0; i < node.lodArray.size(); i++) {
      if (node.lodArray[i] != other.lodArray[i]) {
        return false;
      }
    } */
    return true;
  } else {
    return false;
  }
}
/* bool intersectsNode(const OctreeNode &node, const OctreeNode &other) {
  return containsNode(node, other) || containsNode(other, node);
} */
// returns whether the box of min1 + lod1 intersects the box of min2 + lod2
bool chunksIntersect(const vm::ivec2 &min1, const int lod1, const vm::ivec2 &min2, const int lod2) {
  return min1.x < min2.x + lod2 && min1.x + lod1 > min2.x &&
    min1.y < min2.y + lod2 && min1.y + lod1 > min2.y;
}

//

/* std::array<vm::ivec3, 8> lodOffsets = {
    vm::ivec3{0, 0, 0},
    vm::ivec3{1, 0, 0},
    vm::ivec3{0, 0, 1},
    vm::ivec3{1, 0, 1},
    vm::ivec3{0, 1, 0},
    vm::ivec3{1, 1, 0},
    vm::ivec3{0, 1, 1},
    vm::ivec3{1, 1, 1}
}; */
// OctreeNodeAllocator octreeNodeAllocator;

//

/* OctreeNodeAllocator::OctreeNodeAllocator() {
  for (int i = 0; i < rawNodes.size(); i++) {
    OctreeNode *node = &rawNodes[i];
    nodes.push_back(node);
  }
  this->deleter = [this](void *ptr) -> void {
    this->nodes.push_back((OctreeNode *)ptr);
  };
} */
OctreeNodePtr OctreeNodeAllocator::alloc(const vm::ivec2 &min, int lod) {
  OctreeNode *octreeNode = new OctreeNode(min, lod);
  // octreeNode->min = min;
  // octreeNode->lod = lod;
  return std::shared_ptr<OctreeNode>(octreeNode);
}

OctreeNodePtr OctreeContext::alloc(const vm::ivec2 &min, int lod) {
  auto node = OctreeNodeAllocator::alloc(min, lod);
  if (node) {
    node->min = min;
    node->lod = lod;
    
    uint64_t hash = hashOctreeMinLod(min, lod);
    nodeMap[hash] = node;
  }
  return node;
}
// OctreeContext::OctreeContext() {}
std::vector<OctreeNodePtr> OctreeContext::getLeafNodes() {
  std::vector<OctreeNodePtr> leafNodes;
  for (const auto &iter : nodeMap) {
    auto node = iter.second;
    leafNodes.push_back(node);
  }
  return leafNodes;
}

//

/* OctreeNodePtr getLeafNodeFromPoint(const std::vector<OctreeNodePtr> &leafNodes, const vm::ivec2 &p) {
  for (size_t i = 0; i < leafNodes.size(); i++) {
      auto leafNode = leafNodes[i];
      if (containsPoint(*leafNode, p)) {
          return leafNode;
      }
  }
  return nullptr;
} */
std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>>::iterator getNodeIter(OctreeContext &octreeContext, const vm::ivec2 &min, int lod) {
  auto &nodeMap = octreeContext.nodeMap;

  uint64_t hash = hashOctreeMinLod(min, lod);
  auto iter = nodeMap.find(hash);
  return iter;
}
OctreeNodePtr getNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod) {
  auto &nodeMap = octreeContext.nodeMap;

  auto iter = getNodeIter(octreeContext, min, lod);
  if (iter != nodeMap.end()) {
    return iter->second;
  } else {
    return nullptr;
  }
}
OctreeNodePtr createNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod) {
  auto &nodeMap = octreeContext.nodeMap;

  auto node = OctreeNodeAllocator::alloc(min, lod);
  uint64_t hash = hashOctreeMinLod(min, lod);
  nodeMap[hash] = node;
  return node;
}
OctreeNodePtr getOrCreateNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod) {
  OctreeNodePtr node = getNode(octreeContext, min, lod);
  if (!node) {
      node = createNode(octreeContext, min, lod);
  }
  return node;
}
/* OctreeNodePtr scanForNode(OctreeContext &octreeContext, const vm::ivec2 &min, int minSearchLod, int maxSearchLod) {  
  for (int lod = minSearchLod; lod <= maxSearchLod; lod *= 2) {
    vm::ivec2 snapPosition{
      (min.x / lod) * lod,
      (min.y / lod) * lod
    };
    auto node = getNode(octreeContext, snapPosition, lod);
    if (node) {
      return node;
    }
  }
  return nullptr;
} */
/* OctreeNodePtr findContainingNode(OctreeContext &octreeContext, const vm::ivec2 &min) {
  auto &nodeMap = octreeContext.nodeMap;

  std::vector<OctreeNodePtr> leafNodes;
  for (const auto &iter : nodeMap) {
    auto node = iter.second;
    if (containsPoint(*node, min)) {
      return node;
    }
  }
  return nullptr;
} */
std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>>::iterator findNodeIterAtPoint(OctreeContext &octreeContext, const vm::ivec2 &position) {
  auto &nodeMap = octreeContext.nodeMap;

  // std::vector<OctreeNodePtr> leafNodes;
  for (auto iter = nodeMap.begin(); iter != nodeMap.end(); iter++) {
    auto node = iter->second;
    if (containsPoint(node->min, node->lod, position)) {
      return iter;
    }
  }
  return nodeMap.end();
}
OctreeNodePtr findContainedNode(OctreeContext &octreeContext, const vm::ivec2 &min, const int lod) {
  auto &nodeMap = octreeContext.nodeMap;

  // std::vector<OctreeNodePtr> leafNodes;
  for (const auto &iter : nodeMap) {
    auto node = iter.second;
    if (containsPoint(min, lod, node->min)) {
      return node;
    }
  }
  return nullptr;
}
bool removeNode(OctreeContext &octreeContext, OctreeNodePtr node) {
  auto &nodeMap = octreeContext.nodeMap;

  uint64_t hash = hashOctreeMinLod(node->min, node->lod);
  auto iter = nodeMap.find(hash);
  if (iter != nodeMap.end()) {
    nodeMap.erase(iter);
    return true;
  } else {
    std::cerr << "erase node not found: " << node->min.x << " " << node->min.y << " " << node->lod << std::endl;
    abort();
    return false;
  }
}
/* void ensureChildren(OctreeContext &octreeContext, OctreeNode *parentNode) {
    // auto &nodeMap = octreeContext.nodeMap;

    const vm::ivec3 &lodMin = parentNode->min;
    const int &lod = parentNode->size;

    for (int dx = 0; dx < 2; dx++) {
      for (int dy = 0; dy < 2; dy++) {
        for (int dz = 0; dz < 2; dz++) {
           int childIndex = dx + 2 * (dy + 2 * dz);
           if (parentNode->children[childIndex] == nullptr) {
              parentNode->children[childIndex] = createNode(
                octreeContext,
                lodMin + vm::ivec3{dx, dy, dz} * (lod / 2),
                lod / 2,
                true
              ).get();
              // parentNode->children[childIndex].parent = parentNode;
           }
        }
      }
    }
} */
/* void ensureChildLayer(OctreeContext &octreeContext, OctreeNodePtr parentNode) {
  const vm::ivec2 &parentPosition = parentNode->min;
  int childLod = parentNode->lod / 2;
  
  for (int i = 0; i < parentNode->children.size(); i++) {
    // if (!parentNode->children[i]) {
      vm::ivec2 childPosition{
        parentPosition.x + (i % 2) * childLod,
        parentPosition.y + (i / 2) * childLod
      };
      OctreeNodePtr childNode = getOrCreateNode(
        octreeContext,
        childPosition,
        childLod
      );
      parentNode->children[i] = childNode;
    // }
  }
} */
void splitPointToLod(
  OctreeContext &octreeContext,
  const vm::ivec2 &absolutePosition,
  int targetLod
) {
  auto &nodeMap = octreeContext.nodeMap;

  for (;;) {
    auto iter = findNodeIterAtPoint(octreeContext, absolutePosition);
    if (iter != nodeMap.end()) {
      OctreeNodePtr node = iter->second;
      int parentLod = node->lod;
      
      if (parentLod > targetLod) {
        const vm::ivec2 &parentPosition = node->min;
        int childLod = parentLod / 2;

        // std::cout << "split parent node: " << parentPosition.x << " " << parentPosition.y << " " << parentLod << std::endl;

        nodeMap.erase(iter);

        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            vm::ivec2 min{
              parentPosition.x + dx * childLod,
              parentPosition.y + dy * childLod
            };
            
            OctreeNodePtr childNode = createNode(
              octreeContext,
              min,
              childLod
            );
          }
        }
      } else {
        break;
      }
    } else {
      std::cerr << "could not find split point node: " << absolutePosition.x << " " << absolutePosition.y << std::endl;
      std::cout << "existing nodes: ";
      for (auto iter = nodeMap.begin(); iter != nodeMap.end(); iter++) {
        auto node = iter->second;
        std::cout << node->min.x << " " << node->min.y << " " << node->lod << ", ";
      }
      std::cout << std::endl;
      abort();
    }
  }
}
void initializeLodArrays(OctreeContext &octreeContext) {
  auto &nodeMap = octreeContext.nodeMap;
  
  for (auto iter = nodeMap.begin(); iter != nodeMap.end(); iter++) {
    auto &node = iter->second;
    auto &min = node->min;
    auto &lod = node->lod;
    
    vm::ivec2 bottomNodePosition{
      min.x,
      min.y + lod
    };
    auto bottomNodeIter = findNodeIterAtPoint(octreeContext, bottomNodePosition);

    vm::ivec2 rightNodePosition{
      min.x + lod,
      min.y
    };
    auto rightNodeIter = findNodeIterAtPoint(octreeContext, rightNodePosition);

    node->lodArray[0] = bottomNodeIter != nodeMap.end() ? bottomNodeIter->second->lod : (lod * 2);
    node->lodArray[1] = rightNodeIter != nodeMap.end() ? rightNodeIter->second->lod : (lod * 2);
  }
}
void constructLodTree(
  OctreeContext &octreeContext,
  const vm::ivec2 &currentCoord,
  int lod1Range,
  int minLod,
  int maxLod
) {
  auto &nodeMap = octreeContext.nodeMap;

  // initialize max lod
  vm::ivec2 maxLodCenter{
    (int)std::floor((float)currentCoord.x / (float)maxLod) * maxLod,
    (int)std::floor((float)currentCoord.y / (float)maxLod) * maxLod
  };
  for (int dx = -lod1Range * maxLod; dx <= lod1Range * maxLod; dx += maxLod) {
    for (int dy = -lod1Range * maxLod; dy <= lod1Range * maxLod; dy += maxLod) {
      vm::ivec2 childPosition{
        maxLodCenter.x + dx,
        maxLodCenter.y + dy
      };
      OctreeNodePtr childNode = getOrCreateNode(
        octreeContext,
        childPosition,
        maxLod
      );
    }
  }

  // initialize other lods
  for (int lod = minLod; lod <= maxLod; lod *= 2) {
    for (int dx = -lod1Range * lod; dx <= lod1Range * lod; dx += lod) {
      for (int dz = -lod1Range * lod; dz <= lod1Range * lod; dz += lod) {
        vm::ivec2 currentCoordSnappedToLod{
          (int)std::floor((float)currentCoord.x / (float)lod) * lod,
          (int)std::floor((float)currentCoord.y / (float)lod) * lod
        };
        vm::ivec2 splitPosition{
          currentCoordSnappedToLod.x + dx,
          currentCoordSnappedToLod.y + dz
        };
        splitPointToLod(octreeContext, splitPosition, lod);
      }
    }
  }

  initializeLodArrays(octreeContext);
}
void constructSeedTree(
  OctreeContext &octreeContext,
  const std::vector<vm::ivec2> &maxLodChunkPositions,
  const int maxLod,
  const std::vector<std::pair<vm::ivec2, int>> &lodSplits
) {
  auto &nodeMap = octreeContext.nodeMap;

  // initialize base lods
  for (size_t i = 0; i < maxLodChunkPositions.size(); i++) {
    const vm::ivec2 &nodePos = maxLodChunkPositions[i];
    OctreeNodePtr childNode = getOrCreateNode(
      octreeContext,
      nodePos,
      maxLod
    );
  }

  // initialize other lods
  for (const std::pair<vm::ivec2, int> &lodSplit : lodSplits) {
    const vm::ivec2 &splitPosition = lodSplit.first;
    const int &lod = lodSplit.second;
    splitPointToLod(octreeContext, splitPosition, lod);
  }
}
std::vector<OctreeNodePtr> constructOctreeForLeaf(
  const vm::ivec2 &currentCoord,
  int lod1Range,
  int minLod,
  int maxLod
) {
  OctreeContext octreeContext;

  constructLodTree(
    octreeContext,
    currentCoord,
    lod1Range,
    minLod,
    maxLod
  );

  // collect all nodes in the nodeMap
  std::vector<OctreeNodePtr> leafNodes = octreeContext.getLeafNodes();

  // return
  return leafNodes;
}
/* std::vector<OctreeNodePtr> constructOctreeForSeed(const vm::ivec2 &maxLodCenter, const int maxLod, const std::vector<std::pair<vm::ivec2, int>> &lodSplits) {
  OctreeContext octreeContext;

  constructSeedTree(
    octreeContext,
    currentCoord,
    lod1Range,
    minLod,
    maxLod
  );
  
  std::vector<OctreeNodePtr> leafNodes = octreeContext.getLeafNodes();
  
  return leafNodes;
} */
/* OctreeNodePtr getMaxLodNode(const std::vector<OctreeNodePtr> &newLeafNodes, const std::vector<OctreeNodePtr> &oldLeafNodes, const vm::ivec3 &min) {
    auto newLeafNode = getLeafNodeFromPoint(newLeafNodes, min);
    auto oldLeafNode = getLeafNodeFromPoint(oldLeafNodes, min);
    if (newLeafNode != nullptr && oldLeafNode != nullptr) {
      return newLeafNode->size > oldLeafNode->size ? newLeafNode : oldLeafNode;
    } else if (newLeafNode != nullptr) {
      return newLeafNode;
    } else if (oldLeafNode != nullptr) {
      return oldLeafNode;
    } else {
      return nullptr;
    }
} */
/* std::vector<TrackerTaskPtr> diffLeafNodes(const std::vector<OctreeNodePtr> &newLeafNodes, const std::vector<OctreeNodePtr> &oldLeafNodes) {
  // map from min lod hash to task containing new nodes and old nodes
  std::unordered_map<uint64_t, TrackerTaskPtr> taskMap;

  for (OctreeNodePtr newNode : newLeafNodes) {
    OctreeNodePtr maxLodNode = getMaxLodNode(newLeafNodes, oldLeafNodes, newNode->min);
    const uint64_t hash = hashOctreeMinLod(maxLodNode->min, maxLodNode->size);
    
    TrackerTaskPtr task;
    const auto &iter = taskMap.find(hash);
    if (iter != taskMap.end()) {
      task = iter->second;
    } else {
      TrackerTask *trackerTask = new TrackerTask();
      trackerTask->id = ++nextTrackerId;
      // std::cout << "increment 1 " << trackerTask->id << std::endl;
      // trackerTask->maxLodNode = maxLodNode;
      // trackerTask->type = TrackerTaskType::ADD;
      
      task = std::shared_ptr<TrackerTask>(trackerTask);
      
      taskMap[hash] = task;
    }
    task->newNodes.push_back(newNode);
  }
  for (OctreeNodePtr oldNode : oldLeafNodes) {
    OctreeNodePtr maxLodNode = getMaxLodNode(newLeafNodes, oldLeafNodes, oldNode->min);
    const uint64_t hash = hashOctreeMinLod(maxLodNode->min, maxLodNode->size);
    
    TrackerTaskPtr task;
    const auto &iter = taskMap.find(hash);
    if (iter != taskMap.end()) {
      task = iter->second;
    } else {
      TrackerTask *trackerTask = new TrackerTask();
      trackerTask->id = ++nextTrackerId;
      // std::cout << "increment 2 " << trackerTask->id << std::endl;
      // trackerTask->maxLodNode = maxLodNode;
      // trackerTask->type = TrackerTaskType::REMOVE;
      
      task = std::shared_ptr<TrackerTask>(trackerTask);
      
      taskMap[hash] = task;
    }
    task->oldNodes.push_back(oldNode);
  }

  std::vector<TrackerTaskPtr> tasks;
  for (const auto &iter : taskMap) {
    tasks.push_back(iter.second);
  }
  return tasks;
} */
/* // sort tasks by distance to world position of the central max lod node
std::vector<TrackerTaskPtr> sortTasks(const std::vector<TrackerTaskPtr> &tasks, const vm::vec3 &worldPosition) {
  std::vector<std::pair<TrackerTaskPtr, float>> taskDistances;
  taskDistances.reserve(tasks.size());

  for (const auto &task : tasks) {
    const vm::ivec3 &min = task->maxLodNode->min;
    const int &lod = task->maxLodNode->size;

    vm::vec3 center = vm::vec3{(float)min.x, (float)min.y, (float)min.z} +
      vm::vec3{0.5, 0.5, 0.5} * (float)lod;
    vm::vec3 delta = worldPosition - center;
    float distance = vm::lengthSq(delta);

    taskDistances.push_back(std::pair<TrackerTaskPtr, float>(task, distance));
  }

  std::sort(
    taskDistances.begin(),
    taskDistances.end(),
    [](const std::pair<TrackerTaskPtr, float> &a, const std::pair<TrackerTaskPtr, float> &b) -> bool {
      return a.second < b.second;
    }
  );

  std::vector<TrackerTaskPtr> sortedTasks;
  for (const auto &iter : taskDistances) {
    sortedTasks.push_back(iter.first);
  }
  return sortedTasks;
} */
/* std::pair<std::vector<OctreeNodePtr>, std::vector<TrackerTaskPtr>> updateChunks(const std::vector<OctreeNodePtr> &oldChunks, const std::vector<TrackerTaskPtr> &tasks) {
  std::vector<OctreeNodePtr> newChunks = oldChunks;
  
  // swap old chunks for new chunks
  for (TrackerTaskPtr task : tasks) {
    if (!task->isNop()) {
      const std::vector<OctreeNodePtr> &newNodes = task->newNodes;
      const std::vector<OctreeNodePtr> &oldNodes = task->oldNodes;

      for (OctreeNodePtr oldNode : oldNodes) {
        const auto &iter = std::find_if(
          newChunks.begin(),
          newChunks.end(),
          [&](auto &chunk) -> bool {
            return equalsNode(*chunk, *oldNode);
          }
        );
        if (iter != newChunks.end()) {
          newChunks.erase(iter);
        } else {
          std::cout << "failed to erase old node" << std::endl;
          abort();
        }
      }
      for (OctreeNodePtr newNode : newNodes) {
        newChunks.push_back(newNode);
      }
    }
  }

  // compute extra tasks for outdated chunks
  std::vector<TrackerTaskPtr> extraTasks;

  return std::pair<std::vector<OctreeNodePtr>, std::vector<TrackerTaskPtr>>(
    std::move(newChunks),
    std::move(extraTasks)
  );
} */

//

Tracker::Tracker(PGInstance *inst) :
  inst(inst),
  /* lastCoord{
    INT32_MAX,
    INT32_MAX
  }, */
  epoch(0)
{}
// static methods
vm::ivec2 getCurrentCoord(const vm::vec3 &position, int chunkSize) {
  const int cx = std::floor(position.x / (float)chunkSize);
  const int cz = std::floor(position.z / (float)chunkSize);
  return vm::ivec2{cx, cz};
}
/* bool duplicateTask(const std::vector<TrackerTaskPtr> &tasks) {
  std::unordered_map<int, bool> seen;

  for (TrackerTaskPtr task : tasks) {
    int id = task->id;
    if (seen.find(id) == seen.end()) {
      seen[id] = true;
    } else {
      std::cout << "duplicate task: " <<
        task->maxLodNode->min.x << " " << task->maxLodNode->min.y <<
        id << std::endl;
      abort();
      return true;
    }
  }
  return false;
}
bool duplicateTask(const std::vector<TrackerTaskPtr> &tasks, const std::vector<TrackerTaskPtr> &tasks2) {
  std::unordered_map<int, bool> seen;

  for (TrackerTaskPtr task : tasks) {
    int id = task->id;
    if (seen.find(id) == seen.end()) {
      seen[id] = true;
    } else {
      std::cout << "duplicate task: " <<
        task->maxLodNode->min.x << " " << task->maxLodNode->min.y <<
        id << std::endl;
      abort();
      return true;
    }
  }
  for (TrackerTaskPtr task : tasks2) {
    int id = task->id;
    if (seen.find(id) == seen.end()) {
      seen[id] = true;
    } else {
      std::cout << "duplicate task: " <<
        task->maxLodNode->min.x << " " << task->maxLodNode->min.y << " " << task->maxLodNode->min.z << " " <<
        id << std::endl;
      abort();
      return true;
    }
  }
  return false;
} */
// dynamic methods
// sort nodes by distance to world position of the central max lod node
void Tracker::sortNodes(std::vector<OctreeNodePtr> &nodes) {
  const vm::vec3 &worldPosition = inst->worldPosition;
  const vm::vec3 &cameraPosition = inst->cameraPosition;
  const Quat &cameraQuaternion = inst->cameraQuaternion;
  std::array<float, 16> &projectionMatrix = inst->projectionMatrix;
  
  // compute frustum
  Matrix matrixWorld(
    Vec{
      cameraPosition.x,
      cameraPosition.y,
      cameraPosition.z
    },
    Quat{
      cameraQuaternion.x,
      cameraQuaternion.y,
      cameraQuaternion.z,
      cameraQuaternion.w
    },
    Vec{1, 1, 1}
  );
  Matrix matrixWorldInverse(matrixWorld);
  matrixWorldInverse.invert();
  Frustum frustum = Frustum::fromMatrix(
    Matrix::fromArray(projectionMatrix.data()) *= matrixWorldInverse
  );

  sort<OctreeNodePtr>(nodes, worldPosition, frustum);
}

/* const OctreeNodePtr findLeafNodeForPosition(
  const std::vector<OctreeNodePtr> &nodes,
  const vm::ivec3 &p
) {
  for (auto iter = nodes.begin(); iter != nodes.end(); iter++) {
    const OctreeNodePtr &node = *iter;
    if (containsPoint(*node, p)) {
      return node;
    }
  }
  return nullptr;
} */

DataRequestUpdate Tracker::updateDataRequests(
  const std::unordered_map<uint64_t, DataRequestPtr> &dataRequests,
  const std::vector<OctreeNodePtr> &leafNodes
) {
  DataRequestUpdate dataRequestUpdate;
  dataRequestUpdate.dataRequests = dataRequests;

  // cancel old data requests
  for (auto iter = dataRequestUpdate.dataRequests.begin(); iter != dataRequestUpdate.dataRequests.end();) {
    const uint64_t &hash = iter->first;
    const DataRequestPtr &oldDataRequest = iter->second;

    auto matchingLeafNodeIter = std::find_if(
      leafNodes.begin(),
      leafNodes.end(),
      [&](const OctreeNodePtr &leafNode) -> bool {
        return equalsNodeLodArray(*leafNode, *oldDataRequest->node);
      }
    );
    if (matchingLeafNodeIter != leafNodes.end()) {
      // keep the data request
      dataRequestUpdate.keepDataRequests.push_back(oldDataRequest);

      iter++;
    } else {
      // cancel the data request
      dataRequestUpdate.cancelDataRequests.push_back(oldDataRequest);
      // forget the data request
      auto currentIter = iter;
      auto nextIter = iter;
      nextIter++;
      dataRequestUpdate.dataRequests.erase(currentIter);
      iter = nextIter;
    }
  }

  // add new data requests
  for (auto chunk : leafNodes) {
    uint64_t hash = hashOctreeMinLod(chunk->min, chunk->lod);
    auto dataRequestIter = dataRequestUpdate.dataRequests.find(hash);
    if (dataRequestIter == dataRequestUpdate.dataRequests.end()) {
      DataRequestPtr dataRequest(new DataRequest{chunk});

      dataRequestUpdate.newDataRequests.push_back(dataRequest);

      dataRequestUpdate.dataRequests.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(hash),
        std::forward_as_tuple(dataRequest)
      );
    }
  }

  return dataRequestUpdate;
}
/* std::unordered_map<uint64_t, Dominator> Tracker::updateDominators(
  const std::vector<OctreeNodePtr> &oldRenderedChunks,
  const std::vector<OctreeNodePtr> &newRenderedChunks
) {
  // compute dominators. a dominator is a chunk lod that contains a a transform from old chunks to new chunks.
  std::unordered_map<uint64_t, Dominator> dominators;
  auto addChunk = [&](const OctreeNodePtr &chunk, bool isNew) -> void {
    OctreeNodePtr oldLeafNode = findLeafNodeForPosition(oldRenderedChunks, chunk->min);
    OctreeNodePtr newLeafNode = findLeafNodeForPosition(newRenderedChunks, chunk->min);
    OctreeNodePtr maxLodChunk = nullptr;
    if (oldLeafNode != nullptr && newLeafNode != nullptr) {
      maxLodChunk = oldLeafNode->size > newLeafNode->size ? oldLeafNode : newLeafNode;
    } else if (oldLeafNode != nullptr) {
      maxLodChunk = oldLeafNode;
    } else if (newLeafNode != nullptr) {
      maxLodChunk = newLeafNode;
    } else {
      std::cout << "no chunk match in any leaf set" << std::endl;
      abort();
    }

    uint64_t maxLodHash = hashOctreeMinLod(maxLodChunk->min, maxLodChunk->lod);
    auto dominatorIter = dominators.find(maxLodHash);
    if (dominatorIter == dominators.end()) {
      auto insertResult = dominators.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(maxLodHash),
        std::forward_as_tuple(maxLodChunk)
      );
      dominatorIter = insertResult.first;
    }
    Dominator &dominator = dominatorIter->second;
    if (isNew) {
      dominator.newChunks.push_back(chunk);
    } else {
      dominator.oldChunks.push_back(chunk);
    }
  };
  for (auto chunk : oldRenderedChunks) {
    addChunk(chunk, false);
  }
  for (auto chunk : newRenderedChunks) {
    addChunk(chunk, true);
  }
  return dominators;
} */
TrackerUpdate Tracker::update(const vm::vec3 &position, int lods, int lod1Range) {
  int oldEpoch;
  {
    std::lock_guard<std::mutex> lock(mutex);
    oldEpoch = epoch;
  }

  // new octrees
  int chunkSize = this->inst->chunkSize;
  vm::ivec2 currentCoord = getCurrentCoord(position, chunkSize); // in chunk space
  int lod = 1 << (lods - 1);
  std::vector<OctreeNodePtr> octreeLeafNodes = constructOctreeForLeaf(
    currentCoord,
    lod1Range,
    1,
    lod
  );
  sortNodes(octreeLeafNodes);

  DataRequestUpdate dataRequestUpdate = updateDataRequests(this->dataRequests, octreeLeafNodes);

  //

  // new rendered chunks come from the data requests
  std::vector<OctreeNodePtr> newRenderedChunks(dataRequests.size());
  std::transform(
    dataRequests.begin(),
    dataRequests.end(),
    newRenderedChunks.begin(),
    [](const auto &iter) -> const auto & {
      return iter.second->node;
    }
  );

  //

  {
    std::lock_guard<std::mutex> lock(mutex);
  
    int currentEpoch = epoch;
    if (currentEpoch == oldEpoch) {
      dataRequests = std::move(dataRequestUpdate.dataRequests);
    }
  }

  //

  TrackerUpdate result;

  result.leafNodes = std::move(octreeLeafNodes);

  result.newDataRequests = std::move(dataRequestUpdate.newDataRequests);
  result.keepDataRequests = std::move(dataRequestUpdate.keepDataRequests);
  result.cancelDataRequests = std::move(dataRequestUpdate.cancelDataRequests);
  // result.chunkPosition = node->min / chunkSize;

  //

  return result;
}
/* TrackerUpdate Tracker::update(const vm::vec3 &position) {
  const vm::ivec2 &currentCoord = getCurrentCoord(position);
  TrackerUpdate trackerUpdate = updateCoord(currentCoord);
  return trackerUpdate;
} */