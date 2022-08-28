#include "tracker.h"
#include "octree.h"
#include "sort.h"
#include <iostream>

/*
tracker = await dcWorkerManager.createTracker(3, 2, true); trackerUpdate = await dcWorkerManager.trackerUpdate(tracker, new THREE.Vector3(1, 2, 3));
*/

//

std::atomic<int> nextTrackerId(0);

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
uint8_t *TrackerUpdate::getBuffer() const {
  // std::vector<std::vector<uint8_t>> oldTaskBuffers;
  /* for (const auto &task : oldTasks) {
    oldTaskBuffers.push_back(task->getBuffer());
  } */

  /* std::vector<std::vector<uint8_t>> newTaskBuffers;
  for (const auto &task : newTasks) {
    newTaskBuffers.push_back(task->getBuffer());
  } */

  size_t size = 0;
  // size += sizeof(vm::ivec3); // currentCoord
  // size += sizeof(uint32_t); // numOldTasks
  // size += sizeof(uint32_t); // numNewTasks
  size += sizeof(int32_t); // numLeafNodes
  /* for (auto &buffer : oldTaskBuffers) {
    size += buffer.size();
  }
  for (auto &buffer : newTaskBuffers) {
    size += buffer.size();
  } */
  // leaf nodes
  size += sizeof(vm::ivec3) * leafNodes.size(); // min
  size += sizeof(int) * leafNodes.size(); // lod
  // size += sizeof(int) * leafNodes.size(); // isLeaf
  size += sizeof(int) * leafNodes.size(); // lodArray

  uint8_t *ptr = (uint8_t *)malloc(size);
  int index = 0;
  // numLeafNodes
  *((int32_t *)(ptr + index)) = leafNodes.size();
  
  index += sizeof(int32_t);
  // leaf nodes
  for (auto leafNode : leafNodes) {
    std::memcpy(ptr + index, &leafNode->min, sizeof(vm::ivec3));
    index += sizeof(vm::ivec3);
    *((int *)(ptr + index)) = leafNode->lod;
    index += sizeof(int);
    // *((int *)(ptr + index)) = (leafNode->type == Node_Leaf) ? 1 : 0;
    // index += sizeof(int);
    std::memcpy(ptr + index, &leafNode->lod, sizeof(int));
    index += sizeof(int);
  }
  return ptr;
}

//

std::vector<uint8_t> DataRequest::getBuffer() const {
  size_t size = 0;
  size += sizeof(vm::ivec3); // min
  size += sizeof(int); // size
  size += sizeof(int); // isLeaf
  size += sizeof(int[8]); // lodArray

  std::vector<uint8_t> result(size);
  int index = 0;
  // min
  std::memcpy(result.data() + index, &node->min, sizeof(vm::ivec3));
  index += sizeof(vm::ivec3);
  // size
  *((int *)(result.data() + index)) = node->lod;
  index += sizeof(int);
  // isLeaf
  // *((int *)(result.data() + index)) = (node->type == Node_Leaf) ? 1 : 0;
  // index += sizeof(int);
  // lodArray
  std::memcpy(result.data() + index, &node->lod, sizeof(int));
  index += sizeof(int);
  return result;
}

std::vector<uint8_t> DataRequestUpdate::getBuffer() const {
  return std::vector<uint8_t>(); // XXX
}

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

bool containsPoint(const OctreeNode &node, const vm::ivec2 &p) {
  return p.x >= node.min.x && p.x < node.min.x + node.lod &&
    p.y >= node.min.y && p.y < node.min.y + node.lod;
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
/* bool intersectsNode(const OctreeNode &node, const OctreeNode &other) {
  return containsNode(node, other) || containsNode(other, node);
} */

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

//

OctreeNodePtr getLeafNodeFromPoint(const std::vector<OctreeNodePtr> &leafNodes, const vm::ivec2 &p) {
    for (size_t i = 0; i < leafNodes.size(); i++) {
        auto leafNode = leafNodes[i];
        if (containsPoint(*leafNode, p)) {
            return leafNode;
        }
    }
    return nullptr;
}
OctreeNodePtr getNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod) {
    auto &nodeMap = octreeContext.nodeMap;

    uint64_t hash = hashOctreeMinLod(min, lod);
    auto iter = nodeMap.find(hash);
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
void constructTreeUpwards(OctreeContext &octreeContext, const vm::ivec2 &leafPosition, int minLod, int maxLod) {
    auto &nodeMap = octreeContext.nodeMap;

    OctreeNodePtr rootNode = getOrCreateNode(octreeContext, leafPosition, minLod);
    vm::ivec2 rootNodeMin = rootNode->min;    
    for (int lod = minLod * 2; lod <= maxLod; lod *= 2) {
      rootNodeMin.x = (int)std::floor((float)rootNodeMin.x / (float)lod) * lod;
      rootNodeMin.y = (int)std::floor((float)rootNodeMin.y / (float)lod) * lod;

      vm::ivec2 lodMax = rootNodeMin + vm::ivec2{lod, lod};

      for (int bx = rootNodeMin.x; bx < lodMax.x; bx += lod) {
        for (int bz = rootNodeMin.y; bz < lodMax.y; bz += lod) {
          vm::ivec2 min = {bx, bz};

          // if it's not the node we are rising from, create it
          if (min != rootNodeMin) {
            OctreeNodePtr node = getOrCreateNode(octreeContext, min, lod);
            // node->parent = rootNode;
            // rootNode->children.push_back(node);
          }

          // OctreeNodePtr node = getOrCreateNode(octreeContext, vm::ivec2{bx, by}, lod);
          // node->parent = rootNode;
          // rootNode->children.push_back(node);
        }
      }

      /* const vm::ivec2 &lodCenter = lodMin + (lod / 2);
      const int childIndex = (rootNode->min.x < lodCenter.x ? 0 : 1) +
        (rootNode->min.y < lodCenter.y ? 0 : 2) +
        (rootNode->min.z < lodCenter.z ? 0 : 4);

      OctreeNodePtr parentNode = getOrCreateNode(octreeContext, lodMin, lod);
      if (parentNode->type == Node_Leaf) {
        parentNode->type = Node_Internal;
        for (int i = 0; i < 8; i++) {
          parentNode->children[i] = nullptr;
        }
      } else if (parentNode->type == Node_Internal) {
        // nothing
      } else {
        std::cout << "Invalid node type: " << (int)parentNode->type << std::endl;
        abort();
      }
      
      // set child as we walk upwards
      parentNode->children[childIndex] = rootNode.get();

      // ensure that this internal node has at least leaf nodes
      if (flags & ConstructTreeFlags_Children) {
        ensureChildren(octreeContext, parentNode.get());
      }

      // ensure that the space around is at least this lod
      if (flags & ConstructTreeFlags_Peers) {
        ensurePeers(octreeContext, parentNode.get(), maxLod);
      } */

      // walk upwards
      // rootNode = parentNode;
    }
    // return rootNode;
}
/* // ensure that every neighbor of this lod also is at least at this lod (it could be a lower/more detailed leaf)
void ensurePeers(OctreeContext &octreeContext, OctreeNode *node, int maxLod) {
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dx == 0 && dy == 0 && dz == 0) {
          continue;
        }

        vm::ivec3 neighborMin = node->min + vm::ivec3{dx, dy, dz} * node->size;
        constructTreeUpwards(
          octreeContext,
          neighborMin,
          node->size,
          maxLod,
          ConstructTreeFlags_None
        );
      }
    }
  }
} */

/* void stdoutSpaces(int numSpaces) {
    for (int i = 0; i < numSpaces; i++) {
        printf(" ");
    }
} */
std::vector<OctreeNodePtr> constructOctreeForLeaf(const vm::ivec2 &position, int lod1Range, int maxLod) {
  OctreeContext octreeContext;
  auto &nodeMap = octreeContext.nodeMap;

  constexpr int minLod = 1;

  // sample base leaf nodes to generate octree upwards
  // 1x lod
  const vm::ivec2 &rangeMin = position - vm::ivec2{lod1Range, lod1Range};
  const vm::ivec2 &rangeMax = position + vm::ivec2{lod1Range, lod1Range};
  for (int dx = rangeMin.x; dx <= rangeMax.x; dx++) {
    for (int dy = rangeMin.y; dy <= rangeMax.y; dy++) {
      vm::ivec2 leafPosition = vm::ivec2{dx, dy};
      constructTreeUpwards(
        octreeContext,
        leafPosition,
        minLod,
        maxLod
      );
      /* const int lod = 1;
      std::shared_ptr<OctreeNode> node = std::make_shared<OctreeNode>(new OctreeNode(
        leafPosition,
        lod
      )); */
    }
  }

  std::vector<OctreeNodePtr> lod1Nodes;
  for (const auto &iter : nodeMap) {
    auto node = iter.second;
    if (node->lod == 1) {
      lod1Nodes.push_back(node);
    }
  }

  // collect all nodes in the nodeMap
  std::vector<OctreeNodePtr> leafNodes;
  for (const auto &iter : nodeMap) {
    auto node = iter.second;
    leafNodes.push_back(node);
  }

  // sanity check lod1Nodes for duplicates
  {
    std::unordered_map<uint64_t, OctreeNodePtr> lod1NodeMap;
    for (auto node : lod1Nodes) {
      uint64_t hash = hashOctreeMinLod(node->min, node->lod);
      if (lod1NodeMap.find(hash) != lod1NodeMap.end()) {
        // throw new Error(`Duplicate lod1 node: ${hash}`);
        EM_ASM({
            console.log('Duplicate lod1 node:', $0);
        }, hash);
        abort();
      }
      lod1NodeMap[hash] = node;
    }
  }

  // sanity check that no leaf node contains another leaf node
  for (auto leafNode : leafNodes) {
    for (auto leafNode2 : leafNodes) {
      if (leafNode != leafNode2 && containsNode(*leafNode, *leafNode2)) {
        EM_ASM({
            console.log('Leaf node contains another leaf node 2:', $0, $1);
        }, leafNode->min.x, leafNode->min.y);
      }
    }
  }

  /* // assign lodArray for each node based on the minimum lod of the target point in the world
  for (const auto &iter : nodeMap) {
    auto node = iter.second;

    for (int i = 0; i < 8; i++) {
      const vm::ivec3 &offset = lodOffsets[i];
      int &lod = node->lodArray[i]; // output
      
      const vm::ivec3 &p = node->min + offset * node->size;
      OctreeNodePtr containingLeafNode = getLeafNodeFromPoint(leafNodes, p);
      if (containingLeafNode) {
        lod = containingLeafNode->size;
      } else {
        lod = node->size;
      }
    }
  } */

  // return
  return leafNodes;
}
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

Tracker::Tracker(PGInstance *inst, int lods, int lod1Range) :
  inst(inst),
  lods(lods),
  lod1Range(lod1Range),
  lastCoord{
    INT32_MAX,
    INT32_MAX
  }
{}
// static methods
vm::ivec2 Tracker::getCurrentCoord(const vm::vec3 &position) {
  const int &chunkSize = inst->chunkSize;
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

const OctreeNodePtr findLeafNodeForPosition(
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
}

DataRequestUpdate Tracker::updateDataRequests(
  std::unordered_map<uint64_t, DataRequestPtr> &dataRequests,
  const std::vector<OctreeNodePtr> &leafNodes
) {
  DataRequestUpdate dataRequestUpdate;

  // cancel old data requests
  for (auto iter = dataRequests.begin(); iter != dataRequests.end();) {
    const uint64_t &hash = iter->first;
    const DataRequestPtr &oldDataRequest = iter->second;

    auto matchingLeafNodeIter = std::find_if(
      leafNodes.begin(),
      leafNodes.end(),
      [&](const OctreeNodePtr &leafNode) -> bool {
        return equalsNodeLod(*leafNode, *oldDataRequest->node);
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
      dataRequests.erase(currentIter);
      iter = nextIter;
    }
  }

  // add new data requests
  for (auto chunk : leafNodes) {
    uint64_t hash = hashOctreeMinLod(chunk->min, chunk->lod);
    auto dataRequestIter = dataRequests.find(hash);
    if (dataRequestIter == dataRequests.end()) {
      DataRequestPtr dataRequest(new DataRequest{chunk});

      dataRequestUpdate.newDataRequests.push_back(dataRequest);

      dataRequests.emplace(
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
TrackerUpdate Tracker::updateCoord(const vm::ivec2 &currentCoord) {
  // old octrees
  // const std::vector<OctreeNodePtr> &oldRenderedChunks = this->leafNodes;
  
  // new octrees
  std::vector<OctreeNodePtr> octreeLeafNodes = constructOctreeForLeaf(
    currentCoord,
    this->lod1Range,
    1 << (this->lods - 1)
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

  // compute dominators
  // std::unordered_map<uint64_t, Dominator> dominators = updateDominators(oldRenderedChunks, newRenderedChunks);

  //

  // this->leafNodes = octreeLeafNodes;

  //

  TrackerUpdate result;

  result.leafNodes = std::move(octreeLeafNodes);

  result.newDataRequests = std::move(dataRequestUpdate.newDataRequests);
  result.keepDataRequests = std::move(dataRequestUpdate.keepDataRequests);
  result.cancelDataRequests = std::move(dataRequestUpdate.cancelDataRequests);

  // result.dominators = std::move(dominators);

  //

  return result;
}
TrackerUpdate Tracker::update(const vm::vec3 &position) {
  const vm::ivec2 &currentCoord = getCurrentCoord(position);
  TrackerUpdate trackerUpdate = updateCoord(currentCoord);
  return trackerUpdate;
}
