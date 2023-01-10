#include "tracker.h"
#include "../instance.h"
#include <iostream>

void serializeOctreeNodes(const std::vector<OctreeNodePtr> &datas, uint8_t *ptr, int &index)
{
  // std::cout << "serialize octree nodes" << datas.size() << std::endl;

  *((int32_t *)(ptr + index)) = datas.size();
  index += sizeof(int32_t);

  for (const auto &data : datas)
  {
    std::memcpy(ptr + index, &data->min, sizeof(vm::ivec2));
    index += sizeof(vm::ivec2);

    *((int *)(ptr + index)) = data->lod;
    index += sizeof(int);

    std::memcpy(ptr + index, data->lodArray, sizeof(int[2]));
    index += sizeof(int[2]);
  }
}

void serializeDataRequests(const std::vector<DataRequestPtr> &datas, uint8_t *ptr, int &index)
{
  // std::cout << "serialize data requests" << datas.size() << std::endl;

  *((int32_t *)(ptr + index)) = datas.size();
  index += sizeof(int32_t);

  for (const auto &data : datas)
  {
    std::memcpy(ptr + index, &data->node->min, sizeof(vm::ivec2));
    index += sizeof(vm::ivec2);

    *((int *)(ptr + index)) = data->node->lod;
    index += sizeof(int);

    std::memcpy(ptr + index, data->node->lodArray, sizeof(int[2]));
    index += sizeof(int[2]);

    *((int *)(ptr + index)) = (int)data->node->replacing;
    index += sizeof(int);
  }
}

uint8_t *TrackerUpdate::getBuffer() const
{
  // compute size
  size_t size = 0;

  constexpr size_t leafNodeSize = sizeof(vm::ivec2) + sizeof(int) + sizeof(int[2]);
  constexpr size_t octreeNodeSize = sizeof(vm::ivec2) + sizeof(int) + sizeof(int[2]) + sizeof(int);

  size += sizeof(int32_t); // numLeafNodes
  size += octreeNodeSize * leafNodes.size();

  size += sizeof(int32_t); // numNewDataRequests
  size += octreeNodeSize * newDataRequests.size();

  size += sizeof(int32_t); // numReplacingRequests
  size += octreeNodeSize * replacingRequests.size();

  size += sizeof(int32_t); // numCancelDataRequests
  size += octreeNodeSize * cancelDataRequests.size();

  // size += sizeof(vm::ivec2); // chunkPosition

  // serialize
  uint8_t *ptr = (uint8_t *)malloc(size);
  int index = 0;

  serializeOctreeNodes(leafNodes, ptr, index);
  serializeDataRequests(newDataRequests, ptr, index);
  serializeDataRequests(replacingRequests, ptr, index);
  serializeDataRequests(cancelDataRequests, ptr, index);

  return ptr;
}

bool containsPoint(const vm::ivec2 &min, const int lod, const vm::ivec2 &p)
{
  return p.x >= min.x && p.x < min.x + lod &&
         p.y >= min.y && p.y < min.y + lod;
}
bool containsPoint(const vm::vec2 &min, const int lod, const vm::vec2 &p)
{
  return p.x >= min.x && p.x < min.x + (float)lod &&
         p.y >= min.y && p.y < min.y + (float)lod;
}
bool containsPoint(const OctreeNode &node, const vm::ivec2 &p)
{
  return containsPoint(node.min, node.lod, p);
}
bool containsNode(const OctreeNode &node, const OctreeNode &other)
{
  return containsPoint(node, other.min);
}
bool equalsNodeMin(const OctreeNode &node, const OctreeNode &other)
{
  return node.min == other.min;
}
bool equalsLod(const OctreeNode &node, const OctreeNode &other)
{
  return node.lod == other.lod;
}

bool equalsNodeMinLod(const OctreeNode &node, const OctreeNode &other)
{
  return equalsNodeMin(node, other) && equalsLod(node, other);
}

bool equalsLodArray(const OctreeNode &node, const OctreeNode &other)
{
  return (node.lodArray[0] == other.lodArray[0]) && (node.lodArray[1] == other.lodArray[1]);
}

bool equalsNode(const OctreeNode &node, const OctreeNode &other)
{
  return equalsNodeMinLod(node, other) && equalsLodArray(node, other);
}

bool chunksIntersect(const vm::ivec2 &min1, const int lod1, const vm::ivec2 &min2, const int lod2)
{
  return min1.x < min2.x + lod2 && min1.x + lod1 > min2.x &&
         min1.y < min2.y + lod2 && min1.y + lod1 > min2.y;
}

OctreeNodePtr OctreeNodeAllocator::alloc(const vm::ivec2 &min, int lod)
{
  return std::make_shared<OctreeNode>(min, lod);
}

OctreeNodePtr OctreeContext::alloc(const vm::ivec2 &min, int lod)
{
  auto node = OctreeNodeAllocator::alloc(min, lod);
  if (node)
  {
    node->min = min;
    node->lod = lod;

    const uint64_t hash = hashOctreeMinLod(min, lod);
    nodeMap[hash] = node;
  }
  return node;
}
// OctreeContext::OctreeContext() {}
std::vector<OctreeNodePtr> OctreeContext::getLeafNodes()
{
  std::vector<OctreeNodePtr> leafNodes;
  for (const auto &iter : nodeMap)
  {
    auto node = iter.second;
    leafNodes.push_back(node);
  }
  return leafNodes;
}

std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>>::iterator getNodeIter(OctreeContext &octreeContext, const vm::ivec2 &min, int lod)
{
  auto &nodeMap = octreeContext.nodeMap;

  const uint64_t hash = hashOctreeMinLod(min, lod);
  auto iter = nodeMap.find(hash);
  return iter;
}
OctreeNodePtr getNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod)
{
  auto &nodeMap = octreeContext.nodeMap;

  auto iter = getNodeIter(octreeContext, min, lod);
  if (iter != nodeMap.end())
  {
    return iter->second;
  }
  else
  {
    return nullptr;
  }
}
OctreeNodePtr createNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod)
{
  auto &nodeMap = octreeContext.nodeMap;

  auto node = OctreeNodeAllocator::alloc(min, lod);
  const uint64_t hash = hashOctreeMinLod(min, lod);
  nodeMap[hash] = node;
  return node;
}
OctreeNodePtr getOrCreateNode(OctreeContext &octreeContext, const vm::ivec2 &min, int lod)
{
  OctreeNodePtr node = getNode(octreeContext, min, lod);
  if (!node)
  {
    node = createNode(octreeContext, min, lod);
  }
  return node;
}

std::unordered_map<uint64_t, std::shared_ptr<OctreeNode>>::iterator findNodeIterAtPoint(OctreeContext &octreeContext, const vm::ivec2 &position)
{
  auto &nodeMap = octreeContext.nodeMap;

  // std::vector<OctreeNodePtr> leafNodes;
  for (auto iter = nodeMap.begin(); iter != nodeMap.end(); iter++)
  {
    auto node = iter->second;
    if (containsPoint(node->min, node->lod, position))
    {
      return iter;
    }
  }
  return nodeMap.end();
}
OctreeNodePtr findContainedNode(OctreeContext &octreeContext, const vm::ivec2 &min, const int lod)
{
  auto &nodeMap = octreeContext.nodeMap;

  // std::vector<OctreeNodePtr> leafNodes;
  for (const auto &iter : nodeMap)
  {
    auto node = iter.second;
    if (containsPoint(min, lod, node->min))
    {
      return node;
    }
  }
  return nullptr;
}
bool removeNode(OctreeContext &octreeContext, OctreeNodePtr node)
{
  auto &nodeMap = octreeContext.nodeMap;

  const uint64_t hash = hashOctreeMinLod(node->min, node->lod);
  auto iter = nodeMap.find(hash);
  if (iter != nodeMap.end())
  {
    nodeMap.erase(iter);
    return true;
  }
  else
  {
    std::cerr << "erase node not found: " << node->min.x << " " << node->min.y << " " << node->lod << std::endl;
    abort();
    return false;
  }
}

void splitPointToLod(
    OctreeContext &octreeContext,
    const vm::ivec2 &absolutePosition,
    int targetLod)
{
  auto &nodeMap = octreeContext.nodeMap;

  for (;;)
  {
    auto iter = findNodeIterAtPoint(octreeContext, absolutePosition);
    if (iter != nodeMap.end())
    {
      OctreeNodePtr node = iter->second;
      int parentLod = node->lod;

      if (parentLod > targetLod)
      {
        const vm::ivec2 &parentPosition = node->min;
        int childLod = parentLod / 2;

        // std::cout << "split parent node: " << parentPosition.x << " " << parentPosition.y << " " << parentLod << std::endl;

        nodeMap.erase(iter);

        for (int dx = 0; dx < 2; dx++)
        {
          for (int dy = 0; dy < 2; dy++)
          {
            vm::ivec2 min{
                parentPosition.x + dx * childLod,
                parentPosition.y + dy * childLod};

            OctreeNodePtr childNode = createNode(
                octreeContext,
                min,
                childLod);
          }
        }
      }
      else
      {
        break;
      }
    }
    else
    {
      std::cerr << "could not find split point node: " << absolutePosition.x << " " << absolutePosition.y << std::endl;
      std::cout << "existing nodes: ";
      for (auto iter = nodeMap.begin(); iter != nodeMap.end(); iter++)
      {
        auto node = iter->second;
        std::cout << node->min.x << " " << node->min.y << " " << node->lod << ", ";
      }
      std::cout << std::endl;
      abort();
    }
  }
}
void initializeLodArrays(OctreeContext &octreeContext)
{
  auto &nodeMap = octreeContext.nodeMap;

  for (auto iter = nodeMap.begin(); iter != nodeMap.end(); iter++)
  {
    auto &node = iter->second;
    auto &min = node->min;
    auto &lod = node->lod;

    vm::ivec2 bottomNodePosition{
        min.x,
        min.y + lod};
    auto bottomNodeIter = findNodeIterAtPoint(octreeContext, bottomNodePosition);

    vm::ivec2 rightNodePosition{
        min.x + lod,
        min.y};
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
    int maxLod)
{
  auto &nodeMap = octreeContext.nodeMap;

  int minLodInt = 1 << (minLod - 1);
  int maxLodInt = 1 << (maxLod - 1);

  // initialize max lod
  vm::ivec2 maxLodCenter{
      (int)std::floor((float)currentCoord.x / (float)maxLodInt) * maxLodInt,
      (int)std::floor((float)currentCoord.y / (float)maxLodInt) * maxLodInt};
  for (int dx = -lod1Range * maxLodInt; dx <= lod1Range * maxLodInt; dx += maxLodInt)
  {
    for (int dy = -lod1Range * maxLodInt; dy <= lod1Range * maxLodInt; dy += maxLodInt)
    {
      vm::ivec2 childPosition{
          maxLodCenter.x + dx,
          maxLodCenter.y + dy};
      OctreeNodePtr childNode = getOrCreateNode(
          octreeContext,
          childPosition,
          maxLodInt);
    }
  }

  // initialize other lods
  if (minLod != maxLod)
  {
    for (int lod = minLodInt; lod <= maxLodInt; lod *= 2)
    {
      for (int dx = -lod1Range * lod; dx <= lod1Range * lod; dx += lod)
      {
        for (int dz = -lod1Range * lod; dz <= lod1Range * lod; dz += lod)
        {
          vm::ivec2 currentCoordSnappedToLod{
              (int)std::floor((float)currentCoord.x / (float)lod) * lod,
              (int)std::floor((float)currentCoord.y / (float)lod) * lod};
          vm::ivec2 splitPosition{
              currentCoordSnappedToLod.x + dx,
              currentCoordSnappedToLod.y + dz};
          splitPointToLod(octreeContext, splitPosition, lod);
        }
      }
    }
  }

  initializeLodArrays(octreeContext);
}
void constructSeedTree(
    OctreeContext &octreeContext,
    const std::vector<vm::ivec2> &maxLodChunkPositions,
    const int maxLod,
    const std::vector<std::pair<vm::ivec2, int>> &lodSplits)
{
  auto &nodeMap = octreeContext.nodeMap;

  // initialize base lods
  for (size_t i = 0; i < maxLodChunkPositions.size(); i++)
  {
    const vm::ivec2 &nodePos = maxLodChunkPositions[i];
    OctreeNodePtr childNode = getOrCreateNode(
        octreeContext,
        nodePos,
        maxLod);
  }

  // initialize other lods
  for (const std::pair<vm::ivec2, int> &lodSplit : lodSplits)
  {
    const vm::ivec2 &splitPosition = lodSplit.first;
    const int &lod = lodSplit.second;
    splitPointToLod(octreeContext, splitPosition, lod);
  }
}
std::vector<OctreeNodePtr> constructOctreeForLeaf(
    const vm::ivec2 &currentCoord,
    int lod1Range,
    int minLod,
    int maxLod)
{
  OctreeContext octreeContext;

  constructLodTree(
      octreeContext,
      currentCoord,
      lod1Range,
      minLod,
      maxLod);

  // collect all nodes in the nodeMap
  std::vector<OctreeNodePtr> leafNodes = octreeContext.getLeafNodes();

  // return
  return leafNodes;
}

Tracker::Tracker(PGInstance *inst) : inst(inst)
{
}
// static methods
vm::ivec2 getCurrentCoord(const vm::vec3 &position, int chunkSize)
{
  const int cx = std::floor(position.x / (float)chunkSize);
  const int cz = std::floor(position.z / (float)chunkSize);
  return vm::ivec2{cx, cz};
}

// dynamic methods
// sort nodes by distance to world position of the central max lod node
void Tracker::sortNodes(std::vector<OctreeNodePtr> &nodes)
{
  const vm::vec3 &worldPosition = inst->renderingInfo.worldPosition;
  const vm::vec3 &cameraPosition = inst->renderingInfo.cameraPosition;
  const Quat &cameraQuaternion = inst->renderingInfo.cameraQuaternion;
  std::array<float, 16> &projectionMatrix = inst->renderingInfo.projectionMatrix;

  // compute frustum
  Matrix matrixWorld(
      Vec{
          cameraPosition.x,
          cameraPosition.y,
          cameraPosition.z},
      Quat{
          cameraQuaternion.x,
          cameraQuaternion.y,
          cameraQuaternion.z,
          cameraQuaternion.w},
      Vec{1, 1, 1});
  Matrix matrixWorldInverse(matrixWorld);
  matrixWorldInverse.invert();
  Frustum frustum = Frustum::fromMatrix(
      Matrix::fromArray(projectionMatrix.data()) *= matrixWorldInverse);

  sort<OctreeNodePtr>(nodes, worldPosition, frustum);
}

template <REPLACING_SIZE r>
void addReplacingNode(std::unordered_map<uint64_t, OctreeNode> &replacingNodes, const vm::ivec2 &min, const int &lod)
{
  const OctreeNode node(min, lod, r);
  const uint64_t hash = hashOctreeMinLod(min, lod);
  replacingNodes.emplace(hash, node);
}

DataRequestPtr createNewDataRequest(DataRequestUpdate &dataRequestUpdate, const OctreeNodePtr &node, const bool &isReplacing)
{
  const uint64_t hash = hashOctreeMinLod(node->min, node->lod);
  DataRequestPtr dataRequest = std::make_shared<DataRequest>(node);
  dataRequestUpdate.dataRequests[hash] = dataRequest;
  if (isReplacing)
  {
    dataRequestUpdate.replacingRequests.push_back(dataRequest);
  }
  return dataRequest;
}

void addNewDataRequest(DataRequestUpdate &dataRequestUpdate, const DataRequestPtr &dataRequest, const bool &isReplacing)
{
  // replacing nodes won't get added right away, they get added once the data is ready (handled in Javascript)
  if (!isReplacing)
  {
    dataRequestUpdate.newDataRequests.push_back(dataRequest);
  }
}

void injectReplacingInfo(OctreeNodePtr &node, const auto &replacingNodeIter, const bool &isReplacing)
{
  // if node is replacing, inject replacing info
  if (isReplacing)
  {
    const OctreeNode *nodeInfo = &replacingNodeIter->second;
    node->replacing = nodeInfo->replacing;
  }
}

template <REPLACING_SIZE r>
void handleReplacingRequest(DataRequestUpdate &dataRequestUpdate, std::unordered_map<uint64_t, OctreeNode> &replacingNodes, const OctreeNode &oldNode, const OctreeNode &newNode)
{
  switch (r)
  {
  case REPLACING_SIZE::SAME:
  {
    // replacing chunk, same size
    const vm::ivec2 min = newNode.min;
    const int lod = newNode.lod;

    addReplacingNode<r>(replacingNodes, min, lod);
    break;
  }
  case REPLACING_SIZE::SMALLER:
  {
    // new chunk is bigger than the old chunk
    /* ...........      ...........
       . -  . -  .      .         .
       ...........  ->  .    +    .
       . -  .  - .      .         .
       ...........      ...........  */

    const vm::ivec2 min = newNode.min;
    const int lod = newNode.lod;

    addReplacingNode<r>(replacingNodes, min, lod);
    break;
  }
  case REPLACING_SIZE::BIGGER:
  {
    // old chunk is bigger than the new chunk

    /* ...........      ...........
       .         .      . +  . +  .
       .    -    .  ->  ...........
       .         .      . +  .  + .
       ...........      ........... */

    const int DIMENSION = 2;

    for (int y = 0; y < DIMENSION; y++)
    {
      for (int x = 0; x < DIMENSION; x++)
      {
        const vm::ivec2 min = newNode.min + vm::ivec2{x, y} * newNode.lod;
        const int lod = newNode.lod;

        addReplacingNode<r>(replacingNodes, min, lod);
      }
    }
    break;
  }
  default:
    std::cerr << "unknown replacing type" << std::endl;
    break;
  }
}

void eraseDataRequest(DataRequestUpdate &dataRequestUpdate, auto &iter)
{
  // forget the data request
  auto currentIter = iter;
  auto nextIter = iter;
  nextIter++;
  dataRequestUpdate.dataRequests.erase(currentIter);
  iter = nextIter;
}

template <bool isReplacing>
void cancelDataRequest(DataRequestUpdate &dataRequestUpdate, const DataRequestPtr &oldDataRequest, auto &iter)
{
  if (!isReplacing)
  {
    dataRequestUpdate.cancelDataRequests.push_back(oldDataRequest);
  }
  eraseDataRequest(dataRequestUpdate, iter);
}

void handleMatchingLodNodes(DataRequestUpdate &dataRequestUpdate,
                            std::unordered_map<uint64_t, OctreeNode> &replacingNodes,
                            const OctreeNode &oldNode,
                            const OctreeNode &newNode,
                            auto &iter)
{
  if (equalsLodArray(newNode, oldNode))
  {
    // keep the data request
    iter++;
  }
  else
  {
    // replacing due to lod array change
    handleReplacingRequest<REPLACING_SIZE::SAME>(dataRequestUpdate, replacingNodes, oldNode, newNode);
  }
}

void handleUnbalancedLodNodes(DataRequestUpdate &dataRequestUpdate,
                              std::unordered_map<uint64_t, OctreeNode> &replacingNodes,
                              const OctreeNode &oldNode,
                              const OctreeNode &newNode,
                              auto &iter)
{
  const float LOD_BALANCE = 1.f;
  const float lodRatio = (float)newNode.lod / (float)oldNode.lod;

  if (lodRatio > LOD_BALANCE)
  {
    // replacing smaller chunk
    handleReplacingRequest<REPLACING_SIZE::SMALLER>(dataRequestUpdate, replacingNodes, oldNode, newNode);
  }
  else
  {
    // replacing bigger chunk
    handleReplacingRequest<REPLACING_SIZE::BIGGER>(dataRequestUpdate, replacingNodes, oldNode, newNode);
  }
}
void handleMatchingMinNodes(DataRequestUpdate &dataRequestUpdate,
                            std::unordered_map<uint64_t, OctreeNode> &replacingNodes,
                            const OctreeNode &oldNode,
                            const OctreeNode &newNode,
                            auto &iter)
{
  if (equalsLod(newNode, oldNode))
  {
    handleMatchingLodNodes(dataRequestUpdate, replacingNodes, oldNode, newNode, iter);
  }
  else
  {
    handleUnbalancedLodNodes(dataRequestUpdate, replacingNodes, oldNode, newNode, iter);
  }
}

void createDataRequestFromLeafNode(OctreeNodePtr &node, DataRequestUpdate &dataRequestUpdate, std::unordered_map<uint64_t, OctreeNode> &replacingNodes)
{
  const uint64_t hash = hashOctreeMinLod(node->min, node->lod);
  auto dataRequestIter = dataRequestUpdate.dataRequests.find(hash);
  if (dataRequestIter == dataRequestUpdate.dataRequests.end())
  {
    auto foundReplacingNode = replacingNodes.find(hash);
    bool isReplacing = foundReplacingNode != replacingNodes.end();

    injectReplacingInfo(node, foundReplacingNode, isReplacing);
    const DataRequestPtr dataRequest = createNewDataRequest(dataRequestUpdate, node, isReplacing);
    addNewDataRequest(dataRequestUpdate, dataRequest, isReplacing);
  }
}

DataRequestUpdate Tracker::updateDataRequests(
    const std::unordered_map<uint64_t, DataRequestPtr> &dataRequests,
    const std::vector<OctreeNodePtr> &leafNodes)
{
  DataRequestUpdate dataRequestUpdate;
  dataRequestUpdate.dataRequests = dataRequests;

  std::unordered_map<uint64_t, OctreeNode> replacingNodes;

  // cancel old data requests, and add replacing data requests
  for (auto iter = dataRequestUpdate.dataRequests.begin(); iter != dataRequestUpdate.dataRequests.end();)
  {
    const uint64_t &hash = iter->first;
    const DataRequestPtr &oldDataRequest = iter->second;

    auto matchingMinLeafNodeIter = std::find_if(
        leafNodes.begin(),
        leafNodes.end(),
        [&](const OctreeNodePtr &leafNode) -> bool
        {
          return equalsNodeMin(*leafNode, *oldDataRequest->node);
        });

    if (matchingMinLeafNodeIter != leafNodes.end())
    {
      const OctreeNode newNode = **matchingMinLeafNodeIter;
      const OctreeNode oldNode = *oldDataRequest->node;

      cancelDataRequest<true>(dataRequestUpdate, oldDataRequest, iter);
      handleMatchingMinNodes(dataRequestUpdate, replacingNodes, oldNode, newNode, iter);
    }
    else
    {
      cancelDataRequest<false>(dataRequestUpdate, oldDataRequest, iter);
    }
  }

  // add new data requests
  for (auto node : leafNodes)
  {
    createDataRequestFromLeafNode(node, dataRequestUpdate, replacingNodes);
  }

  return dataRequestUpdate;
}

TrackerUpdate Tracker::update(PGInstance *inst, const vm::vec3 &position, int minLod, int maxLod, int lod1Range)
{
  std::lock_guard<std::mutex> lock(mutex);

  // new octrees
  const int chunkSize = this->inst->heightfieldGenerator.getChunkSize();
  vm::ivec2 currentCoord = getCurrentCoord(position, chunkSize); // in chunk space

  std::vector<OctreeNodePtr> octreeLeafNodes = constructOctreeForLeaf(
      currentCoord,
      lod1Range,
      minLod,
      maxLod);

  sortNodes(octreeLeafNodes);

  DataRequestUpdate dataRequestUpdate = updateDataRequests(this->dataRequests, octreeLeafNodes);

  dataRequests = std::move(dataRequestUpdate.dataRequests);

  TrackerUpdate result;

  result.leafNodes = std::move(octreeLeafNodes);
  result.newDataRequests = std::move(dataRequestUpdate.newDataRequests);
  result.replacingRequests = std::move(dataRequestUpdate.replacingRequests);
  result.cancelDataRequests = std::move(dataRequestUpdate.cancelDataRequests);

  return result;
}