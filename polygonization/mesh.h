#ifndef MESH_H
#define MESH_H

#include "biomes.h"
#include "libs/vectorMath.h"
#include "generation/heightfield.h"

#include <vector>
#include <array>
#include <map>
#include <stdint.h>
#include <cstring>

//

class OctreeNode;
typedef std::shared_ptr<OctreeNode> OctreeNodePtr;

//

typedef std::vector<vm::vec3> PositionBuffer;
typedef std::vector<vm::ivec2> PositionInt2DBuffer;
typedef std::vector<vm::vec3> NormalBuffer;
typedef std::vector<vm::vec2> UvBuffer;
typedef std::vector<vm::ivec4> BiomesBuffer;
typedef std::vector<vm::vec4> BiomesWeightBuffer;
typedef std::vector<vm::ivec4> MaterialsBuffer;
typedef std::vector<vm::vec4> MaterialsWeightBuffer;
typedef std::vector<vm::vec4> MaterialsSplatBuffer;
typedef std::vector<std::array<UV, 2>> BiomesUvsBuffer;
typedef std::vector<uint32_t> IndexBuffer;
typedef std::vector<int> BiomeBuffer;
// typedef std::vector<float> SeedBuffer;
typedef std::vector<float> FactorBuffer;
typedef std::vector<uint8_t> LightBuffer;
typedef unsigned char PeekBuffer[15];

//

class TerrainGeometry {
public:
    PositionBuffer positions;
    NormalBuffer normals;
    IndexBuffer indices;

    BiomesBuffer biomes;
    BiomesWeightBuffer biomesWeights;
    BiomesUvsBuffer biomesUvs1;
    BiomesUvsBuffer biomesUvs2;

    MaterialsBuffer materials;
    MaterialsWeightBuffer materialsWeights;

    // SeedBuffer seeds;

    // LightBuffer skylights;
    // LightBuffer aos;

    // PeekBuffer peeks;

    void pushPointMetadata(const Heightfield &fieldValue);

    uint8_t *getBuffer() const;
};

//

class WaterGeometry {
public:
    PositionBuffer positions;
    NormalBuffer normals;
    IndexBuffer indices;

    FactorBuffer factors;

    MaterialsBuffer materials;
    MaterialsWeightBuffer materialsWeights;
    // PeekBuffer peeks;

    void pushPointMetadata(const Waterfield &fieldValue);

    uint8_t *getBuffer() const;
};

//

class SplatInstance {
public:
    int instanceId;
    std::vector<float> ps;
    std::vector<float> qs;
};
class SplatInstanceGeometry {
public:
    std::map<int, SplatInstance> instances;

    uint8_t *getBuffer() const;
};

class VegetationGeometry : public SplatInstanceGeometry {
    // nothing
};

class RockGeometry : public SplatInstanceGeometry {
    // nothing
};
//

class MaterialAwareSplatInstance : public SplatInstance {
public:
    // terrain materials
    MaterialsSplatBuffer materials; 
    MaterialsSplatBuffer materialsWeights; 
};

class GrassSplatInstance : public MaterialAwareSplatInstance {
public:
    MaterialsSplatBuffer grassProps; 
};

class MaterialAwareSplatInstanceGeometry {
public:
    std::map<int, MaterialAwareSplatInstance> instances;

    uint8_t *getBuffer() const;
};

class GrassGeometry {
public:
    std::map<int, GrassSplatInstance> instances;

    uint8_t *getBuffer() const;
};

class PoiGeometry {
public:
    std::vector<float> ps;
    std::vector<int32_t> instances;

    uint8_t *getBuffer() const;
};

//

class BarrierGeometry {
public:
    PositionBuffer positions;
    NormalBuffer normals;
    IndexBuffer indices;
    
    UvBuffer uvs;
    PositionInt2DBuffer positions2D;

    std::vector<OctreeNodePtr> leafNodes;
    vm::ivec2 leafNodesMin;
    vm::ivec2 leafNodesMax;
    std::vector<int> leafNodesIndex;

    uint8_t *getBuffer() const;
};

//

class HeightfieldGeometry {
public:
    std::vector<vm::vec4> heightfieldImage;

    uint8_t *getBuffer() const;
};

#endif // MESH_H