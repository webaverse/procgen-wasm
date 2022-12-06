#ifndef _CACHE_H_
#define _CACHE_H_

#include <array>
#include <unordered_map>
#include <mutex>
#include "../task/sync.h"
#include "../constants.h"
#include "../utils/util.h"
#include <emscripten.h>

//

class PGInstance;

//

inline int modulo(int x, int N)
{
  return (x % N + N) % N;
}

//

typedef std::array<uint8_t, 4> BiomesArray;
typedef std::array<float, 4> BiomesWeightsArray;
typedef std::array<uint8_t, 4> MaterialsArray;
typedef std::array<float, 4> MaterialsWeightsArray;
typedef std::array<uint8_t, 4> LiquidsArray;
typedef std::array<float, 4> LiquidsWeightsArray;

struct CachedField {
  float height = 0.0f;

  float liquidHeight = MIN_WORLD_HEIGHT;
  float liquidFactor = 0.f;

  float wetness = 0.f;
  float grass = 0.f;
};

class Heightfield
{
public:
  CachedField field;

  vm::vec3 normal{0.f, 1.f, 0.f};

  BiomesArray biomes;
  BiomesWeightsArray biomesWeights;

  LiquidsArray liquids;
  LiquidsWeightsArray liquidsWeights;

  MaterialsArray materials;
  MaterialsWeightsArray materialsWeights;

  float getHeight() const
  {
    return field.height;
  }
  bool hasWater() const
  {
    return field.liquidFactor > 0.f;
  }
  float getSlope() const
  {
    return std::max(0.f, 1.f - normal.y);
  }
  uint8_t getDominantBiome() const
  {
    // because biomes are sorted by weight the dominant biome is the first one
    return biomes[0];
  }
  static bool acceptIndices(
      const Heightfield &a,
      const Heightfield &b,
      const Heightfield &c)
  {
    return true;
  }
};

class Waterfield : public Heightfield
{
public:
  float getHeight() const
  {
    return field.liquidHeight;
  }
  bool acceptIndex() const
  {
    return field.liquidFactor > 0.f;
  }
  static bool acceptIndices(
      const Waterfield &a,
      const Waterfield &b,
      const Waterfield &c)
  {
    return a.acceptIndex() && b.acceptIndex() && c.acceptIndex();
  }
};

class HeightfieldSampler
{
public:
  const vm::vec2 &worldPositionXZ;
  const int &lod;
  const int chunkSize;
  const int chunkSizeP2;
  const std::vector<Heightfield> &heightfields;

  HeightfieldSampler(
      const vm::vec2 &worldPositionXZ,
      const int &lod,
      const int chunkSize,
      const std::vector<Heightfield> &heightfields) : worldPositionXZ(worldPositionXZ),
                                                      lod(lod),
                                                      chunkSize(chunkSize),
                                                      chunkSizeP2(chunkSize + 2),
                                                      heightfields(heightfields)
  {
  }

  float getHeight(float x, float z)
  {
    const vm::vec2 location = getLocalPosition(x, z);

    float result = bilinear<HeightfieldSampler, float>(location, chunkSize, *this);
    return result;
  }

  float getWorldHeight(float x, float z)
  {
    return getHeight(x, z) - (float)WORLD_BASE_HEIGHT;
  }

  Heightfield getHeightfield(float x, float z)
  {
    const vm::vec2 location = getLocalPosition(x, z);

    float rx = std::floor(location.x);
    float ry = std::floor(location.y);

    int ix = (int)rx;
    int iy = (int)ry;

    return getHeightfieldByLocalPosition(ix, iy);
  }

  float get(int x, int z)
  {
    const Heightfield heightfield = getHeightfieldByLocalPosition(x, z);
    const float &height = heightfield.field.height;
    return height;
  }

private:
  vm::vec2 getLocalPosition(float x, float z)
  {
    vm::vec2 location{
        x - worldPositionXZ.x,
        z - worldPositionXZ.y};
    location /= (float)lod;
    return location;
  };
  Heightfield getHeightfieldByLocalPosition(int x, int z)
  {
    const int dx = x + 1;
    const int dz = z + 1;
    const int index = dx + dz * chunkSizeP2;
    const Heightfield &heightfield = heightfields[index];
    return heightfield;
  }
};

class BiomeNoiseField
{
public:
  float temperature;
  float humidity;
};

class LiquidNoiseField
{
public:
  float ocean;
  float river;
};

//

class MaterialWeightAccumulator
{
public:
  MaterialWeightAccumulator() : seen(false), weight(0.f){};

  void addWeight(const float &w)
  {
    seen = true;
    weight += w;
  }

  bool getSeen()
  {
    return seen;
  }

  float getWeight()
  {
    return weight;
  }

private:
  bool seen;
  float weight;
};

#endif // _CACHE_H_