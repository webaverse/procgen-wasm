#ifndef _CACHE_H_
#define _CACHE_H_

#include <array>
#include <unordered_map>
#include <mutex>
#include "../sync.h"
#include "../constants.h"
#include <emscripten.h>

//

class PGInstance;

//

inline int modulo(int x, int N){
    return (x % N + N) % N;
}

//

typedef std::array<uint8_t, 4> MaterialsArray;
typedef std::array<float, 4> MaterialsWeightsArray;

class Heightfield {
public:
    float height;

    float liquidHeight = MIN_WORLD_HEIGHT;
    float liquidFactor = 0.f;

    vm::vec3 normal{0.f, 1.f, 0.f};

    std::array<uint8_t, 4> biomes;
    std::array<uint8_t, 4> biomesWeights;

    // TODO : send these to the shader for blending liquid shaders
    std::array<uint8_t, 4> liquids;
    std::array<uint8_t, 4> liquidsWeights;

    MaterialsArray materials;
    MaterialsWeightsArray materialsWeights;

    float getHeight() const {
      return height;
    }
    float getSlope() const
    {
      return std::max(0.f, 1.f - normal.y);
    }
    static bool acceptIndices(
      const Heightfield &a,
      const Heightfield &b,
      const Heightfield &c
    ) {
      return true;
    }
};

class Waterfield : public Heightfield {
public:
    float getHeight() const {
      return liquidHeight;
    }
    bool acceptIndex() const {
      // return true;
      return liquidFactor > 0.f;
    }
    static bool acceptIndices(
      const Waterfield &a,
      const Waterfield &b,
      const Waterfield &c
    ) {
      return a.acceptIndex() && b.acceptIndex() && c.acceptIndex();
    }
};

class BiomeNoiseField {
public:
    float temperature;
    float humidity;
};

class LiquidNoiseField {
public:
    float ocean;
    float river;
};

//

class MaterialWeightAccumulator {
public:
  MaterialWeightAccumulator() : seen(false), weight(0.f) {};

  void addWeight(const float &w){
    seen = true;
    weight += w;
  }

  bool getSeen(){
    return seen;
  }

  float getWeight(){
    return weight;
  }

private:
  bool seen;
  float weight;
};

#endif // _CACHE_H_