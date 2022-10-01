#ifndef _CACHE_H_
#define _CACHE_H_

#include <array>
#include <unordered_map>
#include <mutex>
#include "sync.h"
#include "constants.h"
#include <emscripten.h>

//

class PGInstance;

//

inline int modulo(int x, int N){
    return (x % N + N) % N;
}

//

struct MaterialCountWeightPair {
  std::pair<int, float> pair;

  MaterialCountWeightPair() : pair(0, 0.f){};

  void addWeight(float weight){
    pair.first += 1;
    pair.second += weight;
  }
};
typedef std::array<uint8_t, 4> MaterialsArray;
typedef std::array<float, 4> MaterialsWeightsArray;

class Heightfield {
public:
    float heightField;

    std::array<unsigned char, 4> biomesVectorField;
    std::array<unsigned char, 4> biomesWeightsVectorField;

    float waterFactor;

    MaterialsArray materials;
    MaterialsWeightsArray materialsWeights;

    float getHeight() const {
      return heightField;
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
      return WORLD_BASE_HEIGHT;
    }
    bool acceptIndex() const {
      return waterFactor > 0;
    }
    static bool acceptIndices(
      const Waterfield &a,
      const Waterfield &b,
      const Waterfield &c
    ) {
      return a.acceptIndex() && b.acceptIndex() && c.acceptIndex();
    }
};

class NoiseField {
public:
    float temperature;
    float humidity;
    float ocean;
    float river;
};

//

template<typename T>
class HashValue {
public:
    uint32_t hash;
    T value;

    HashValue() : hash(0xFFFFFFFF) {}
    HashValue(uint32_t hash, T value) : hash(hash), value(value) {}
};

template <typename T>
class ChunkCacheValue {
public:
    union {
      HashValue<T> hashValue;
      std::atomic<uint64_t> data;
    };

    static_assert(sizeof(hashValue) == sizeof(data), "wrong size");

    ChunkCacheValue() : hashValue() {}
    HashValue<T> get() {
        uint64_t localData = data.load();
        uint64_t *localDataP = &localData;
        HashValue<T> &hashValue = *((HashValue<T> *)localDataP);
        return hashValue;
    }
    void set(uint32_t hash, const T &value) {
      HashValue<T> hashValue{hash, value};
      uint64_t &localData = *((uint64_t *)&hashValue);
      data.store(localData);
    }
};
template <>
class ChunkCacheValue<Heightfield> {
public:
  union {
    HashValue<float> hashValue1;
    std::atomic<uint64_t> data1;
  };
  union {
    HashValue<std::array<unsigned char, 4>> hashValue2;
    std::atomic<uint64_t> data2;
  };
  union {
    HashValue<std::array<unsigned char, 4>> hashValue3;
    std::atomic<uint64_t> data3;
  };

  static_assert(sizeof(hashValue1) == sizeof(data1), "wrong size");
  static_assert(sizeof(hashValue2) == sizeof(data2), "wrong size");
  static_assert(sizeof(hashValue3) == sizeof(data3), "wrong size");

  ChunkCacheValue() : hashValue1(), hashValue2(), hashValue3() {}
  HashValue<Heightfield> get() {
      uint64_t localData1 = data1.load();
      uint64_t localData2 = data2.load();
      uint64_t localData3 = data3.load();

      uint64_t *localDataP1 = &localData1;
      uint64_t *localDataP2 = &localData2;
      uint64_t *localDataP3 = &localData3;

      HashValue<float> &hashValue1 = *((HashValue<float> *)localDataP1);
      HashValue<std::array<unsigned char, 4>> &hashValue2 = *((HashValue<std::array<unsigned char, 4>> *)localDataP2);
      HashValue<std::array<unsigned char, 4>> &hashValue3 = *((HashValue<std::array<unsigned char, 4>> *)localDataP3);

      if (hashValue1.hash == hashValue2.hash && hashValue1.hash == hashValue3.hash) {
        return HashValue<Heightfield>{
          hashValue1.hash,
          Heightfield{hashValue1.value, hashValue2.value, hashValue3.value}
        };
      } else {
        return HashValue<Heightfield>{
          0,
          Heightfield()
        };
      }
  }
  void set(uint32_t hash, const Heightfield &value) {
    HashValue<float> localHashValue1{hash, value.heightField};
    HashValue<std::array<unsigned char, 4>> localHashValue2{hash, value.biomesVectorField};
    HashValue<std::array<unsigned char, 4>> localHashValue3{hash, value.biomesWeightsVectorField};

    uint64_t &localData1 = *((uint64_t *)&localHashValue1);
    uint64_t &localData2 = *((uint64_t *)&localHashValue2);
    uint64_t &localData3 = *((uint64_t *)&localHashValue3);

    data1.store(localData1);
    data2.store(localData2);
    data3.store(localData3);
  }
};

template<>
class ChunkCacheValue<NoiseField> {
public:
    union {
        HashValue<float> hashValue1;
        std::atomic<uint64_t> data1;
    };
    union {
        HashValue<float> hashValue2;
        std::atomic<uint64_t> data2;
    };
    union {
        HashValue<float> hashValue3;
        std::atomic<uint64_t> data3;
    };
    union {
        HashValue<float> hashValue4;
        std::atomic<uint64_t> data4;
    };

    static_assert(sizeof(hashValue1) == sizeof(data1), "wrong size");
    static_assert(sizeof(hashValue2) == sizeof(data2), "wrong size");
    static_assert(sizeof(hashValue3) == sizeof(data3), "wrong size");
    static_assert(sizeof(hashValue4) == sizeof(data4), "wrong size");

    ChunkCacheValue() : hashValue1(), hashValue2(), hashValue3(), hashValue4() {}
    HashValue<NoiseField> get() {
        uint64_t localData1 = data1.load();
        uint64_t localData2 = data2.load();
        uint64_t localData3 = data3.load();
        uint64_t localData4 = data4.load();

        uint64_t *localDataP1 = &localData1;
        uint64_t *localDataP2 = &localData2;
        uint64_t *localDataP3 = &localData3;
        uint64_t *localDataP4 = &localData4;

        HashValue<float> &hashValue1 = *((HashValue<float> *)localDataP1);
        HashValue<float> &hashValue2 = *((HashValue<float> *)localDataP2);
        HashValue<float> &hashValue3 = *((HashValue<float> *)localDataP3);
        HashValue<float> &hashValue4 = *((HashValue<float> *)localDataP4);

        if (hashValue1.hash == hashValue2.hash && hashValue1.hash == hashValue3.hash && hashValue1.hash == hashValue4.hash) {
            return HashValue<NoiseField>{
                hashValue1.hash,
                NoiseField{hashValue1.value, hashValue2.value, hashValue3.value, hashValue4.value}
            };
        } else {
            return HashValue<NoiseField>{
                0,
                NoiseField()
            };
        }
    }
    void set(uint32_t hash, const NoiseField &value) {
        HashValue<float> localHashValue1{hash, value.temperature};
        HashValue<float> localHashValue2{hash, value.humidity};
        HashValue<float> localHashValue3{hash, value.ocean};
        HashValue<float> localHashValue4{hash, value.river};

        uint64_t &localData1 = *((uint64_t *)&localHashValue1);
        uint64_t &localData2 = *((uint64_t *)&localHashValue2);
        uint64_t &localData3 = *((uint64_t *)&localHashValue3);
        uint64_t &localData4 = *((uint64_t *)&localHashValue4);

        data1.store(localData1);
        data2.store(localData2);
        data3.store(localData3);
        data4.store(localData4);
    }
};
//

uint32_t getCacheIndexWorld(int x, int y);
uint32_t getCacheIndexWorld(int x, int y, int z);

//

uint32_t getCacheIndexLocal(int x, int y);
uint32_t getCacheIndexLocal(int x, int y, int z);

#endif // _CACHE_H_