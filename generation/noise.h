#ifndef NOISE_H
#define NOISE_H

#include "../libs/FastNoise.h"
#include "../libs/Worley.hpp"
#include "../libs/vectorMath.h"
#include "../biomes.h"
#include "heightfield.h"
#include <iostream>

using namespace vm;

bool getNoiseVisibility(float value, float min, float max = 1.f);

class Noise
{
public:
  FastNoise fastNoise;

  explicit Noise(int s = 0, double frequency = 0.01, int octaves = 1);
  ~Noise();

  double in2DRaw(double x, double y);
  double in3DRaw(double x, double y, double z);

  double in2D(double x, double y);
  double in3D(double x, double y, double z);

  double in2DBidirectional(double x, double y);
  double in3DBidirectional(double x, double y, double z);
};

class UberNoise
{
public:
  UberNoise();
  ~UberNoise();

  float humidityNoise(float x, float z);
  float temperatureNoise(float x, float z);
  float grassMaterialNoise(float x, float z, float wetness);
  float wetnessNoise(float x, float z);
  float coldnessNoise(float x, float z);
  float desertNoise(float x, float z);
  float mountainNoise(float x, float z);
  float iceMountainNoise(float x, float z);
  float hashNoise(float x, float z);
  float simplexNoise(float x, float z);
  float waterDepthNoise(float x, float z);
  float oceanNoise(float x, float z);
  float riverNoise(float x, float z, float ocean);
  bool waterVisibilityNoise(float x, float z);

  vm::vec3 flowNoise(float x, float z);

  template <uint8_t T>
  bool instanceVisibility(float x, float z, const Heightfield &heightfield)
  {
    switch (T)
    {
    case (uint8_t)INSTANCE::TREE:
    {
      return heightfield.field.tree;
    }
    case (uint8_t)INSTANCE::FLOWER:
    {
      return heightfield.field.flower;
    }
    case (uint8_t)INSTANCE::GRASS:
    {
      const float grass = heightfield.field.grass;
      return getNoiseVisibility(grass, GRASS_THRESHOLD);
    }
    case (uint8_t)INSTANCE::ROCK:
      return heightfield.field.rock;
    case (uint8_t)INSTANCE::STONE:
      return heightfield.field.stone;
    default:
      return false;
    }
  }

  template <uint8_t T>
  float scaleNoise(float x, float z, const Heightfield &heightfield)
  {
    const float hash = (hashNoise(x, z) * 2.0 - 1.0);

    switch (T)
    {
    case (uint8_t)INSTANCE::TREE:
    {
      const float wetness = heightfield.field.wetness;
      return TREE_BASE_SCALE + wetness + TREE_SCALE_RANGE * hash;
    }
    case (uint8_t)INSTANCE::FLOWER:
    {
      const float wetness = heightfield.field.wetness;
      return FLOWER_BASE_SCALE + wetness + FLOWER_SCALE_RANGE * hash;
    }
    case (uint8_t)INSTANCE::GRASS:
    {
      return GRASS_BASE_SCALE + GRASS_SCALE_RANGE * hash;
    }
    case (uint8_t)INSTANCE::ROCK:
    {
      const float dryness = 1.0 - heightfield.field.wetness;
      return ROCK_BASE_SCALE + dryness + ROCK_SCALE_RANGE * hash;
    }
    case (uint8_t)INSTANCE::STONE:
    {
      const float dryness = 1.0 - heightfield.field.wetness;
      return STONE_BASE_SCALE + dryness + STONE_SCALE_RANGE * hash;
    }
    default:
    {
      return 0.f;
    }
    }
  }

  template <uint8_t T>
  vm::vec3 rotationNoise(float x, float z, const Heightfield &heightfield)
  {
    vm::vec3 rotation;

    const float TAU = 2.0f * M_PI;
    const float hash1 = hashNoise(x, z) * 2.0f - 1.0f;
    const float hash2 = hashNoise(x + 20.f, z - 30.f) * 2.0f - 1.0f;
    const float hash3 = hashNoise(x - 10.f, z + 5.f) * 2.0f - 1.0f;

    switch (T)
    {
    case (uint8_t)INSTANCE::TREE:
    {
      rotation = vm::vec3{0.f, hash1, 0.f};
      break;
    }
    case (uint8_t)INSTANCE::FLOWER:
    {
      rotation = vm::vec3{hash1 / 1.5f, hash2, hash3 / 1.3f} / TAU;
      break;
    }
    case (uint8_t)INSTANCE::GRASS:
    {
      rotation = vm::vec3{hash3 / 1.7f, hash1 / 1.5f, hash3 / 1.6f} / TAU;
      break;
    }
    case (uint8_t)INSTANCE::ROCK:
    {
      rotation = vm::vec3{hash1 / 1.5f, hash3, hash2 / 1.3f};
      break;
    }
    case (uint8_t)INSTANCE::STONE:
    {
      rotation = vm::vec3{hash2 / 1.1f, hash3 / 1.5f, hash1 / 1.5f};
      break;
    }
    default:
    {
      rotation = vm::vec3{0.f, hash1, 0.f};
      break;
    }
    }

    return rotation * TAU;
  }

  template <uint8_t T>
  vm::vec3 colorNoise(float x, float z, const Heightfield &heightfield)
  {
    vm::vec3 color;

    const float simplexm15 = simplexNoise(x * 15.f ,z * 15.f ) * 2.0f - 1.0f;

    switch (T)
    {
    case (uint8_t)INSTANCE::TREE:
    {
      color = vm::vec3{0.f, 0.f, 0.f};
      break;
    }
    case (uint8_t)INSTANCE::FLOWER:
    {
      color = vm::vec3{0.f, 0.f, 0.f};
      break;
    }
    case (uint8_t)INSTANCE::GRASS:
    {
      const float colorVariationNoise = vm::clamp(GRASS_COLOR_VARIATION_BASE + simplexm15 * GRASS_COLOR_VARIATION_RANGE, 0.0f, 1.f + GRASS_COLOR_VARIATION_RANGE);
      const float randomBladeFactor = heightfield.field.hash * 2.f - 1.f;
      color = (vm::vec3{colorVariationNoise / 1.4f, colorVariationNoise / 1.6f, colorVariationNoise / 1.75f} +
               vm::vec3{randomBladeFactor / 8.f, randomBladeFactor / 10.f, randomBladeFactor / 12.f});
      break;
    }
    case (uint8_t)INSTANCE::ROCK:
    {
      color = vm::vec3{0.f, 0.f, 0.f};
      break;
    }
    case (uint8_t)INSTANCE::STONE:
    {
      color = vm::vec3{0.f, 0.f, 0.f};
      break;
    }
    default:
    {
      color = vm::vec3{0.f, 0.f, 0.f};
      break;
    }
    }

    return color;
  }

  bool flowerVisibility(float x, float z, float grass);
  bool grassVisibility(float x, float z, float wetness);
  bool treeVisibility(float x, float z, float wetness);
  bool rockVisibility(float x, float z);
  bool stoneVisibility(float x, float z);
};

#endif