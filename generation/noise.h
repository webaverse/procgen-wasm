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
  float stiffnessNoise(float x, float z);
  float desertNoise(float x, float z);
  float mountainNoise(float x, float z);
  float iceMountainNoise(float x, float z);
  float hashNoise(float x, float z);
  float simplexNoise(float x, float z);
  float waterDepthNoise(float x, float z);
  float oceanNoise(float x, float z);
  float riverNoise(float x, float z, float ocean);
  bool waterVisibilityNoise(float x, float z);

  template <uint8_t T>
  bool instanceVisibility(float x, float z, const Heightfield &heightfield)
  {
    switch (T)
    {
    case (uint8_t)INSTANCE::TREE:
    {
      const float wetness = heightfield.field.wetness;
      return treeVisibility(x, z, wetness);
    }
    case (uint8_t)INSTANCE::FLOWER:
      {
        const float grass = heightfield.field.grass;
        return flowerVisibility(x, z, grass);
      }
    case (uint8_t)INSTANCE::GRASS:
    {
      const float grass = heightfield.field.grass;
      return getNoiseVisibility(grass, GRASS_THRESHOLD);
    }
    case (uint8_t)INSTANCE::ROCK:
      return rockVisibility(x, z);
    case (uint8_t)INSTANCE::STONE:
      return stoneVisibility(x, z);
    default:
      std::cerr << "Unknown instance type" << std::endl;
      break;
    }

    return false;
  }

  bool flowerVisibility(float x, float z, float grass);
  bool grassVisibility(float x, float z, float wetness);
  bool treeVisibility(float x, float z, float wetness);
  bool rockVisibility(float x, float z);
  bool stoneVisibility(float x, float z);
};

#endif