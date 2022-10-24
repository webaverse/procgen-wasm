#ifndef NOISE_H
#define NOISE_H

#include "FastNoise.h"
#include "Worley.hpp"
#include "vectorMath.h"

using namespace vm;

class Noise {
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
  float grassMaterialNoise(float x, float z);
  float grassObjectNoise(float x, float z);
  float treeObjectNoise(float x, float z);
  float stiffnessNoise(float x, float z);
  float desertNoise(float x, float z);
  float mountainNoise(float x, float z);
  float iceMountainNoise(float x, float z);
  float oceanNoise(float x, float z);
};

#endif