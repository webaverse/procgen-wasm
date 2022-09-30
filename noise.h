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

  float in2DWarp(const vec2 &position);
};

#endif