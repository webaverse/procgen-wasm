#ifndef NOISE_H
#define NOISE_H

#include "FastNoise.h"
#include "Worley.hpp"

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

#endif