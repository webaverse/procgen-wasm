#include "noise.h"
#include "vectorMath.h"

Noise::Noise(int s, double frequency, int octaves) : fastNoise(s) {
  fastNoise.SetFrequency(frequency);
  fastNoise.SetFractalOctaves(octaves);
}

Noise::~Noise() {}

double Noise::in2DRaw(double x, double y) {
  return fastNoise.GetSimplexFractal(x, y);
}
double Noise::in3DRaw(double x, double y, double z) {
  return fastNoise.GetSimplexFractal(x, y, z);
}

double Noise::in2D(double x, double y) {
  return vm::clamp((1.0 + in2DRaw(x, y)) / 2.0, 0., 1.);
}
double Noise::in3D(double x, double y, double z) {
  return vm::clamp((1.0 + in3DRaw(x, y, z)) / 2.0, 0., 1.);
}
double Noise::in2DBidirectional(double x, double y) {
  return vm::clamp(in2DRaw(x, y), -1., 1.);
}
double Noise::in3DBidirectional(double x, double y, double z) {
  return vm::clamp(in3DRaw(x, y, z), -1., 1.);
}


UberNoise::UberNoise(int s, double frequency, int octaves) : Noise(s, frequency, octaves){};
UberNoise::~UberNoise() {};
float UberNoise::FBM(const vm::vec2 &position)
{
  const int octaves = 4;
	const float frequency = 0.54f;
	const float lacunarity = 3.f;
	const float persistence = 0.7f;

	const float SCALE = 1.f / 128.f;
	vm::vec2 p = position * SCALE;
	double noise = 0.f;

	float amplitude = 1.f;
	p *= frequency;

	for (int i = 0; i < octaves; i++)
	{
		noise += in2D(position.x, position.y) * amplitude;
		p *= lacunarity;	
		amplitude *= persistence;	
	}

	return noise;
}
float UberNoise::in2DTwist(double x, double y)
{
  vm::vec2 p{(float)x, (float)y};
  vm::vec2 q = vm::vec2{FBM(p + vm::vec2{0.0, 0.0}),
                        FBM(p + vm::vec2{7.2, 2.3})};

  // vm::vec2 r = vm::vec2{FBM(p + vm::vec2{50.f, 50.f} * q + vm::vec2{1.7, 9.2}),
  //                       FBM(p + vm::vec2{50.f, 50.f} * q + vm::vec2{8.3, 2.8})};

  return FBM(p + vm::vec2{40.f, 40.f} * q);
  // return FBM(p);
  // return in2D(p.x, p.y);
}