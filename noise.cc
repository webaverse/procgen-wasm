#include "noise.h"
#include "vectorMath.h"
#include "glsl.h"

using namespace vm;

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
  return clamp((1.0 + in2DRaw(x, y)) / 2.0, 0., 1.);
}
double Noise::in3D(double x, double y, double z) {
  return clamp((1.0 + in3DRaw(x, y, z)) / 2.0, 0., 1.);
}
double Noise::in2DBidirectional(double x, double y) {
  return clamp(in2DRaw(x, y), -1., 1.);
}
double Noise::in3DBidirectional(double x, double y, double z) {
  return clamp(in3DRaw(x, y, z), -1., 1.);
}


UberNoise::UberNoise(){};
UberNoise::~UberNoise() {};

float UberNoise::humidityNoise(float x, float z)
{
  return GLSL::humidityNoise(vm::vec2(x, z));
}
float UberNoise::temperatureNoise(float x, float z)
{
  return GLSL::temperatureNoise(vm::vec2(x, z));
}
float UberNoise::wetnessNoise(float x, float z)
{
  return GLSL::wetnessNoise(vm::vec2(x, z));
}
float UberNoise::stiffnessNoise(float x, float z)
{
  return GLSL::stiffnessNoise(vm::vec2(x, z));  
}
float UberNoise::desertNoise(float x, float z)
{
  return GLSL::desertNoise(vm::vec2(x, z));
}
float UberNoise::mountainNoise(float x, float z)
{
  return GLSL::mountainNoise(vm::vec2(x, z));
}
float UberNoise::iceMountainNoise(float x, float z)
{
  return GLSL::iceMountainNoise(vm::vec2(x, z));
}
float UberNoise::oceanNoise(float x, float z)
{
  return GLSL::oceanNoise(vm::vec2(x, z));
}