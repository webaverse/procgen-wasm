#include "noise.h"
#include "../libs/vectorMath.h"
#include "glsl.h"

using namespace vm;

bool getNoiseVisibility(float value, float min, float max)
{
    return value >= min && value <= max;
}

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
float UberNoise::grassMaterialNoise(float x, float z, float wetness)
{
  return GLSL::grassMaterialNoise(vm::vec2(x, z), wetness);
}
float UberNoise::wetnessNoise(float x, float z)
{
  return GLSL::wetnessNoise(vm::vec2(x, z));
}
bool UberNoise::flowerVisibility(float x, float z, float grass)
{
  return GLSL::flowerVisibility(vm::vec2(x, z), grass);
}
bool UberNoise::grassVisibility(float x, float z, float wetness)
{
  return GLSL::grassVisibility(vm::vec2(x, z), wetness);
}
bool UberNoise::treeVisibility(float x, float z, float wetness)
{
  return GLSL::treeVisibility(vm::vec2(x, z), wetness);
}
float UberNoise::coldnessNoise(float x, float z)
{
  return GLSL::coldnessNoise(vm::vec2(x, z));  
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
float UberNoise::hashNoise(float x, float z)
{
  return GLSL::hashNoise(vm::vec2(x, z));
}
float UberNoise::simplexNoise(float x, float z)
{
  return GLSL::simplexNoise(vm::vec2(x, z));
}
float UberNoise::waterDepthNoise(float x, float z)
{
  return GLSL::waterDepthNoise(vm::vec2(x, z));
}
float UberNoise::oceanNoise(float x, float z)
{
  return GLSL::oceanNoise(vm::vec2(x, z));
}
float UberNoise::riverNoise(float x, float z, float ocean)
{
  return GLSL::riverNoise(vm::vec2(x, z), ocean);
}
bool UberNoise::waterVisibilityNoise(float x, float z)
{
  return GLSL::waterVisibilityNoise(vm::vec2(x, z));
}
bool UberNoise::rockVisibility(float x, float z)
{
  return GLSL::rockVisibility(vm::vec2(x, z));
}
bool UberNoise::stoneVisibility(float x, float z)
{
  return GLSL::stoneVisibility(vm::vec2(x, z));
}

vm::vec3 UberNoise::flowNoise(float x, float z)
{
  return GLSL::flowNoise(vm::vec2(x, z));
}