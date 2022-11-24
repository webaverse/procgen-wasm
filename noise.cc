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
float UberNoise::grassMaterialNoise(float x, float z)
{
  return GLSL::grassMaterialNoise(vm::vec2(x, z));
}
bool UberNoise::grassVisibility(float x, float z)
{
  return GLSL::grassVisibility(vm::vec2(x, z));
}
bool UberNoise::treeVisibility(float x, float z)
{
  return GLSL::treeVisibility(vm::vec2(x, z));
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
template<typename T>
bool UberNoise::instanceVisibility(float x, float z)
{
  switch (T)
  {
  case INSTANCE::TREE:
    return treeVisibility(x, z); 
  case INSTANCE::ROCK:
    return rockVisibility(x, z); 
  case INSTANCE::GRASS:
    return grassVisibility(x, z); 
  default:
    std::cerr << "Unknown instance type" << std::endl;
    break;
  }
}
bool UberNoise::rockVisibility(float x, float z)
{
  return GLSL::rockVisibility(vm::vec2(x, z));
}