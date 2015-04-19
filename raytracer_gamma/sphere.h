// Ray struct header file

#include <vec.h>
#include <material.h>
#include <ray.h>
#include <common.h>
#include <stdbool.h>

struct Sphere
{
  Vec origin;
  cl_float radius;
  struct Material material;
};
// EO Struct



bool raySphere(struct Sphere * sphere, struct Ray *ray, cl_float *t);