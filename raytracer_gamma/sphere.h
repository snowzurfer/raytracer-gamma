#ifndef _SPHERE_H
#define _SPHERE_H

// Ray struct header file

#include "vec.h"
#include "material.h"

struct Sphere
{
  Vec pos;
  float radius;
  struct Material material;
};
// EO Struct


#endif