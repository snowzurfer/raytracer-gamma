#ifndef _INTERSECTION_H
#define _INTERSECTION_H

#include "vec.h"
#include "sphere.h"

struct Intersection {
  // Object which intersects with the ray
  struct Sphere object;

  // Point of intersection
  Vec point;
  // Normal at intersection
  Vec normal;

  // Distance non-sqrted
  float squaredDist;
};

#endif