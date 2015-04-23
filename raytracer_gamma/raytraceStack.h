#ifndef _RAYTRACE_STACK_H
#define _RAYTRACE_STACK_H

// Includes
#include "vec.h"
#include "ray.h"


struct rtSnapshot 
{
  struct Ray ray;
  int traceDepth;
  int stage;
};

#endif