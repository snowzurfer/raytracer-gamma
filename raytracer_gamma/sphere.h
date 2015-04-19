// Ray struct header file

#include <CL/cl.h>

struct Sphere
{
  // Default ctor
  Sphere() {
    origin.s[0] = origin.s[1] = origin.s[3] = 0.f;
    radius = 1.f;
  };

  // Ctor 2
  Sphere(const cl_float3 &Origin, const cl_float &Radius) :
    origin(Origin), radius(Radius) {};


  cl_float3 origin;
  cl_float radius;
};
// EO Struct