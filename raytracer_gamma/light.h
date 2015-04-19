
// Includes
#include <CL/cl.h>

struct Light
{
  // Default ctor
  Light() {};

  Light(const cl_float3 &Pos, const cl_float3 &Dir, const cl_float4 &Colour) :
    pos(Pos), dir(Dir), colour(Colour)
  {
  };

  cl_float3 pos;
  cl_float3 dir;
  cl_float4 colour;

};
// EO Struct