// Ray struct header file

#include <CL/cl.h>

namespace rtg {

  struct Ray
  {
    // Default ctor
    Ray() {
      origin.s[0] = origin.s[1] = origin.s[3] = 0.f;
      dir.s[0] = dir.s[1] = dir.s[3] = 0.f;
    };

    // Ctor 2
    Ray(cl_float3 &Origin, cl_float3 &Direction) :
      origin(Origin), dir(Direction) {};


    cl_float3 origin;
    cl_float3 dir;
  };
  // EO Struct

}
// EO Namespace