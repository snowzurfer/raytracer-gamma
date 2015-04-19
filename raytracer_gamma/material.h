
// Includes
#include <CL/cl.h>

struct Material
{
  // Default ctor
  Material() {};

  // Ctor 1
  Material(const cl_float4 &Colour, const cl_float &Reflectivity,
    const cl_float &Refractivity) :
    colour(Colour), reflectivity(Reflectivity),
    refractivity(Refractivity) 
  {
  };


  cl_float4 colour;
  cl_float reflectivity;
  cl_float refractivity;
};
// EO Struct