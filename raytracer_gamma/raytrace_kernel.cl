#define GPU_KERNEL

#include "vec.h"
#include "raytracer.h"

__kernel void raytrace( 
   __global struct Sphere *spheres, 
   __private const unsigned int sphNum,       
   __global struct Light *lights,
   __private const unsigned int lgtNum,
   __private const unsigned int width,
   __private const unsigned int height,
   __global float4 *dst)     
{                                                                     
   int i = get_global_id(0);                                          
   
	
} 