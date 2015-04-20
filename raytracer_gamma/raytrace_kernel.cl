#define GPU_KERNEL

#include "vec.h"
#include "raytracer.h"

__kernel void raytrace( 
	__global struct Sphere *spheres, 
	__private const unsigned int kSphNum,       
	__global struct Light *lights,
	__private const unsigned int kLgtNum,
	__private const unsigned int kWidth,
	__private const unsigned int kHeight,
	__private const float kZoom,
	__private const float kAliasFactor,
	__global float4 *dst)     
{                            
	// Retrieve the global ID of the kernel
	const unsigned gid = get_global_id(0);                                          
   
	// Calculate inverse of aliasFactor
	const float kAliasFactorInv = 1.f / kAliasFactor;
	// Calculate total size of samples to be taken
	const float kSamplesTot = kAliasFactor * kAliasFactor;
	// Also its inverse
	const float kSamplesTotinv = 1.f / kSamplesTot;
	
	// Calculate world position of pixel being currently worked on
	const float kPxWorldX = (((float)gid) - (kWidth * 0.5f));
	const float kPxWorldY = ( (kHeight *0.5f) - ((float)gid));
	
	// The ray to be shot. The vantage point (camera) is at the origin,
	// and its intensity is maximum
	struct Ray ray; vinit(ray.origin, 0.f, 0.f, 0.f); vinit(ray.intensity, 1.f, 1.f, 1.f);
	
	// The colour of the pixel to be computed
	float4 pixelCol = (float4)(0.f);
	
	// Mock background material
	struct Material bgMaterial;
	Vec black; vinit(black, 0.f, 0.f, 0.f);
	setMatteGlossBalance(&bgMaterial, 0.f, &black, &black);
	
	// For each sample to be taken
	/*for(int i = 0; i < kAliasFactor; ++i) {
		for(int j = 0; j < kAliasFactor; ++j) {
			// Calculate the direction of the ray
			float x = kPxWorldX + (float)(((float)j) * kAliasFactorInv);
			float y = kPxWorldY + (float)(((float)i) * kAliasFactorInv);
			
			// Set the ray's dir and normalise it
			vinit(ray.dir, x, y, kZoom); vnorm(ray.dir);
			
			// Raytrace for the current sample
			Vec currentSampleCol = rayTrace(spheres, kSphNum, lights, kLgtNum,
				&ray, &bgMaterial, 0);
			
			float4 float4Samplecol; 
			float4Samplecol.s[0] = currentSampleCol.x;
			float4Samplecol.s[1] = currentSampleCol.y;
			float4Samplecol.s[2] = currentSampleCol.z;
			float4Samplecol.s[3] = 1.f;
			
			// Compute the average
			pixelCol += (float4Samplecol * kSamplesTotinv);
		}
	}*/
	
	// Write result in destination buffer
	dst[gid] = pixelCol;
} 