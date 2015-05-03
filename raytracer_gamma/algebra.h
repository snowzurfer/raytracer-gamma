
#ifndef _ALGEBRA_H
#define _ALGEBRA_H

#include <vec.h>


#ifndef GPU_KERNEL
#include <math.h>
#endif

// Tolerance used in floating point comparisons
#define TOL    (0.001f)   

bool isZero(float x);

// Returns the number of roots found and if the equation
// cannot be solved returns zero.
//
// It also solves linear equations.
//
// The result is stored in *roots
int solveQuadratic(float a, float b, float c, float *roots);

// Computes the maximum colour value for the given pixel buffer
float maxColourValuePixelBuffer(const Vec *bfr, size_t n);

#endif