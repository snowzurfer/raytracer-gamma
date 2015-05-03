#include <algebra.h>

bool isZero(float x) {
  return (fabs(x) < TOL);
}

// Returns the number of roots found and if the equation
// cannot be solved returns zero.
//
// It also solves linear equations.
//
// The result is stored in *roots
int solveQuadratic(float a, float b, float c, float *roots) {
  // If the a component is zero
  if (isZero(a))  {
    // If also the b component is zero
    if (isZero(b)) {
      // The equation devolves to: c = 0, where the variable x has vanished!
      // Cannot divide by zero, so there is no solution.
      return 0;
    }
    else {
      // Simple linear equation: bx + c = 0, so x = -c/b.
      roots[0] = -c / b;

      // There is a single solution.
      return 1;
    }
  }
  // The equation is quadratic
  else {
    // Calculate the radicand of the quadratic equation solution formula.
    const float radicand = (b*b) - (4.f * a * c);

    // If the radicand is zero
    if (isZero(radicand)) {
      // Both roots have the same value: -b / 2a.
      roots[0] = -b / (2.f * a);

      // Return the number of solutions
      return 1;
    }
    else {
      // There are two distinct real roots.
      const float root = sqrt(radicand);
      const float denom = 2.0f * a;

      // Store them in the output variable
      roots[0] = (-b + root) / denom;
      roots[1] = (-b - root) / denom;

      // Return the number of solutions
      return 2;
    }
  }
}

// Computes the maximum colour value for the given pixel buffer
float maxColourValuePixelBuffer(const Vec *bfr, size_t n) {
  // Max value to be returned
  float max = 0.f;

  // For each pixel
  for (size_t i = 0; i < n; ++i) {
    if (bfr[i].x > max) {
      max = bfr[i].x;
    }
    if (bfr[i].y > max) {
      max = bfr[i].y;
    }
    if (bfr[i].z > max) {
      max = bfr[i].z;
    }
  }
  // Check to see if the image is solid black to avoid division by zero.
  // Since it is all black there isn't a reason to scale it.
  if (max == 0.f) {
    max = 1.f;
  }

  return max;
}