
// Includes
#include <sphere.h>
#include <ray.h>
#include <CL/cl.h>

bool raySphere(struct Sphere * sphere, struct Ray *ray, float *t) {
  // Calculate the coefficients of the quadratic equation
  //     au^2 + bu + c = 0.
  // Solving this equation gives us the value of u
  // for any intersection points.
  Vec displacement;
  vsub(displacement, ray->dir, sphere->origin);

  const cl_float a = vdot(ray->dir, ray->dir);
  const cl_float b = 2.0f * vdot(ray->dir, displacement);
  const cl_float c = (vdot(displacement, displacement)) - (sphere->radius * sphere->radius);

  // Calculate the radicand of the quadratic equation solution formula.
  // The radicand must be non-negative for there to be real solutions.
  const cl_float radicand = (b*b) - (4.0f * a * c);
  if (radicand >= 0.0) {
    // There are two intersection solutions, one involving
    // +sqrt(radicand), the other -sqrt(radicand).
    // Check both because there are weird special cases,
    // like the camera being inside the sphere,
    // or the sphere being behind the camera (invisible).
    const cl_float root = sqrt(radicand);
    const cl_float denom = 2.0f * a;
    const cl_float u[2] = {
      (-b + root) / denom,
      (-b - root) / denom
    };

    // Select the smallest param intersection value
    int smaller = 0;
    if (u[0] > kEPSILON && u[0] < u[1]) {
      *t = u[0];
      return true;
    }
    else if (u[1] > kEPSILON && u[1] < u[0]) {
      *t = u[1];
      return true;
    }
    

    return false;
  }

}