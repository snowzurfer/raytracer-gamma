
#include <vec.h>
#include <CL/cl.h>

// EPSILON is a tolerance value for floating point roundoff error.
// It is used in many calculations where we want to err
// on a certain side of a threshold, such as determining
// whether or not a point is inside a solid or not,
// or whether a point is at least a minimum distance
// away from another point.
const float kEPSILON = 1.0e-6f;

// Structs definition
struct Ray
{
  Vec origin;
  Vec dir;
};
// EO Struct

struct Material
{
  Vec matteColour;
  Vec glossColour;
  float opacity;
};

struct Sphere
{
  Vec pos;
  cl_float radius;
  struct Material material;
};

struct Light
{
  Vec pos;
  /*Vec dir*/;
  Vec col;
};

void setMatMatte(struct Material *m, const Vec *col) {
  vassign(m->matteColour, *col);
}
void setMatGloss(struct Material *m, const Vec *col) {
  vassign(m->glossColour, *col);
}
void setMatOpacity(struct Material *m, const float opacity) {
  m->opacity = opacity;
}
void setMatteGlossBalance(struct Material *m, const float glossFactor, 
  const Vec *matte, const Vec *gloss) 
{
  // Scale the glossy and matte parts of reflected light 
  // so that an object does not reflect more light than hits it.
  Vec newMatte; vsmul(newMatte, (1.0 - glossFactor), *matte);
  setMatMatte(m, matte);
  Vec newGloss; vsmul(newGloss, glossFactor, *gloss);
  setMatGloss(m, gloss);
}
// Intersection of a sphere with a ray; it returns if the collision
// was found and the parameter for the distance from the ray's
// origin to the intersection
bool raySphere(struct Sphere * sphere, struct Ray *ray, float *t) {
  // Calculate the coefficients of the quadratic equation
  //     au^2 + bu + c = 0.
  // Solving this equation gives us the value of u
  // for any intersection points.
  Vec displacement;
  vsub(displacement, ray->dir, sphere->pos);

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