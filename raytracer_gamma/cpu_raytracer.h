#ifndef _CPU_RAYTRACER_H
#define _CPU_RAYTRACER_H



#include "vec.h"
#include "algebra.h"
#include "ray.h"
#include "material.h"
#include "raytrace_stack.h"
#include "sphere.h"
#include "intersection.h"
#include "light.h"


// EPSILON is a tolerance value for floating point roundoff error.
// It is used in many calculations where we want to err
// on a certain side of a threshold, such as determining
// whether or not a point is inside a solid or not,
// or whether a point is at least a minimum distance
// away from another point.







// Functions which raytrace on CPU; mostly similar to the
// ones used for the kernel
// Determine the matte reflection contribution to the illumination
// of an intersection point
struct Ray calculateReflection(
struct Intersection *intersection,
struct Ray incidentRay);

// Determine the matte refraction contribution to the illumination
// of an intersection point
struct Ray calculateRefraction(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Light *lights, const unsigned int lgtNum,
struct Intersection *intersection,
struct Ray incidentRay,
  int traceDepth,
  const struct Material *refractiveMaterial,
struct Material *targetMaterial,
  float *outReflectionFactor);


void setMatMatte(struct Material *m, const Vec *col);
void setMatGloss(struct Material *m, const Vec *col);
void setMatOpacity(struct Material *m, const float opacity);
void setMatteGlossBalance(struct Material *m, const float glossFactor,
  const Vec *matte, const Vec *gloss);
void setMatRefractivityIndex(struct Material *m, const float refIndex);



// Intersection of a sphere with a ray; it returns if the collision
// was found and the parameter for the distance from the ray's
// origin to the intersection
bool raySphere(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere * sphere, struct Ray *ray, float *t);

// Calculates eventual intersection of ray with the scene and returns
// the intersection a populated intersection object
bool calcIntersection(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
struct Ray *ray,
#ifdef GPU_KERNEL
  OCL_PRIVATE_BUFFER
#endif
struct Intersection *intersection);

// Checks to see if the ray intersects with any of the objects and if so
// it returns immediately
/*bool intersectsExit(struct Sphere *spheres, const unsigned int sphNum,
struct Ray *ray, struct Intersection *intersection)
{
float minT = kMaxRenderDist;
bool found = false;

for (int i = 0; i < sphNum; ++i) {
// Find closest intersection
float t;

// If there is an intersection
if (raySphere(&spheres[i], ray, &t)) {
// If the newly found parameter is smaller than the
// last smallest
if (t < minT) {

// Compute intersection position
vsmul(intersection->point, t, ray->dir);
vadd(intersection->point, ray->origin, intersection->point);

// Compute normal
vsub(intersection->normal, intersection->point, spheres[i].pos);
vnorm(intersection->normal);

// Distance squared
Vec dist; vsmul(dist, t, ray->dir);
intersection->squaredDist = vdot(dist, dist);

found = true;
}
}
}

// If at least one intersection has been found
return found;
}*/

bool isSignificant(const Vec *colour);


// Returns the index to the sphere which contains the point, if any
int primaryContainer(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
  const Vec *pt);

bool hasClearLineOfSight(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
  const Vec *ptA, const Vec *ptB);

// Determine the matte reflection contribution to the illumination
// of an intersection point
Vec calculateMatte(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
  const struct Light *lights, const unsigned int lgtNum,
  const struct Intersection *intersection);

// Calculate the reflection using fresnel's equations
float polarisedReflection(
  float n1,		// The source material's index of refraction
  float n2,		// The target material's index of refraction
  float cosA1,	// Incident or outgoing ray's angle cosine
  float cosA2);	// Incident or outgoing ray's angle cosine







Vec rayTrace(
#ifdef GPU_KERNEL
OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Light *lights, const unsigned int lgtNum,
struct Ray ray, struct Material refractiveMaterial,
  int traceDepth);



// Determine the matte reflection contribution to the illumination
// of an intersection point
struct Ray calculateRefraction(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Light *lights, const unsigned int lgtNum,
struct Intersection *intersection,
struct Ray incidentRay,
  int traceDepth,
  const struct Material *refractiveMaterial,
struct Material *targetMaterial,
  float *outReflectionFactor);

struct Ray calculateReflection(
struct Intersection *intersection,
struct Ray incidentRay);


// Class

#include <raytracer.h>

namespace rtg {

  class CPURaytracer : public Raytracer
  {
  public:
    // Ctor
    CPURaytracer(const unsigned int imgWidth,
      const unsigned int imgHeight,
      const float aliasFactor);

    // Dtor
    ~CPURaytracer();


    // Setup method
    void setup(struct Sphere *spheres,
      const unsigned int sphNum,
    struct Light *lights,
      const unsigned int lgtNum);

    void raytrace();

    void cleanup();

    void readResult(Vec *destBuffer);

  private:
    // Pointer to a buffer to hold the results of the raytracing
    float *imagePtr_;

    // Pointers to the scene
    struct Sphere *spheres_;
    struct Light *lights_;

    // Number of elements
    unsigned int sphNum_, lgtNum_;
  };

}
// EO Namespace


#endif