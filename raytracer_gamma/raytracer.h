#include "vec.h"
#include "algebra.h"

// EPSILON is a tolerance value for floating point roundoff error.
// It is used in many calculations where we want to err
// on a certain side of a threshold, such as determining
// whether or not a point is inside a solid or not,
// or whether a point is at least a minimum distance
// away from another point.

// Structs definition
struct Ray
{
  Vec origin;
  Vec dir;
  Vec intensity;
};
// EO Struct

struct Material
{
  Vec matteColour;
  Vec glossColour;
  float opacity;
  float refractiveIndex;
};

struct Sphere
{
  Vec pos;
  float radius;
  struct Material material;
};

struct Light
{
  Vec pos;
  /*Vec dir*/;
  Vec col;
};

struct Intersection {
  // Object which intersects with the ray
  struct Sphere object;

  // Point of intersection
  Vec point;
  // Normal at intersection
  Vec normal;

  // Distance non-sqrted
  float squaredDist;
};

// Determine the matte reflection contribution to the illumination
// of an intersection point
Vec calculateRefraction(
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
  float &outReflectionFactor);


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
  setMatMatte(m, &newMatte);
  Vec newGloss; vsmul(newGloss, glossFactor, *gloss);
  setMatGloss(m, &newGloss);
}
void setMatRefractivityIndex(struct Material *m, const float refIndex) {
  m->refractiveIndex = refIndex;
}

// Intersection of a sphere with a ray; it returns if the collision
// was found and the parameter for the distance from the ray's
// origin to the intersection
bool raySphere(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere * sphere, struct Ray *ray, float *t) {
  const float kEPSILON = 1.0e-6f;

  // Calculate the coefficients of the quadratic equation
  //     au^2 + bu + c = 0.
  // Solving this equation gives us the value of u
  // for any intersection points.
  Vec displacement;
  vsub(displacement, ray->origin, sphere->pos);

  const float a = vdot(ray->dir, ray->dir);
  const float b = 2.0f * (vdot(ray->dir, displacement));
  const float c = (vdot(displacement, displacement)) - (sphere->radius * sphere->radius);

  // Calculate the radicand of the quadratic equation solution formula.
  // The radicand must be non-negative for there to be real solutions.
  const float radicand = (b*b) - (4.0f * a * c);
  if (radicand >= 0.0) {
    // There are two intersection solutions, one involving
    // +sqrt(radicand), the other -sqrt(radicand).
    // Check both because there are weird special cases,
    // like the camera being inside the sphere,
    // or the sphere being behind the camera (invisible).
    const float root = sqrt(radicand);
    const float denom = 2.0f * a;
    const float u[2] = {
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
  }

  return false;
}

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
struct Intersection *intersection)
{
  const float kMaxRenderDist = 1000.f;

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
        intersection->object = *(spheres + i);
      }
    }
  }

  // If at least one intersection has been found
  return found;
}

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

bool isSignificant(const Vec *colour) {
  const float kMinOpticalIntesity = 0.001f;

  return (colour->x >= kMinOpticalIntesity) ||
    (colour->y >= kMinOpticalIntesity) ||
    (colour->z >= kMinOpticalIntesity);
}

bool hasClearLineOfSight(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
  const Vec *ptA, const Vec *ptB) {
  // Calculate direction from A to B
  Vec dir; vsub(dir, *ptA, *ptB);
  const float squaredDistGap = vdot(dir, dir);

  // Construct a ray
  struct Ray ray;
  ray.dir = dir;
  ray.origin = *ptA;

  // Intersection to be returned from the function
  struct Intersection closest;
  // For each object in the scene
  if (calcIntersection(spheres, sphNum, &ray, &closest)) {

    // We found the closest intersection, but it is only
    // a blocker if it is closer to point1 than point2 is.
    // If the closest intersection is farther away than
    // point2, there is nothing on this object blocking
    // the line of sight.

    if (closest.squaredDist < squaredDistGap) {
      // We found a surface that is definitely blocking
      // the line of sight.  No need to keep looking!

      return false;
    }
  }

  // We couldn't find any objects blocking the path
  return true;
}

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
  const struct Intersection *intersection)
{
  // Colour to be computed and returned
  Vec colourSum; vinit(colourSum, 0.f, 0.f, 0.f);

  // For each light source
  for (int i = 0; i < lgtNum; ++i) {
    // Create a copy of the current light source
    const struct Light light = *(lights + i);

    // If there isn't any object in the way between the intersection point
    // and the light source
    if (hasClearLineOfSight(spheres, sphNum, &intersection->point,
      &light.pos))
    {
      // Distance vector from intersection point to light source
      Vec dist; vsub(dist, light.pos, intersection->point);
      // Direction
      Vec dir = dist; vnorm(dir);

      // Calculate the incidence of the light ray with the surface normal
      const float incidence = vdot(intersection->normal, dir);

      // If the dot prod is negative it means that is hitting the
      // surface from the inside of the object.
      // If zero it's grazing the edge of the object.
      // Only when the dot prod is positive this light source contributes
      if (incidence > 0.f) {

        // Calculate contribution based on the inverse-square law
        // and the incidence value
        const float distMagSquared = vdot(dist, dist);
        const float intensity = incidence / distMagSquared;

        // Scale light's colour depending on the intensity
        Vec scaledLightColour; vsmul(scaledLightColour, intensity, light.col);

        // Add it to the total sum
        vadd(colourSum, colourSum, scaledLightColour);
      }

    }
  }

  return colourSum;
}

// Calculate the reflection using fresnel's equations
float polarisedReflection(
  float n1,		// The source material's index of refraction
  float n2,		// The target material's index of refraction
  float cosA1,	// Incident or outgoing ray's angle cosine
  float cosA2)	// Incident or outgoing ray's angle cosine
{
  const float kEPSILON = 1.0e-6f;

  const float left = n1 * cosA1;
  const float right = n2 * cosA2;
  double numerator = left - right;
  double denominator = left + right;

  // Square the denominator
  denominator *= denominator;

  // If the denominator is approx zero
  if (denominator < kEPSILON) {
    // Assume complete reflection
    return 1.f;
  }

  // Compute reflection
  float reflection = (numerator * numerator) / denominator;

  // If the reflection value is more than one
  if (reflection > 1.f) {
    // Cap it
    reflection = 1.f;
  }

  return reflection;

}






Vec rayTrace(
#ifdef GPU_KERNEL
OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Light *lights, const unsigned int lgtNum,
struct Ray ray, struct Material *refractiveMaterial,
  int traceDepth)
{
  const int kMaxTraceDepth = 2;

  // Colour to be computed and returned
  Vec colourSum; vinit(colourSum, 0.f, 0.f, 0.f);

  // Intersection to be returned from the function
  struct Intersection intersection;
  if (calcIntersection(spheres, sphNum, &ray, &intersection)) {
    // Check for end of recursion condition
    if (traceDepth <= kMaxTraceDepth) {
      // If the ray still has significant intensity
      if (isSignificant(&ray.intensity)) {
        // Calculate the opacity and transparency available
        // of the light ray.
        const float opacity = intersection.object.material.opacity;
        const float transparency = 1.f - opacity;

        // If the object is opaque
        if (opacity > 0.f) {
          // Calculate matte colour
          Vec calcTemp; vmul(calcTemp, ray.intensity, intersection.object.material.matteColour);
          vsmul(calcTemp, opacity, calcTemp);
          Vec matteCalcResult = calculateMatte(spheres, sphNum, lights,
            lgtNum, &intersection);

          vmul(calcTemp, matteCalcResult, calcTemp);

          vadd(colourSum, calcTemp, colourSum);

        }

        // Variable to hold the fraction of light which is reflected.
        // It's calculated by the calculateRefraction function if there
        // is transparency.
        float refractiveReflectionFactor = 0.f;

        // If there is transparency
        if (transparency > 0.f) {
          // Calculate the refraction
          Vec refrCalcResult = calculateRefraction(
            spheres,
            sphNum,
            lights,
            lgtNum,
            &intersection,
            ray,
            traceDepth + 1,
            refractiveMaterial,
            refractiveReflectionFactor);

          vadd(colourSum, refrCalcResult, colourSum);
        }
      }
    }
  }
  else {
    // Return background colour
    vmul(colourSum, ray.intensity, refractiveMaterial->matteColour)
  }

  return colourSum;
}



// Determine the matte reflection contribution to the illumination
// of an intersection point
Vec calculateRefraction(
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
  float &outReflectionFactor)
{
  // Store the cos of the angle between the incident ray and the surface normal
  float cosA1 = vdot(incidentRay.dir, intersection->normal);
  float sinA1 = 0.f;

  // if the cos is less than one
  if (cosA1 <= -1.0) {
    // The incident ray points in the opposite direction as the normal
    // vector and therefore the ray is entering the solid exactly
    // perpendicularly to the surface.
    // Hence clamp to lower limit
    cosA1 = -1.f;
    sinA1 = 0.f;
  }
  else if (cosA1 >= +1.f) {
    // The incident ray points in the same direction as the normal
    // vector and therefore the ray is entering the solid exactly
    // perpendicularly to the surface.
    // Hence clamp to upper limit
    cosA1 = 1.f;
    sinA1 = 0.f;
  }
  else {
    // The incident ray is entering the solid at some positive angle;
    // Hence calculate the sine of such angle using the trig
    // identity cos^2 + sin^2 = 1
    sinA1 = sqrt(1.0 - (cosA1 * cosA1));
  }

  // Calculate the ratio between the source and the target material 
  // refractive indices. This is necessary for snell's law
  const float refIndRatio =
    refractiveMaterial->refractiveIndex /
    intersection->object.material.refractiveIndex;

  // Compute the sine of the refracted angle using Snell's law:
  // the sine of the refracted angle with respect to the normal
  const float sinA2 = refIndRatio * sinA1;

  // If the refracted angle is more than 90 degrees
  if (sinA2 <= -1.f || sinA2 >= 1.f) {
    // The ray experiences total internal reflection, hence
    // the refracted ray doesn't exist.

    // Inform the caller that the ray experiences total reflection
    outReflectionFactor = 1.f;
    // Return a zero contribution
    Vec contrib; vinit(contrib, 0.f, 0.f, 0.f);
    return contrib;
  }

  // Since there is some refracted light, determine its direction.
  // A quadratic equation needs to be solved for this reason
  float roots[2];
  const int numSolutions = solveQuadratic(
    1.f,
    (2.f * cosA1),
    (1.f - (1.f / (refIndRatio * refIndRatio))),
    roots);

  // Since there usually are 2 solutions for the roots, the correct
  // one needs to be determined.
  // This solution is the one which causes the light ray to bend
  // the smallest angle when comparing its direction with the one
  // of the incident ray.
  // It can be found by finding the refracted ray with the largest
  // positive dot product.

  // Setup the max alignment value for the dot product 
  float maxAlignment = -0.1;

  // Setup a vector to hold the direction of the refracted ray
  Vec refractionDir = { 0.f, 0.f, 0.f };

  // For each solution
  for (int i = 0; i < numSolutions; ++i) {
    // Calculate the refraction for this solution
    Vec normalTimesSolution;
    vsmul(normalTimesSolution, roots[i], intersection->normal);
    Vec currentSolRef; vadd(currentSolRef, incidentRay.dir, normalTimesSolution);

    float alignment = vdot(incidentRay.dir, currentSolRef);

    // If the alignment for this solution is larger than the
    // largest calculated before
    if (alignment > maxAlignment) {
      // Set this to be the largest
      maxAlignment = alignment;
      vassign(refractionDir, currentSolRef);
    }
  }

  // TODO Do check for maxAlignment <= 0.f here

  // Determine the cosine of the angle of the refracted ray 
  float cosA2 = sqrt(1.f - (sinA2 * sinA2));
  // If the cosine is less than zero
  if (cosA2 < 0.f) {
    // Set its polarity to match the one of cosA1
    cosA2 = -cosA2;
  }

  // Compute the fraction of the light which is reflected
  // using the Fresnel's equations, assuming uniform
  // polarisation of light.
  const float Rs = polarisedReflection(
    refractiveMaterial->refractiveIndex,
    intersection->object.material.refractiveIndex,
    cosA1,
    cosA2);
  const float Rp = polarisedReflection(
    refractiveMaterial->refractiveIndex,
    intersection->object.material.refractiveIndex,
    cosA2,
    cosA1);

  // Compute total reflection
  outReflectionFactor = ((Rs + Rp) *0.5);

  // The fraction of light which is not reflected contributes
  // to refraction and therefore the incoming light ray's intensity
  // is diminished by this factor


  // Construct the refracted ray
  struct Ray refractedRay;
  vsmul(refractedRay.intensity, (1.f - outReflectionFactor), incidentRay.intensity);
  vassign(refractedRay.origin, intersection->point);
  vassign(refractedRay.dir, refractionDir);

  return rayTrace(
    spheres,
    sphNum,
    lights,
    lgtNum,
    refractedRay,
    &intersection->object.material,
    traceDepth);
}