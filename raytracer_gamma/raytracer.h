
#include <vec.h>

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
  struct Sphere *object;

  // Point of intersection
  Vec point;
  // Normal at intersection
  Vec normal;

  // Distance non-sqrted
  float squaredDist;
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
  setMatMatte(m, &newMatte);
  Vec newGloss; vsmul(newGloss, glossFactor, *gloss);
  setMatGloss(m, &newGloss);
}
// Intersection of a sphere with a ray; it returns if the collision
// was found and the parameter for the distance from the ray's
// origin to the intersection
bool raySphere(struct Sphere * sphere, struct Ray *ray, float *t) {
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
bool calcIntersection(struct Sphere *spheres, const unsigned int sphNum,
struct Ray *ray, struct Intersection *intersection)
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
        intersection->object = (spheres + i);
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

bool hasClearLineOfSight(struct Sphere *spheres, const unsigned int sphNum,
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
Vec calculateMatte(struct Sphere *spheres, const unsigned int sphNum,
  const struct Light *lights, const unsigned int lgtNum,
  const struct Intersection *intersection)
{
  // Colour to be computed and returned
  Vec colourSum; vinit(colourSum, 0.f, 0.f, 0.f);

  // For each light source
  for (int i = 0; i < lgtNum; ++i) {
    // Create a reference to the current light source
    const struct Light *light = (lights + i);

    // If there isn't any object in the way between the intersection point
    // and the light source
    if (hasClearLineOfSight(spheres, sphNum, &intersection->point,
      &light->pos))
    {
      // Distance vector from intersection point to light source
      Vec dist; vsub(dist, light->pos, intersection->point);
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
        Vec scaledLightColour; vsmul(scaledLightColour, intensity, light->col);

        // Add it to the total sum
        vadd(colourSum, colourSum, scaledLightColour);
      }

    }
  }

  return colourSum;
}

Vec rayTrace(struct Sphere *spheres, const unsigned int sphNum,
            struct Light *lights, const unsigned int lgtNum,
            struct Ray *ray, struct Material *refractiveMaterial,
            int traceDepth)
{
  const int kMaxTraceDepth = 2;

  // Colour to be computed and returned
  Vec colourSum; vinit(colourSum, 0.f, 0.f, 0.f);

  // Intersection to be returned from the function
  struct Intersection intersection;
  if (calcIntersection(spheres, sphNum, ray, &intersection)) {
    // Check for end of recursion condition
    if (traceDepth <= kMaxTraceDepth) {
      // If the ray still has significant intensity
      if (isSignificant(&ray->intensity)) {
        // Calculate the opacity and transparency available
        // of the light ray.
        const float opacity = intersection.object->material.opacity;
        const float transparency = 1.f - opacity;

        // If the object is opaque
        if (opacity > 0.f) {
          // Calculate matte colour
          Vec calcTemp; vmul(calcTemp, ray->intensity, intersection.object->material.matteColour);
          vsmul(calcTemp, opacity, calcTemp);
          Vec matteCalcResult = calculateMatte(spheres, sphNum, lights,
            lgtNum, &intersection);

          vmul(calcTemp, matteCalcResult, calcTemp);

          vadd(colourSum, calcTemp, colourSum);

        }
      }
    }
  }
  else {
    // Return background colour
    vmul(colourSum, ray->intensity, refractiveMaterial->matteColour)
  }

  return colourSum;
}