
// Includes
#include <cpu_raytracer.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>

// Functions which raytrace on CPU; mostly similar to the
// ones used for the kernel
// Determine the matte reflection contribution to the illumination
// of an intersection point
struct Ray calculateReflection(
struct Intersection *intersection,
struct Ray *incidentRay);

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
  const float kEPSILON = 1.0e-5f;

  // The result of the test
  bool result = false;

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

    // Smallest value found
    float smallestT = 10000.f;
    // For each root
    for (int i = 0; i < 2; ++i) {
      // If the root is more than a threshold value
      if (u[i] > kEPSILON) {
        // If the root is larger than the previous
        // smallest
        if (u[i] < smallestT) {
          // Set it to to be the smallest
          smallestT = u[i];
          // Report success of collision
          result = true;

        }
      }
    }

    // Set the root to be returned
    *t = smallestT;
  }

  return result;
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

        // Set the newly smallest
        minT = t;
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


// Returns the index to the sphere which contains the point, if any
int primaryContainer(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
  const Vec *pt) {

  const float kEPSILON = 1.0e-6f;

  // For each sphere in the scene
  for (int i = 0; i < sphNum; ++i) {
    // Add a little bit to the actual radius to be more tolerant
    // of rounding errors that would incorrectly exclude a 
    // point that should be inside the sphere.
    const float r = spheres[i].radius + kEPSILON;

    // A point is inside the sphere if the square of its distance 
    // from the center is within the square of the radius.
    Vec dist; vsub(dist, *pt, spheres[i].pos);
    if ((vdot(dist, dist) <= (r * r))) {
      return i;
    }
  }

  return -1;
}

bool hasClearLineOfSight(
#ifdef GPU_KERNEL
  OCL_GLOBAL_BUFFER
#endif
struct Sphere *spheres, const unsigned int sphNum,
  const Vec *ptA, const Vec *ptB) {
  // Calculate direction from A to B
  Vec dir; vsub(dir, *ptB, *ptA);
  const float squaredDistGap = vdot(dir, dir);

  // Construct a ray
  struct Ray ray;
  ray.dir = dir;
  vnorm(ray.dir);
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
struct Ray ray, struct Material refractiveMaterial,
  int traceDepth)
{
  const int kMaxTraceDepth = RTSTACK_MAXSIZE - 1;

  // Colour to be computed and returned
  Vec colourSum; vinit(colourSum, 0.f, 0.f, 0.f);

  // Create a stack to hold the
  RtStack snapshotsStack;
  rtStackInit(&snapshotsStack);

  // Create a new snapshot
  RtSnapshot currSnapshot;
  currSnapshot.ray = ray;
  currSnapshot.traceDepth = traceDepth;
  currSnapshot.stage = 0;
  currSnapshot.colour = colourSum;
  currSnapshot.refractiveMat = refractiveMaterial;

  // Push it into the stack
  rtStackPush(&snapshotsStack, &currSnapshot);

  bool loop = !(rtStackIsEmpty(&snapshotsStack));
  // While the stack is not empty
  while (loop) {
    // Read top / pop an element from the stack
    currSnapshot = *rtStackTop(&snapshotsStack);
    rtStackPop(&snapshotsStack);

    // Depending on the stage of the recursion
    // (There is code to be processed both before a recursion call
    // happens and after it has returned)
    switch (currSnapshot.stage) {
      // Before recursion for refraction
      case 0: {
        if (calcIntersection(spheres, sphNum, &currSnapshot.ray,
          &currSnapshot.intersection)) {
          // Check for end of recursion condition
          if (currSnapshot.traceDepth <= kMaxTraceDepth) {
            // If the ray still has significant intensity
            if (isSignificant(&currSnapshot.ray.intensity)) {
              // Calculate the opacity and transparency available
              // of the light ray.
              const float opacity =
                currSnapshot.intersection.object.material.opacity;
              const float transparency = 1.f - opacity;

              // If the object is opaque
              if (opacity > 0.f) {
                // Calculate matte colour
                Vec calcTemp; vmul(calcTemp, currSnapshot.ray.intensity,
                  currSnapshot.intersection.object.material.matteColour);
                vsmul(calcTemp, opacity, calcTemp);
                Vec matteCalcResult = calculateMatte(spheres, sphNum, lights,
                  lgtNum, &currSnapshot.intersection);

                if (isSignificant(&matteCalcResult)) {
                  int lol = 0;
                }

                vmul(calcTemp, matteCalcResult, calcTemp);

                vadd(currSnapshot.colour, calcTemp, currSnapshot.colour);

              }

              // Variable to hold the fraction of light which is reflected.
              // It's calculated by the calculateRefraction function if there
              // is transparency.
              float refractiveReflectionFactor = 0.f;



              // If there is transparency
              if (transparency > 0.f) {
                // Calculate a ray to pass in the calculate refraction method
                struct Ray refractionRay;
                vassign(refractionRay.dir, currSnapshot.ray.dir);
                vsmul(refractionRay.intensity, transparency, currSnapshot.ray.intensity);
                vassign(refractionRay.origin, currSnapshot.ray.origin);

                struct Material targetMaterial;

                // Calculate the refraction
                struct Ray refractedRay = calculateRefraction(
                  spheres,
                  sphNum,
                  lights,
                  lgtNum,
                  &currSnapshot.intersection,
                  refractionRay,
                  currSnapshot.traceDepth,
                  &currSnapshot.refractiveMat,
                  &targetMaterial,
                  &refractiveReflectionFactor);

                // Store the state of the current snapshot
                currSnapshot.refractiveReflectionFactor =
                  refractiveReflectionFactor;
                currSnapshot.stage = 1;

                // Push the current state
                rtStackPush(&snapshotsStack, &currSnapshot);

                // Create a new snapshot for simulating recursion
                RtSnapshot newSnapshot;
                newSnapshot.ray = refractedRay;
                newSnapshot.traceDepth = traceDepth + 1;
                newSnapshot.stage = 0;
                vinit(newSnapshot.colour, 0.f, 0.f, 0.f);
                newSnapshot.refractiveMat = targetMaterial;

                // Push the newly created snapshot
                rtStackPush(&snapshotsStack, &newSnapshot);

                // Execute a new loop
                break;
              }

              vassign(colourSum, currSnapshot.colour);
            }
          }
        }
        else {
          // Return background colour
          vmul(currSnapshot.colour, currSnapshot.ray.intensity,
            currSnapshot.refractiveMat.matteColour);
        }


        break;
      }
        // After refraction recursion
      case 1: {
        vadd(currSnapshot.colour, colourSum, currSnapshot.colour);

        // Two sources of shiny reflection must be considered:
        // a. Reflection caused by refraction
        // b. Glossy part

        // a.
        // The refractive part causes reflection of all colours equally.
        // The components of the colour are diminished based on the 
        // transparency available as calculated by calculateRefraction.
        Vec reflectionCol = { 1.f, 1.f, 1.f };
        float transparency = 1.f - currSnapshot.intersection.object.material.opacity;
        float prod = transparency * currSnapshot.refractiveReflectionFactor;
        vsmul(reflectionCol, prod, reflectionCol);


        // Add the glossy part of the reflection. This contribution
        // is diminished by the part of the light which wasn't available
        // for refraction (and therefore, reflection)
        Vec glossColContrib;
        vsmul(glossColContrib, currSnapshot.refractiveMat.opacity,
          currSnapshot.intersection.object.material.glossColour);
        vadd(reflectionCol, reflectionCol, glossColContrib);

        // Multiply by the intensity of the ray
        vmul(reflectionCol, currSnapshot.ray.intensity, reflectionCol);



        // If the contribution is significant
        //if (isSignificant(&reflectionCol)) {
        //  // Compute a ray to pass in the function
        //  //struct Ray reflectionRay;
        //  //vassign(reflectionRay.dir, currSnapshot.ray.dir);
        //  //vassign(reflectionRay.intensity, reflectionCol);
        //  //vassign(reflectionRay.origin, currSnapshot.ray.origin);


        //  //// Calculate the reflected ray
        //  //struct Ray reflectedRay = calculateReflection(
        //  //  &currSnapshot.intersection,
        //  //  &reflectionRay);

        //  //// Store the state of the current snapshot
        //  //currSnapshot.stage = 2;

        //  //// Push the current state
        //  //rtStackPush(&snapshotsStack, &currSnapshot);

        //  //// Create a new snaphot for simulating recursion
        //  ////RtSnapshot newSnapshot;
        //  //currSnapshot.ray = reflectedRay;
        //  //currSnapshot.traceDepth = traceDepth + 1;
        //  //currSnapshot.stage = 0;
        //  //vinit(currSnapshot.colour, 0.f, 0.f, 0.f);
        //  //// currSnapshot.refractiveMat;

        //  //// Push the newly created snapshot
        //  ////rtStackPush(&snapshotsStack, &currSnapshot);

        //  //// Execute a new loop
        //  continue;
        //  break;
        //}

        vassign(colourSum, currSnapshot.colour);
        break;
      }
      case 2: {
        // Add the result of the reflection to the total
        vadd(currSnapshot.colour, colourSum, currSnapshot.colour);

        // One iteration is finished, therefore save the result into the
        // main variable
        vassign(colourSum, currSnapshot.colour);

        printf("case2\n");

        break;
      }

    }

    // Assign the computation of the current colour to the sum
    vassign(colourSum, currSnapshot.colour);

    loop = !(rtStackIsEmpty(&snapshotsStack));
  }

  return colourSum;
}



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
  float *outReflectionFactor)
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


  // Find the refractive index of the material after the collision point
  const float kSmallShift = 0.01f;
  {
    Vec testPt; vsmul(testPt, kSmallShift, incidentRay.dir);
    vadd(testPt, testPt, intersection->point);
    int containerNum = primaryContainer(spheres, sphNum, &testPt);

    struct Material bgMaterial;
    Vec black; vinit(black, 0.f, 0.f, 0.f);
    setMatteGlossBalance(&bgMaterial, 0.f, &black, &black);
    setMatRefractivityIndex(&bgMaterial, 1.00f);

    // If the point was within a sphere
    if (containerNum != -1) {
      // Read the material of that sphere
      *targetMaterial = spheres[containerNum].material;
    }
    else {
      // Use the ambient's material
      *targetMaterial = bgMaterial;
    }
  }

  // Calculate the ratio between the source and the target material 
  // refractive indices. This is necessary for snell's law
  const float refIndRatio =
    refractiveMaterial->refractiveIndex /
    targetMaterial->refractiveIndex;

  // Compute the sine of the refracted angle using Snell's law:
  // the sine of the refracted angle with respect to the normal
  const float sinA2 = refIndRatio * sinA1;

  // If the refracted angle is more than 90 degrees
  if (sinA2 <= -1.f || sinA2 >= 1.f) {
    // The ray experiences total internal reflection, hence
    // the refracted ray doesn't exist.

    // Inform the caller that the ray experiences total reflection
    *outReflectionFactor = 1.f;
    // Return a zero contribution
    Vec contrib; vinit(contrib, 0.f, 0.f, 0.f);

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
  if (cosA1 < 0.f) {
    // Set its polarity to match the one of cosA1
    cosA2 = -cosA2;
  }

  // Compute the fraction of the light which is reflected
  // using the Fresnel's equations, assuming uniform
  // polarisation of light.
  const float Rs = polarisedReflection(
    refractiveMaterial->refractiveIndex,
    targetMaterial->refractiveIndex,
    cosA1,
    cosA2);
  const float Rp = polarisedReflection(
    refractiveMaterial->refractiveIndex,
    targetMaterial->refractiveIndex,
    cosA2,
    cosA1);

  // Compute total reflection
  *outReflectionFactor = ((Rs + Rp) *0.5);

  // The fraction of light which is not reflected contributes
  // to refraction and therefore the incoming light ray's intensity
  // is diminished by this factor


  // Construct the refracted ray
  struct Ray refractedRay;
  vsmul(refractedRay.intensity, (1.f - *outReflectionFactor), incidentRay.intensity);
  vassign(refractedRay.origin, intersection->point);
  vassign(refractedRay.dir, refractionDir);
  // Shift the ray by a little amount from the surface it collided with
  /*Vec smallShift; vsmul(smallShift, kSmallShift, refractedRay.dir);
  vadd(refractedRay.origin, refractedRay.origin, smallShift);*/

  return refractedRay;
}

struct Ray calculateReflection(
struct Intersection *intersection,
struct Ray *incidentRay)
{

  // Calculate the direction of the reflected ray
  const float perp = 2.f * (vdot(incidentRay->dir, intersection->normal));
  Vec perpTimesNormal;  vsmul(perpTimesNormal, perp, intersection->normal);
  Vec reflectedDir;
  vsub(reflectedDir, incidentRay->dir, perpTimesNormal);
  // Normalise it
  vnorm(reflectedDir);

  // Compute the new reflected ray
  const float kSmallShift = 0.01f;
  struct Ray reflectedRay;
  vassign(reflectedRay.dir, reflectedDir);
  vassign(reflectedRay.origin, intersection->point);
  vassign(reflectedRay.intensity, incidentRay->intensity);
  // Shift the ray by a little amount from the surface it collided with
  Vec smallShift; vsmul(smallShift, kSmallShift, reflectedRay.dir);
  vadd(reflectedRay.origin, reflectedRay.origin, smallShift);

  // Trace the ray in the new direction
  return reflectedRay;
}

namespace rtg {

  CPURaytracer::CPURaytracer(const unsigned int imgWidth,
    const unsigned int imgHeight,
    const float aliasFactor) :
    Raytracer(imgWidth, imgHeight, aliasFactor),
    imagePtr_(NULL), spheres_(NULL),
    lights_(NULL)
  {}


  CPURaytracer::~CPURaytracer() {}

  void CPURaytracer::setup(struct Sphere *spheres,
    const unsigned int sphNum,
  struct Light *lights,
    const unsigned int lgtNum)
  {
    // Allocate memory for the pixels buffer
    unsigned int globalWorkSize = imgH_ * imgW_;
    imagePtr_ = (float *)malloc(globalWorkSize * sizeof(Vec));

    // Save pointers for later use
    spheres_ = spheres;
    lights_ = lights;

    sphNum_ = sphNum;
    lgtNum_ = lgtNum;
  }

  void CPURaytracer::raytrace() {
    unsigned int globalWorkSize = imgH_ * imgW_;

    // Screen in world coordinates
    const float kImageWorldWidth = imgW_ * 0.02; // Pic width / 50
    const float kImageWorldHeight = imgH_ * 0.02; // Pic height / 50

    // Amount to increase each step for the ray direction
    const float kRayXStep = kImageWorldWidth / ((float)imgW_);
    const float kRayYStep = kImageWorldHeight / ((float)imgH_);
    const float aspectRatio = kImageWorldWidth / kImageWorldHeight;

    // Variables holding the current step in world coordinates
    float rayX = 0.f, rayY = 0.f;

    int pixelsCounter = 0;

    // Calculate size of an alias step in world coordinates
    const float kAliasFactorStepInv = kRayXStep / aliasFactor_;
    // Calculate total size of samples to be taken
    const float kSamplesTot = aliasFactor_ * aliasFactor_;
    // Also its inverse
    const float kSamplesTotinv = 1.f / kSamplesTot;

    for (unsigned int y = 0; y < globalWorkSize; ++y, pixelsCounter += 3) {
      // Retrieve the global ID of the kernel
      const unsigned gid = y;

      
      // Calculate world position of pixel being currently worked on
      const float kPxWorldX = ((((float)(gid % imgW_) - 
        (imgW_ * 0.5f))) * kRayXStep);
      const float kPxWorldY = ((imgH_ *0.5f) - ((float)(gid / imgW_))) * kRayYStep;

      // The ray to be shot. The vantage point (camera) is at the origin,
      // and its intensity is maximum
      struct Ray ray; vinit(ray.origin, 0.f, 0.f, 0.f); vinit(ray.intensity, 1.f, 1.f, 1.f);

      // The colour of the pixel to be computed
      Vec pixelCol = { 0.f, 0.f, 0.f };

      // Mock background material
      struct Material bgMaterial;
      Vec black; vinit(black, 0.f, 0.f, 0.f);
      setMatteGlossBalance(&bgMaterial, 0.f, &black, &black);
      setMatRefractivityIndex(&bgMaterial, 1.00f);

      // For each sample to be taken
      for (int i = 0; i < aliasFactor_; ++i) {
        for (int j = 0; j < aliasFactor_; ++j) {
          // Calculate the direction of the ray
          float x = (kPxWorldX + (float)(((float)j) * kAliasFactorStepInv)) * aspectRatio;
          float y = (kPxWorldY + (float)(((float)i) * kAliasFactorStepInv));

          // Set the ray's dir and normalise it
          float zoomFactor = -4.f;
          vinit(ray.dir, x, y, zoomFactor); vnorm(ray.dir);

          // Raytrace for the current sample
          Vec currentSampleCol = rayTrace(
            spheres_, 
            sphNum_, 
            lights_, 
            lgtNum_,
            ray, 
            bgMaterial, 
            0);

          if (currentSampleCol.x > 0.f) {
            int lol = 0;
          }

          vsmul(currentSampleCol, kSamplesTotinv, currentSampleCol);

          // Compute the average
          vadd(pixelCol, pixelCol, currentSampleCol);
        }
      }

      // Write result in destination buffer
      *(imagePtr_ + pixelsCounter) = pixelCol.x;
      *(imagePtr_ + pixelsCounter + 1) = pixelCol.y;
      *(imagePtr_ + pixelsCounter + 2) = pixelCol.z;
    }

  }

  void CPURaytracer::cleanup() {
    free(imagePtr_);
  }

  void CPURaytracer::readResult(Vec *destBuffer) {
    memcpy(destBuffer, (Vec *)imagePtr_, sizeof(Vec) * imgH_ * imgW_);
  }

}
// EO Namespace