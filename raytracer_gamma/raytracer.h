#ifndef _RAYTRACER_H
#define _RAYTRACER_H

// Includes
#include <light.h>
#include <sphere.h>

namespace rtg {

  class Raytracer
  {
  public:
    // Ctor
    Raytracer(const unsigned int imgWidth,
      const unsigned int imgHeight,
      const float aliasFactor) :
      imgW_(imgWidth), imgH_(imgHeight),
      aliasFactor_(aliasFactor)
    {};

    // Dtor
    virtual ~Raytracer() {};

    // Setup memory before raytracing
    virtual void setup(struct Sphere *spheres,
      const unsigned int sphNum,
      struct Light *lights,
      const unsigned int lgtNum) = 0;

    // Raytrace
    virtual void raytrace() = 0;

    // Cleanup buffers once they've been used
    virtual void cleanup() = 0;

    // Read the result and give it to the host
    virtual void readResult(Vec *destBuffer) = 0;

  protected:
    // Image sizes
    const unsigned int imgW_, imgH_;

    // Alias factor
    const float aliasFactor_;
  };


}
// EO Namespace


#endif