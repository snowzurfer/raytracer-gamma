#ifndef _GPU_RAYTRACER_H
#define _GPU_RAYTRACER_H

// Includes
#include <raytracer.h>
#include <CL/cl.h>
#include <device_info.h>
#include <err_code.h>

namespace rtg {

  class GPURaytracer : public Raytracer
  {
  public:
    // Ctor
    GPURaytracer(const unsigned int imgWidth,
      const unsigned int imgHeight,
      const float aliasFactor);

    // Dtor
    ~GPURaytracer() {};

    // Setup method
    void setup(struct Sphere *spheres, 
      const unsigned int sphNum, 
      struct Light *lights,
      const unsigned int lgtNum);

    void raytrace();

    void cleanup();

    void readResult(Vec *destBuffer);

  private:
    // Context for the GPU
    cl_context gpuContext_;

    // Command queue
    cl_command_queue commandsGPU_;

    // OpenCL program compiled from source
    cl_program program_;

    // Kernel created from the program
    cl_kernel koRTG_;

    // Global work size
    size_t globalWorkSize_;

    // Local work size
    size_t localWorkSize_;

    // Device buffers
    cl_mem dSpheres_, dLights_, dPixelBuffer_;
  };
}
// EO Namespace

#endif