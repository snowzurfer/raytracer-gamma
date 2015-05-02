
// Includes
#include <gpu_raytracer.h>


#include <cstdio>
#include <CL/cl.h>
#include <cstdlib>
#include <fstream>


namespace rtg {

  GPURaytracer::GPURaytracer(const unsigned int imgWidth,
    const unsigned int imgHeight, const float aliasFactor) :
    Raytracer(imgWidth, imgHeight, aliasFactor) 
  {

  }

  void GPURaytracer::setup(struct Sphere *hSpheres,
    const unsigned int sphNum, struct Light *hLights,
    const unsigned int lgtNum)
  {
    // Error code returned from openCL calls
    int err;


    cl_uint numPlatforms;

    // Find number of platforms
    err = clGetPlatformIDs(0, NULL, &numPlatforms);
    checkError(err, "Finding platforms");
    if (numPlatforms == 0) {
      printf("Found 0 platforms!\n");

      return;
    }

    
    // Get all platforms
    cl_platform_id *platform = (cl_platform_id *)malloc(sizeof(cl_platform_id)* numPlatforms);
    err = clGetPlatformIDs(numPlatforms, platform, NULL);
    checkError(err, "Getting platforms");



    // Define an ID for the device
    cl_device_id deviceId = 0;
    // Secure a GPU
    for (int i = 0; i < numPlatforms; i++) {
      err = clGetDeviceIDs(platform[i], CL_DEVICE_TYPE_GPU, 1, &deviceId, NULL);
      if (err == CL_SUCCESS) {
        break;
      }
    }

    // Once a device has been obtained, print out its info
    err = output_device_info(deviceId);
    checkError(err, "Printing device output");



    // Create a context for the GPU
    gpuContext_ = clCreateContext(NULL, 1, &deviceId, NULL, NULL, &err);
    checkError(err, "Creating context");




    // Create a command queue
    commandsGPU_ = clCreateCommandQueue(gpuContext_, deviceId, NULL, &err);
    checkError(err, "Creating command queue");




    // Load the kernel code
    std::ifstream sourceFstream("raytrace_kernel.cl");
    std::string source((std::istreambuf_iterator<char>(sourceFstream)),
      std::istreambuf_iterator<char>());

    // Create a program from the source
    const char* str = source.c_str();
    program_ = clCreateProgramWithSource(gpuContext_, 1, &str, NULL, &err);
    checkError(err, "Creating program");




    // Compile the program
    err = clBuildProgram(program_, 0, NULL, "-I C:\Drive\Alberto\Projects\Code\C++\raytracer_gamma\raytracer_gamma", NULL, NULL);
    // If there were compilation errors
    if (err != CL_SUCCESS) {
      // Print out compilation log
      size_t len;
      char buffer[2048];

      printf("Error: Failed to build program executable!\n%s\n", err_code(err));
      clGetProgramBuildInfo(program_, deviceId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      printf("%s\n", buffer);

      return;
    }




    // Create the kernel
    koRTG_ = clCreateKernel(program_, "raytrace", &err);
    checkError(err, "Creating kernel");




    // Create the list of spheres and lights in device memory
    cl_mem dSpheres_ = clCreateBuffer(gpuContext_, CL_MEM_READ_WRITE,
      sizeof(struct Sphere) * sphNum, NULL, &err);
    checkError(err, "Creating buffer for spheres");
    cl_mem dLights_ = clCreateBuffer(gpuContext_, CL_MEM_READ_WRITE,
      sizeof(struct Light) * lgtNum, NULL, &err);
    checkError(err, "Creating buffer for lights");
    cl_mem dPixelBuffer_ = clCreateBuffer(gpuContext_, CL_MEM_WRITE_ONLY,
      imgW_ * imgH_ * sizeof(Vec), NULL, &err);
    checkError(err, "Creating buffer for pixels");

    // Write data from host into device memory (fill the buffers with
    // the host arrays)
    err = clEnqueueWriteBuffer(commandsGPU_, dSpheres_, CL_TRUE, 0,
      sizeof(struct Sphere) * sphNum, hSpheres, 0, NULL, NULL);
    checkError(err, "Copying hSperes in dSpheres");
    err = clEnqueueWriteBuffer(commandsGPU_, dLights_, CL_TRUE, 0,
      sizeof(struct Light) * lgtNum, hLights, 0, NULL, NULL);
    checkError(err, "Copying hLights into dLights");

    cl_int   status;
    cl_uint maxDims;
    cl_uint recommWorkSize;
    cl_event events[2];
    size_t maxWorkGroupSize;

    /**
    * Query device capabilities.
    */
    status = clGetKernelWorkGroupInfo(
      koRTG_,
      deviceId,
      CL_KERNEL_WORK_GROUP_SIZE,
      sizeof(cl_uint),
      &recommWorkSize,
      NULL);
    checkError(err, "Getting recommended work group size");
    printf("Recommended workgroup size for this kernel: %d\n\n", recommWorkSize);



    // Set the work item size to be as the one recommended by
    // the implementation depending on the kernel itself
    localWorkSize_ = recommWorkSize;

    // If the global work size is not a multiple of the local
    // work size
    if (globalWorkSize_ % localWorkSize_ != 0) {
      // Pad the global work size to bea multiple of the local
      // work size
      globalWorkSize_ = (globalWorkSize_ / localWorkSize_ + 1) * localWorkSize_;
    }

    float zoomFactor = -4.f;

    // Set kernel arguments
    err = clSetKernelArg(koRTG_, 0, sizeof(cl_mem), &dSpheres_);
    err |= clSetKernelArg(koRTG_, 1, sizeof(unsigned int), &sphNum);
    err |= clSetKernelArg(koRTG_, 2, sizeof(cl_mem), &dLights_);
    err |= clSetKernelArg(koRTG_, 3, sizeof(unsigned int), &lgtNum);
    err |= clSetKernelArg(koRTG_, 4, sizeof(unsigned int), &imgH_);
    err |= clSetKernelArg(koRTG_, 5, sizeof(unsigned int), &imgW_);
    err |= clSetKernelArg(koRTG_, 6, sizeof(float), &zoomFactor);
    err |= clSetKernelArg(koRTG_, 7, sizeof(float), &aliasFactor_);
    err |= clSetKernelArg(koRTG_, 8, sizeof(cl_mem), &dPixelBuffer_);
    err |= clSetKernelArg(koRTG_, 9, sizeof(struct Sphere) * sphNum, NULL);
    err |= clSetKernelArg(koRTG_, 10, sizeof(struct Light) * lgtNum, NULL);
    checkError(err, "Setting kernel arguments");
  }

  void GPURaytracer::raytrace() {
    // Error code returned from openCL calls
    int err;

    printf("Enqueueing kernel...\t");

    err = clEnqueueNDRangeKernel(commandsGPU_, koRTG_, 1, NULL,
      &globalWorkSize_, &localWorkSize_, 0, NULL, NULL);
    checkError(err, "Enqueueing kernel");

    // Wait for the commands in the queue to be executed
    err = clFinish(commandsGPU_);
    checkError(err, "Waiting for commands to finish");
  }

  void GPURaytracer::cleanup() {
    clReleaseMemObject(dPixelBuffer_);
    clReleaseMemObject(dLights_);
    clReleaseMemObject(dSpheres_);
    clReleaseProgram(program_);
    clReleaseKernel(koRTG_);
    clReleaseCommandQueue(commandsGPU_);
    clReleaseContext(gpuContext_);
  }

  void GPURaytracer::readResult(Vec *destBuffer) {
    // Error code returned from openCL calls
    int err;

    err = clEnqueueReadBuffer(commandsGPU_, dPixelBuffer_, CL_TRUE, 0,
      imgW_ * imgH_ * sizeof(Vec), destBuffer, 0, NULL, NULL);
    // If the reading operation didn't complete successfully
    if (err != CL_SUCCESS) {
      printf("Error: Failed to read output buffer!\n%s\n", err_code(err));

      // Exit
      exit(1);
    }
  }
}
// EO Namespace