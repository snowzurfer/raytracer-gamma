// Hello world with OpenCL


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <device_info.h>
#include <err_code.h>
#include <fstream>
#include <algorithm>
#include <raytracer.h>
#include <chrono>
#include <ppm.h>
#include <string>
#include <iostream>


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif

#define LENGTH (1024)    // length of vectors a, b, and c


// Colours
const Vec whiteCol = { 8.f, 8.f, 8.f };
const Vec lowerWhite = { 0.5f, 0.5f, 0.5f };
const Vec redCol = { 0.8f, 0.1f, 0.1f };
const Vec greenCol = { 0.1f, 0.8f, 0.1f };
const Vec col1 = { 0.01f, 0.8f, 0.01f };

// Setup the default scene
void setupBaseScene(struct Sphere *spheres, struct Light *lights) {
  // Setup materials
  struct Material ballMaterial1; // White
  Vec bm1Gloss; vassign(bm1Gloss, greenCol);
  Vec bm1Matte; vassign(bm1Matte, greenCol);
  setMatOpacity(&ballMaterial1, 0.8f);
  setMatteGlossBalance(&ballMaterial1, 0.2f, &bm1Matte, &bm1Gloss);
  setMatRefractivityIndex(&ballMaterial1, 1.5500f);

  struct Material ballMaterial2; // Red
  Vec bm2Gloss; vassign(bm2Gloss, redCol);
  Vec bm2Matte; vassign(bm2Matte, redCol);
  setMatOpacity(&ballMaterial2, 0.3f);
  setMatteGlossBalance(&ballMaterial2, 0.95f, &bm2Matte, &bm2Gloss);
  setMatRefractivityIndex(&ballMaterial2, 1.5500f);

  struct Material ballMaterial3; // Red
  Vec bm3Gloss; vassign(bm3Gloss, col1);
  Vec bm3Matte; vassign(bm3Matte, col1);
  setMatOpacity(&ballMaterial3, 0.6f);
  setMatteGlossBalance(&ballMaterial3, 0.0, &bm3Matte, &bm3Gloss);
  setMatRefractivityIndex(&ballMaterial3, 1.5500f);

  // Setup spheres
  spheres[0].material = ballMaterial1;
  vinit(spheres[0].pos, -9.f, 0.f, -13.f);
  spheres[0].radius = 5.f;
  spheres[1].material = ballMaterial2;
  vinit(spheres[1].pos, -4.f, 1.5f, -5.f);
  spheres[1].radius = 2.f;
  spheres[2].material = ballMaterial3;
  vinit(spheres[2].pos, 1.f, -1.f, -7.f);
  spheres[2].radius = 3.f;

  // Setup light sources
  vinit(lights[0].pos, -45.f, 10.f, 85.f);
  vassign(lights[0].col, lowerWhite);
  vinit(lights[1].pos, 20.f, 60.f, -5.f);
  vassign(lights[1].col, lowerWhite);
}


int main(int argc, char** argv)
{

  // Read the parameters
  if (argc == 1) {
    // Use default settings
  }
  else if (argc == 3) {
    // Read the passed parameters

    // Read the first parameter

  }
  else {
    // Something is wrong; communicate it and exit
  }

  // Allocate spheres on the host
  unsigned int sphNum = 3;
  struct Sphere *hSpheres =
    (struct Sphere *)calloc(sphNum, sizeof(struct Sphere));

  // Allocate lights on the host
  unsigned int lgtNum = 2;
  struct Light *hLights =
    (struct Light *)calloc(lgtNum, sizeof(struct Light));

  // Setup the default scene
  setupBaseScene(hSpheres, hLights);


  

  // Define the scene
  const unsigned int kScreenWidth = 800;
  const unsigned int kScreenHeight = 600;
  float zoomFactor = -4.f;
  float aliasFactor = 3.f;

  // Define the constants for OCL
  size_t globalWorkSize = kScreenWidth * kScreenHeight;
  size_t localWorkSize = 64;

  // Error code returned from openCL calls
  int err;


  cl_uint numPlatforms;

  // Find number of platforms
  err = clGetPlatformIDs(0, NULL, &numPlatforms);
  checkError(err, "Finding platforms");
  if (numPlatforms == 0) {
    printf("Found 0 platforms!\n");
    return EXIT_FAILURE;
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
  cl_context gpuContext;
  gpuContext = clCreateContext(NULL, 1, &deviceId, NULL, NULL, &err);
  checkError(err, "Creating context");




  // Create a command queue
  cl_command_queue commandsGPU;
  commandsGPU = clCreateCommandQueue(gpuContext, deviceId, NULL, &err);
  checkError(err, "Creating command queue");




  // Load the kernel code
  std::ifstream sourceFstream("raytrace_kernel.cl");
  std::string source((std::istreambuf_iterator<char>(sourceFstream)),
    std::istreambuf_iterator<char>());

  // Create a program from the source
  const char* str = source.c_str();
  cl_program program;
  program = clCreateProgramWithSource(gpuContext, 1, &str, NULL, &err);
  checkError(err, "Creating program");




  // Compile the program
  err = clBuildProgram(program, 0, NULL, "-I C:\Drive\Alberto\Projects\Code\C++\raytracer_gamma\raytracer_gamma", NULL, NULL);
    // If there were compilation errors
    if (err != CL_SUCCESS) {
      // Print out compilation log
      size_t len;
      char buffer[2048];

      printf("Error: Failed to build program executable!\n%s\n", err_code(err));
      clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      printf("%s\n", buffer);

      // Exit
      return EXIT_FAILURE;
    }




  // Create the kernel
  cl_kernel koRTG;
  koRTG = clCreateKernel(program, "raytrace", &err);
  checkError(err, "Creating kernel");




  // Create the list of spheres and lights in device memory
  cl_mem dSpheres = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE,
    sizeof(struct Sphere) * sphNum, NULL, &err);
  checkError(err, "Creating buffer for spheres");
  cl_mem dLights = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE,
    sizeof(struct Light) * lgtNum, NULL, &err);
  checkError(err, "Creating buffer for lights");
  cl_mem dPixelBuffer = clCreateBuffer(gpuContext, CL_MEM_WRITE_ONLY,
    kScreenWidth * kScreenHeight * sizeof(Vec), NULL, &err);
  checkError(err, "Creating buffer for pixels");

  // Write data from host into device memory (fill the buffers with
  // the host arrays)
  err = clEnqueueWriteBuffer(commandsGPU, dSpheres, CL_TRUE, 0,
    sizeof(struct Sphere) * sphNum, hSpheres, 0, NULL, NULL);
  checkError(err, "Copying hSperes in dSpheres");
  err = clEnqueueWriteBuffer(commandsGPU, dLights, CL_TRUE, 0,
    sizeof(struct Light) * lgtNum, hLights, 0, NULL, NULL);
  checkError(err, "Copying hLights into dLights");

  cl_int   status;
  cl_uint maxDims;
  cl_uint recommWorkSize;
  cl_event events[2];
  size_t maxWorkGroupSize;

  /**
  * Query device capabilities. Maximum
  * work item dimensions and the maximmum
  * work item sizes
  */
  err = clGetDeviceInfo(
	  deviceId,
	  CL_DEVICE_MAX_WORK_GROUP_SIZE,
	  sizeof(size_t),
	  (void*)&maxWorkGroupSize,
	  NULL);
  checkError(err, "Getting max work group size");

  status = clGetKernelWorkGroupInfo(
    koRTG,
    deviceId,
    CL_KERNEL_WORK_GROUP_SIZE,
    sizeof(cl_uint),
    &recommWorkSize,
    NULL);
  checkError(err, "Getting recommended work group size");
  printf("Recommended workgroup size for this kernel: %d\n\n", recommWorkSize);


  
  // Set the work item size to be as the one recommended by
  // the implementation depending on the kernel itself
  localWorkSize = recommWorkSize;

  // If the global work size is not a multiple of the local
  // work size
  if (globalWorkSize % localWorkSize != 0) {
    // Pad the global work size to bea multiple of the local
    // work size
    globalWorkSize = (globalWorkSize / localWorkSize + 1) * localWorkSize;
  }



  // Set kernel arguments
  err = clSetKernelArg(koRTG, 0, sizeof(cl_mem), &dSpheres);
  err |= clSetKernelArg(koRTG, 1, sizeof(unsigned int), &sphNum);
  err |= clSetKernelArg(koRTG, 2, sizeof(cl_mem), &dLights);
  err |= clSetKernelArg(koRTG, 3, sizeof(unsigned int), &lgtNum);
  err |= clSetKernelArg(koRTG, 4, sizeof(unsigned int), &kScreenWidth);
  err |= clSetKernelArg(koRTG, 5, sizeof(unsigned int), &kScreenHeight);
  err |= clSetKernelArg(koRTG, 6, sizeof(float), &zoomFactor);
  err |= clSetKernelArg(koRTG, 7, sizeof(float), &aliasFactor);
  err |= clSetKernelArg(koRTG, 8, sizeof(cl_mem), &dPixelBuffer);
  err |= clSetKernelArg(koRTG, 9, sizeof(struct Sphere) * sphNum, NULL);
  err |= clSetKernelArg(koRTG, 10, sizeof(struct Light) * lgtNum, NULL);
  checkError(err, "Setting kernel arguments");

  
  printf("Enqueueing kernel...\t");

  // Start counting the time between kernel enqueuing and completion
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

  // Execute the kernel over the entire range of our 1d input data set
  // letting the OpenCL runtime choose the work-group size
  err = clEnqueueNDRangeKernel(commandsGPU, koRTG, 1, NULL,
    &globalWorkSize, &localWorkSize, 0, NULL, NULL);
  checkError(err, "Enqueueing kernel");

  // Wait for the commands in the queue to be executed
  err = clFinish(commandsGPU);
  checkError(err, "Waiting for commands to finish");

  // Read the time after the kernel has executed
  std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();

  // Compute the duration
  double kernelExecTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

  

  // Print the duration
  printf("Kernel executed!\n");
  printf("Exec time: %.5f ms\n\nSaving PPM...\t", kernelExecTime);



  // Image
  //float *imagePtr = (float *)malloc(globalWorkSize * sizeof(Vec));



  //// Screen in world coordinates
  //const float kImageWorldWidth = 16.f;
  //const float kImageWorldHeight = 12.f;

  //// Amount to increase each step for the ray direction
  //const float kRayXStep = kImageWorldWidth / ((float)kScreenWidth);
  //const float kRayYStep = kImageWorldHeight / ((float)kScreenHeight);
  //const float aspectRatio = kImageWorldWidth / kImageWorldHeight;

  //// Variables holding the current step in world coordinates
  //float rayX = 0.f, rayY = 0.f;

  //int pixelsCounter = 0;

  //// Calculate size of an alias step in world coordinates
  //const float kAliasFactorStepInv = kRayXStep / aliasFactor;
  //// Calculate total size of samples to be taken
  //const float kSamplesTot = aliasFactor * aliasFactor;
  //// Also its inverse
  //const float kSamplesTotinv = 1.f / kSamplesTot;

  //for (int y = 0; y < kScreenWidth * kScreenHeight; ++y, pixelsCounter += 3) {
  //  // Retrieve the global ID of the kernel
  //  const unsigned gid = y;

  //  

  //  // Calculate world position of pixel being currently worked on
  //  const float kPxWorldX = ((((float)(gid % kScreenWidth) - 
  //    (kScreenWidth * 0.5f))) * kRayXStep);
  //  const float kPxWorldY = ((kScreenHeight *0.5f) - ((float)(gid / kScreenWidth))) * kRayYStep;

  //  // The ray to be shot. The vantage point (camera) is at the origin,
  //  // and its intensity is maximum
  //  struct Ray ray; vinit(ray.origin, 0.f, 0.f, 0.f); vinit(ray.intensity, 1.f, 1.f, 1.f);

  //  // The colour of the pixel to be computed
  //  Vec pixelCol = { 0.f, 0.f, 0.f };

  //  // Mock background material
  //  struct Material bgMaterial;
  //  Vec black; vinit(black, 0.f, 0.f, 0.f);
  //  setMatteGlossBalance(&bgMaterial, 0.f, &black, &black);
  //  setMatRefractivityIndex(&bgMaterial, 1.00f);

  //  // For each sample to be taken
  //  for (int i = 0; i < aliasFactor; ++i) {
  //    for (int j = 0; j < aliasFactor; ++j) {
  //      // Calculate the direction of the ray
  //      float x = (kPxWorldX + (float)(((float)j) * kAliasFactorStepInv)) * aspectRatio;
  //      float y = (kPxWorldY + (float)(((float)i) * kAliasFactorStepInv));

  //      // Set the ray's dir and normalise it
  //      vinit(ray.dir, x, y, zoomFactor); vnorm(ray.dir);

  //      // Raytrace for the current sample
  //      Vec currentSampleCol = rayTrace(hSpheres, sphNum, hLights, lgtNum,
  //        ray, bgMaterial, 0);

  //      vsmul(currentSampleCol, kSamplesTotinv, currentSampleCol);

  //      // Compute the average
  //      vadd(pixelCol, pixelCol, currentSampleCol);
  //    }
  //  }

  //  // Write result in destination buffer
  //  *(imagePtr + pixelsCounter) = pixelCol.x;
  //  *(imagePtr + pixelsCounter + 1) = pixelCol.y;
  //  *(imagePtr + pixelsCounter + 2) = pixelCol.z;
  //}

  // Create a buffer to hold the result of the computation on the device
  Vec *pixelsIntermediate = (Vec *)calloc(kScreenHeight * kScreenWidth, sizeof(Vec));
  //Vec *pixelsIntermediate = (Vec *)(imagePtr);

  // Read the results back from the device into the host
  err = clEnqueueReadBuffer(commandsGPU, dPixelBuffer, CL_TRUE, 0,
	  kScreenWidth * kScreenHeight * sizeof(Vec), pixelsIntermediate, 0, NULL, NULL);
  // If the reading operation didn't complete successfully
  if (err != CL_SUCCESS) {
    printf("Error: Failed to read output buffer!\n%s\n", err_code(err));

    // Exit
    exit(1);
  }

  // Calculate the maximum colour value across the whole picture
  float maxColourValue = maxColourValuePixelBuffer(pixelsIntermediate,
    kScreenWidth * kScreenHeight);


  using rtg::RGB;
  // Cast the buffer to the type accepted by the savePPM function
  RGB *pixels = (RGB *)(pixelsIntermediate);
 

  // Cleanup
  clReleaseMemObject(dPixelBuffer);
  clReleaseMemObject(dLights);
  clReleaseMemObject(dSpheres);
  clReleaseProgram(program);
  clReleaseKernel(koRTG);
  clReleaseCommandQueue(commandsGPU);
  clReleaseContext(gpuContext);
  // ... Also on host

  free(hLights);
  free(hSpheres);
  free(platform);


  // Try to save a PPM picture
  rtg::savePPM(pixels, "testPPM.ppm", kScreenWidth, kScreenHeight, maxColourValue);
  free(pixelsIntermediate);
  //free(imagePtr);
 
  printf("Saved PPM!\nEnter a char, then ENTER to exit...\n");

  getchar();

  return 0;
}