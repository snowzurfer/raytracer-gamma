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


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif

#define LENGTH (1024)    // length of vectors a, b, and c

//extern double wtime();       // returns time since some fixed past point (wtime.c)
//extern int output_device_info(cl_device_id);


// Structure representing a pixel
struct RGB
{
  RGB() : r(0), g(0), b(0)  {}
  RGB(float c) : r(c), g(c), b(c) {}
  RGB(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {}
  bool operator != (const RGB &c) const { return c.r != r || c.g != g || c.b != b; }
  RGB& operator *= (const RGB &rgb) { r *= rgb.r, g *= rgb.g, b *= rgb.b; return *this; }
  RGB& operator += (const RGB &rgb) { r += rgb.r, g += rgb.g, b += rgb.b; return *this; }

  // RGB Values
  float r, g, b;
};


void savePPM(const RGB *pixels, 
  const char *filename, 
  const int width, 
  const int height,
  const float maxColourVal) {
  // Check for wrong parameters
  if (width == 0 || height == 0) {
    fprintf(stderr, "Can't save an empty image\n"); 
    return; 
  }

  // Create a filestream to write the image
  std::ofstream ofs;
  try {
    // Open/create the file
    ofs.open(filename, std::ios::binary);

    // Catch errors
    if (ofs.fail()) {
      throw("Can't open output file");
    }
    
    // Write headers of PPM
    ofs << "P6\n" << width << " " << height << "\n255\n";
    // Create temp variables to store 8-bit values of each color component
    unsigned char r, g, b;
    // Loop over each pixel in the image, clamp it and convert to byte format
    for (int i = 0; i < width * height; ++i) {
      r = static_cast<unsigned char>
        (std::min(1.f, pixels[i].r) * 255 / maxColourVal);
      g = static_cast<unsigned char>
        (std::min(1.f, pixels[i].g) * 255 / maxColourVal);
      b = static_cast<unsigned char>
        (std::min(1.f, pixels[i].b) * 255 / maxColourVal);


      // Write the values
      ofs << r << g << b;
    }

    // Close the file
    ofs.close();
  }
  // If there were errors, catch them
  catch (const char *err) {
    fprintf(stderr, "%s\n", err);
    ofs.close();
  }
};


int main(int argc, char** argv)
{
  // Error code returned from openCL calls
  int err;

  // A, B and C arrays
  float *hA = (float *)calloc(LENGTH, sizeof(float));
  float *hB = (float *)calloc(LENGTH, sizeof(float));
  float *hC = (float *)calloc(LENGTH, sizeof(float));

  // Define the scene
  const unsigned int kScreenWidth = 800;
  const unsigned int kScreenHeight = 600;
  float zoomFactor = -4.f;
  float aliasFactor = 3.f;

  size_t globalWorkSize = kScreenWidth * kScreenHeight;
  size_t localWorkSize = 256;

  // Colours
  Vec whiteCol;
  vinit(whiteCol, 8.f, 8.f, 8.f);
  Vec lowerWhite;
  vinit(lowerWhite, 0.5f, 0.5f, 0.5f);
  Vec redCol;
  vinit(redCol, 0.8f, 1.f, 0.7f);
  Vec greenCol;
  vinit(greenCol, 0.4f, 0.5f, 0.7f);
  Vec col1;
  vinit(col1, 0.01f, 0.8f, 0.01f);

  // Setup materials
  struct Material ballMaterial1; // White
  Vec bm1Gloss; vassign(bm1Gloss, redCol);
  Vec bm1Matte; vassign(bm1Matte, greenCol);
  setMatOpacity(&ballMaterial1, 0.8f);
  setMatteGlossBalance(&ballMaterial1, 0.2f, &bm1Matte, &bm1Gloss);
  setMatRefractivityIndex(&ballMaterial1, 1.5500f);

  struct Material ballMaterial2; // Red
  Vec bm2Gloss; vassign(bm2Gloss, redCol);
  Vec bm2Matte; vassign(bm2Matte, greenCol);
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
  unsigned int sphNum = 3;
  struct Sphere *hSpheres =
    (struct Sphere *)calloc(sphNum, sizeof(struct Sphere));
  hSpheres[0].material = ballMaterial1;
  vinit(hSpheres[0].pos, -9.f, 0.f, -13.f);
  hSpheres[0].radius = 5.f;
  hSpheres[1].material = ballMaterial2;
  vinit(hSpheres[1].pos, -4.f, 1.5f, -5.f);
  hSpheres[1].radius = 2.f;
  hSpheres[2].material = ballMaterial3;
  vinit(hSpheres[2].pos, 1.f, -1.f, -7.f);
  hSpheres[2].radius = 3.f;

  // Setup light sources
  unsigned int lgtNum = 2;
  struct Light *hLights =
    (struct Light *)calloc(lgtNum, sizeof(struct Light));
  vinit(hLights[0].pos, -45.f, 10.f, 85.f);
  vassign(hLights[0].col, lowerWhite);
  vinit(hLights[1].pos, 20.f, 60.f, -5.f);
  vassign(hLights[1].col, lowerWhite);



  // Fill vectors a and b with random float values
  int count = LENGTH;
  for (int i = 0; i < count; i++){
    hA[i] = rand() / (float)RAND_MAX;
    hB[i] = rand() / (float)RAND_MAX;
  }




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
  cl_event events[2];
  size_t maxWorkGroupSize;

  /**
  * Query device capabilities. Maximum
  * work item dimensions and the maximmum
  * work item sizes
  */
  status = clGetDeviceInfo(
	  deviceId,
	  CL_DEVICE_MAX_WORK_GROUP_SIZE,
	  sizeof(size_t),
	  (void*)&maxWorkGroupSize,
	  NULL);
  if (status != CL_SUCCESS)
  {
	  fprintf(stderr, "Error: Getting Device Info. (clGetDeviceInfo)\n");
	  return 1;
  }

  status = clGetDeviceInfo(
	  deviceId,
	  CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
	  sizeof(cl_uint),
	  (void*)&maxDims,
	  NULL);
  if (status != CL_SUCCESS)
  {
	  fprintf(stderr, "Error: Getting Device Info. (clGetDeviceInfo)\n");
	  return 1;
  }

  localWorkSize = maxWorkGroupSize;

  if (globalWorkSize % localWorkSize != 0) {
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
  printf("Exec time: %.5f ms", kernelExecTime);



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

  // Cast the buffer to the type accepted by the savePPM function
  RGB *pixels = (RGB *)(pixelsIntermediate);

  // Print execution time
  /*rtime = wtime() - rtime;
  printf("\nThe kernel ran in %lf seconds\n", rtime);*/
 

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
  free(hA);
  free(hB);
  free(hC);
  free(platform);


  // Try to save a PPM picture
  savePPM(pixels, "testPPM.ppm", kScreenWidth, kScreenHeight, maxColourValue);
  free(pixelsIntermediate);
  //free(imagePtr);
 
  getchar();

  return 0;
}