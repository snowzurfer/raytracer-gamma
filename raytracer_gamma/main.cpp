// Hello world with OpenCL


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "err_code.h"
#include <fstream>
#include <algorithm>
#include <raytracer.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif

#define TOL    (0.001)   // tolerance used in floating point comparisons
#define LENGTH (1024)    // length of vectors a, b, and c

//extern double wtime();       // returns time since some fixed past point (wtime.c)
extern int output_device_info(cl_device_id);


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


void savePPM(const RGB *pixels, const char *filename, const int width, const int height) {
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
      r = static_cast<unsigned char>(std::min(1.f, pixels[i].r) * 255);
      g = static_cast<unsigned char>(std::min(1.f, pixels[i].g) * 255);
      b = static_cast<unsigned char>(std::min(1.f, pixels[i].b) * 255);

      if (r > 1) {

        printf("red: %d \n", r);
      }

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
  float zoomFactor = -1.5f;
  float aliasFactor = 1.f;

  size_t globalWorkSize = kScreenWidth * kScreenHeight;

  // Colours
  Vec whiteCol;
  vinit(whiteCol, 1.f, 1.f, 1.f);
 
  // Setup materials
  struct Material ballMaterial1; // White
  Vec bm1Gloss; vassign(bm1Gloss, whiteCol);
  Vec bm1Matte; vassign(bm1Matte, whiteCol);
  setMatOpacity(&ballMaterial1, 1.f);
  setMatteGlossBalance(&ballMaterial1, 0.6f, &bm1Matte, &bm1Gloss);

  // Setup spheres
  unsigned int sphNum = 1;
  struct Sphere *hSpheres = 
    (struct Sphere *)calloc(sphNum, sizeof(struct Sphere));
  hSpheres[0].material = ballMaterial1;
  vinit(hSpheres[0].pos, 0.f, 0.f, -10.f);
  hSpheres[0].radius = 3.f;

  // Setup light sources
  unsigned int lgtNum = 1;
  struct Light *hLights = 
    (struct Light *)calloc(lgtNum, sizeof(struct Light));
  vinit(hLights[0].pos, 0.f, 6.f, -4.f);
  vassign(hLights[0].col, whiteCol);




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
  //err = output_device_info(deviceId);
  //checkError(err, "Printing device output");



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
  err = clBuildProgram(program, 0, NULL, "-I C:\Drive\Alberto\Projects\Code\C++\raytracer_gamma\raytracer_gamma ", NULL, NULL);
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
  err = clEnqueueWriteBuffer(commandsGPU, dSpheres, CL_FALSE, 0,
    sizeof(struct Sphere) * sphNum, hSpheres, 0, NULL, NULL);
  checkError(err, "Copying hSperes in dSpheres");
  err = clEnqueueWriteBuffer(commandsGPU, dLights, CL_FALSE, 0,
    sizeof(struct Light) * lgtNum, hLights, 0, NULL, NULL);
  checkError(err, "Copying hLights into dLights");



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
  checkError(err, "Setting kernel arguments");

  /*double rtime = wtime();*/

  // Execute the kernel over the entire range of our 1d input data set
  // letting the OpenCL runtime choose the work-group size
  err = clEnqueueNDRangeKernel(commandsGPU, koRTG, 1, NULL,
    &globalWorkSize, NULL, 0, NULL, NULL);
  checkError(err, "Enqueueing kernel");

  // Wait for the commands in the queue to be executed
  err = clFinish(commandsGPU);
  checkError(err, "Waiting for commands to finish");

  printf("Size of vec3: %d", sizeof(Vec));



  // Image
  float *imagePtr = (float *)malloc(globalWorkSize * sizeof(Vec));



  // Screen in world coordinates
  const float kImageWorldWidth = 16.f;
  const float kImageWorldHeight = 12.f;

  // Amount to increase each step for the ray direction
  const float kRayXStep = kImageWorldWidth / ((float)kScreenWidth);
  const float kRayYStep = kImageWorldHeight / ((float)kScreenHeight);

  // Variables holding the current step in world coordinates
  //float rayX = 0.f, rayY = 0.f;

  /*int pixelsCounter = 0;

  for (int y = 0; y < kScreenWidth * kScreenHeight; ++y, pixelsCounter += 3) {
      // Retrieve the global ID of the kernel
      const unsigned gid = y;

      // Calculate inverse of aliasFactor
      const float kAliasFactorInv = 1.f / aliasFactor;
      // Calculate total size of samples to be taken
      const float kSamplesTot = aliasFactor * aliasFactor;
      // Also its inverse
      const cl_float kSamplesTotinv = 1.f / kSamplesTot;

      // Calculate world position of pixel being currently worked on
      const float kPxWorldX = (((float)(gid % kScreenWidth) - (kScreenWidth * 0.5f))) * kRayXStep;
      const float kPxWorldY = ((kScreenHeight *0.5f) - ((float)(gid / kScreenWidth))) * kRayYStep;

      // The ray to be shot. The vantage point (camera) is at the origin,
      // and its intensity is maximum
      struct Ray ray; vinit(ray.origin, 0.f, 0.f, 0.f); vinit(ray.intensity, 1.f, 1.f, 1.f);

      // The colour of the pixel to be computed
      Vec pixelCol = { 0.f, 0.f, 0.f };

      // Mock background material
      struct Material bgMaterial;
      Vec black; vinit(black, 0.f, 0.f, 0.f);
      setMatteGlossBalance(&bgMaterial, 0.f, &black, &black);

      // For each sample to be taken
      for (int i = 0; i < aliasFactor; ++i) {
        for (int j = 0; j < aliasFactor; ++j) {
          // Calculate the direction of the ray
          float x = kPxWorldX + (float)(((float)j) * kAliasFactorInv);
          float y = kPxWorldY + (float)(((float)i) * kAliasFactorInv);

          // Set the ray's dir and normalise it
          vinit(ray.dir, x, y, zoomFactor); vnorm(ray.dir);

          // Raytrace for the current sample
          Vec currentSampleCol = rayTrace(hSpheres, sphNum, hLights, lgtNum,
            &ray, &bgMaterial, 0);

          vsmul(currentSampleCol, kSamplesTotinv, currentSampleCol);

          // Compute the average
          vadd(pixelCol, pixelCol, currentSampleCol);
        }
      }

      // Write result in destination buffer
      *(imagePtr + pixelsCounter)  = pixelCol.x;
      *(imagePtr + pixelsCounter + 1) = pixelCol.y;
      *(imagePtr + pixelsCounter + 2) = pixelCol.z;
    }*/
  







  // Print execution time
  /*rtime = wtime() - rtime;
  printf("\nThe kernel ran in %lf seconds\n", rtime);*/

  // Read back the results from the device memory
  // Create a buffer of pixels
  

  err = clEnqueueReadBuffer(commandsGPU, dPixelBuffer, CL_TRUE, 0,
    kScreenWidth * kScreenHeight * sizeof(Vec), imagePtr, 0, NULL, NULL);
  // If the reading operation didn't complete successfully
  if (err != CL_SUCCESS) {
    printf("Error: Failed to read output array!\n%s\n", err_code(err));

    // Exit
    exit(1);
  }

  // Test the results

  

  // Summarise results
 /* printf("C = A+B:  %d out of %d results were correct.\n", correctResNum,
    count);*/

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

  fprintf(stdout, "Size of RGB: %d \n", sizeof(RGB));

  RGB *pixels = (RGB *)calloc(kScreenHeight * kScreenWidth, sizeof(RGB));
  
  memcpy(pixels, imagePtr, (kScreenHeight * kScreenWidth * sizeof(RGB)));


  float temp1;

  int correctResNum = 0;

  printf("tot: %d \n", globalWorkSize * sizeof(Vec));

  /*for (int i = 0; i < (globalWorkSize * 12); i += 3) {
    temp1 = imagePtr[i];
    float temp2 = imagePtr[1 + i];
    float temp3 = imagePtr[i + 2];

    // assign element i of a+b to tmp
    // compute deviation of expected and output result
    if (temp1 > 0.f || temp2 > 0.f || temp3 > 0.f) {  // correct if square deviation is less than tolerance squared
      correctResNum++;
    }
    else if (temp1 != 0.f && temp2 != 0.f && temp3 != 0.f) {
      printf(" temp1 %f temp2 %f temp3 %f \n", temp1,
        temp2, temp3);
    }
  }*/

  printf("correct resz: %d \n", correctResNum);

  RGB temp;
  for (int i = 0; i < globalWorkSize; i++) {
    temp = pixels[i];     // assign element i of a+b to tmp
                            // compute deviation of expected and output result
    if (temp.b > 0.f || temp.r > 0.f || temp.g > 0.f) {  // correct if square deviation is less than tolerance squared
      correctResNum++;
      printf(" temp.r %f temp.g %f temp.b %f \n", temp.r,
        temp.g, temp.b);
    }
    else {
      //printf(" temp.r %f temp.g %f temp.b %f \n", temp.r,
        //temp.g, temp.b);
    }
  }

  free(imagePtr);

  //int imageCounter
  //for (int i = 0; i < (kScreenHeight * kScreenWidth); i++) {
  //  pixels[i].r = imagePtr[i]; /* red */
  //  pixels[i].g = imagePtr[i + 1];  /* green */
  //  pixels[i].b = (float)((i % 256) / 256.f);  /* blue */
  //}

  savePPM(pixels, "testPPM.ppm", kScreenWidth, kScreenHeight);
  free(pixels);

 

  return 0;
}