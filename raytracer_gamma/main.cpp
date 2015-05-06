// Hello world with OpenCL


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <device_info.h>
#include <err_code.h>
#include <fstream>
#include <algorithm>
#include <cpu_raytracer.h>
#include <chrono>
#include <ppm.h>
#include <string>
#include <iostream>
#include <raytracer.h>
#include <gpu_raytracer.h>
#include <cpu_raytracer.h>


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif

#define LENGTH (1024)    // length of vectors a, b, and c


// Colours
const Vec whiteCol = { 0.8f, 0.8f, 0.8f };
const Vec lowerWhite = { 0.5f, 0.5f, 0.5f };
const Vec redCol = { 0.8f, 0.0f, 0.0f };
const Vec greenCol = { 0.0f, 0.8f, 0.0f };
const Vec blueCol = { 0.0f, 0.0f, 0.8f };
const Vec goldCol = { 1.f, 0.843137f, 0.0f };

// Setup the default scene
void setupBaseScene(struct Sphere *spheres, struct Light *lights) {
  // Setup materials
  struct Material ballMaterial1; // Green
  Vec bm1Gloss; vassign(bm1Gloss, greenCol);
  Vec bm1Matte; vassign(bm1Matte, greenCol);
  setMatOpacity(&ballMaterial1, 0.35f);
  setMatteGlossBalance(&ballMaterial1, 0.8f, &bm1Matte, &bm1Gloss);
  setMatRefractivityIndex(&ballMaterial1, 1.5500f);

  struct Material ballMaterial2; // Red
  Vec bm2Gloss; vassign(bm2Gloss, redCol);
  Vec bm2Matte; vassign(bm2Matte, redCol);
  setMatOpacity(&ballMaterial2, 0.9f);
  setMatteGlossBalance(&ballMaterial2, 0.01f, &bm2Matte, &bm2Gloss);
  setMatRefractivityIndex(&ballMaterial2, 1.5500f);

  struct Material ballMaterial3; // White
  Vec bm3Gloss; vassign(bm3Gloss, whiteCol);
  Vec bm3Matte; vassign(bm3Matte, whiteCol);
  setMatOpacity(&ballMaterial3, 0.6f);
  setMatteGlossBalance(&ballMaterial3, 0.9f, &bm3Matte, &bm3Gloss);
  setMatRefractivityIndex(&ballMaterial3, 1.5500f);

  struct Material ballMaterial4; // Blue
  Vec bm4Gloss; vassign(bm4Gloss, goldCol);
  Vec bm4Matte; vassign(bm4Matte, blueCol);
  setMatOpacity(&ballMaterial4, 0.8f);
  setMatteGlossBalance(&ballMaterial4, 0.9f, &bm4Matte, &bm4Gloss);
  setMatRefractivityIndex(&ballMaterial4, 1.5500f);

  // Setup spheres
  spheres[0].material = ballMaterial2;
  vinit(spheres[0].pos, -10.f, 0.f, -13.f);
  spheres[0].radius = 5.f;
  spheres[1].material = ballMaterial1;
  vinit(spheres[1].pos, -3.f, 1.5f, -5.f);
  spheres[1].radius = 2.f;
  spheres[2].material = ballMaterial3;
  vinit(spheres[2].pos, 2.f, -1.f, -14.f);
  spheres[2].radius = 5.f;
  spheres[3].material = ballMaterial4;
  vinit(spheres[3].pos, 0.f, -11.f, -13.5f);
  spheres[3].radius = 4.5f;

  // Setup light sources
  vinit(lights[0].pos, -45.f, 10.f, 85.f);
  vassign(lights[0].col, whiteCol);
  vinit(lights[1].pos, -45.f, 1.f, 85.f);
  vassign(lights[1].col, whiteCol);
  vinit(lights[2].pos, 45.f, 5.f, 70.f);
  vassign(lights[2].col, whiteCol);
}


int main(int argc, char** argv)
{

  

  // Define the scene
  const unsigned int kImgWidth = 640;
  const unsigned int kImgHeight = 480;
  float zoomFactor = -4.f;
  float aliasFactor = 1.f;

  rtg::Raytracer *raytracer = NULL;

  // Read the parameters
  if (argc == 1) {
    // Use default settings
    raytracer =
      new rtg::CPURaytracer(
      kImgWidth,
      kImgHeight,
      aliasFactor);
  }
  else if (argc == 2) {
    // Read the parameter

    // Read the first parameter
    std::string first = argv[1];

    // Depending on its type
    if (first == "-cpu") {
      raytracer =
        new rtg::CPURaytracer(
        kImgWidth,
        kImgHeight,
        aliasFactor);
    }
  }
  else if (argc == 3) {
    // Read the passed parameters

    // Read the first parameter
    std::string first = argv[1];

    // Depending on its type
    if (first == "-ocl") {
      // Read the second parameter
      std::string second = argv[2];

      // Depending on its type

      // If default mode
      if (second == "-d") {
        raytracer =
          new rtg::GPURaytracer(
          kImgWidth, 
          kImgHeight, 
          aliasFactor,
          rtg::kModeDefault);
      }
      // If one line per work item mode
      else if (second == "-l") {
        raytracer =
          new rtg::GPURaytracer(
          kImgWidth,
          kImgHeight,
          aliasFactor,
          rtg::kModeLine);
      }
    }
    else if (first == "-cpu") {
      // Read the second parameter
      std::string second = argv[2];

      // Depending on its type

      // If default mode
      if (second == "-d") {
        raytracer =
          new rtg::CPURaytracer(
          kImgWidth,
          kImgHeight,
          aliasFactor);
      }

    }

    

  }
  else {
    // Something is wrong; communicate it and exit
  }

  // Allocate spheres on the host
  unsigned int sphNum = 4;
  struct Sphere *hSpheres =
    (struct Sphere *)calloc(sphNum, sizeof(struct Sphere));

  // Allocate lights on the host
  unsigned int lgtNum = 3;
  struct Light *hLights =
    (struct Light *)calloc(lgtNum, sizeof(struct Light));

  // Setup the default scene
  setupBaseScene(hSpheres, hLights);


  // Setup the raytracer. This will also print out information 
  // if in OCL mode
  raytracer->setup(hSpheres, sphNum, hLights, lgtNum);

  
  
  printf("Launching raytracing...\n");

  // Start counting the time between raytracing start and completion
  std::chrono::steady_clock::time_point startTime = 
    std::chrono::steady_clock::now();

  raytracer->raytrace();

  // Read the time after the raytracer has executed
  std::chrono::steady_clock::time_point endTime = 
    std::chrono::steady_clock::now();

  // Compute the duration
  double raytraceExecTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

  

  // Print the duration
  printf(" Raytracing completed!\n");
  printf("Exec time: %.5f ms\n\nSaving PPM...\t", raytraceExecTime);



  // Create a buffer to hold the result of the computation on the device
  Vec *pixelsIntermediate = (Vec *)calloc(kImgHeight * kImgWidth, sizeof(Vec));
  //Vec *pixelsIntermediate = (Vec *)(imagePtr);

  raytracer->readResult(pixelsIntermediate);

  // Calculate the maximum colour value across the whole picture
  float maxColourValue = maxColourValuePixelBuffer(pixelsIntermediate,
    kImgWidth * kImgHeight);


  using rtg::RGB;
  // Cast the buffer to the type accepted by the savePPM function
  RGB *pixels = (RGB *)(pixelsIntermediate);
 

  // Cleanup
  // ... Also on host

  free(hLights);
  free(hSpheres);
  delete raytracer;


  // Try to save a PPM picture
  rtg::savePPM(pixels, "testPPM.ppm", kImgWidth, kImgHeight, maxColourValue);
  free(pixelsIntermediate);
  //free(imagePtr);
 
  printf("Saved PPM!\nEnter a char, then ENTER to exit...\n");

  getchar();

  return 0;
}