
// Includes
#include <ppm.h>
#include <cstdio>
#include <fstream>
#include <algorithm>

namespace rtg {

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

      // Throw errors
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

}
// EO Namespace