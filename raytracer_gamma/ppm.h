#ifndef _PPM_H
#define _PPM_H


namespace rtg {

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


  // Function which saves the given pixels buffer to a
  // PPM file
  void savePPM(const RGB *pixels,
    const char *filename,
    const int width,
    const int height,
    const float maxColourVal);

}
// EO Namespace


#endif