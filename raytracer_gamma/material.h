#ifndef _MATERIAL_H
#define _MATERIAL_H


// Includes
#include "vec.h"

struct Material
{
  Vec matteColour;
  Vec glossColour;
  float opacity;
  float refractiveIndex;
};

#endif // !_MATERIAL_H