#ifndef _RAYTRACE_STACK_H
#define _RAYTRACE_STACK_H

// Includes
#include "vec.h"
#include "ray.h"
#include "material.h"
#include "intersection.h"

#define RTSTACK_MAXSIZE 6

// RtSnapshot
typedef struct
{
  struct Ray ray;
  int traceDepth;
  int stage;
  Vec colour;
  struct Intersection intersection;
  struct Material refractiveMat;
  float refractiveReflectionFactor;
} RtSnapshot;

// RtStack
typedef struct 
{
  // Elements contained by the stack
  RtSnapshot elements[RTSTACK_MAXSIZE];
  
  // Top element in the stack
  int top;

  // Size of the stack
  int maxSize;
} RtStack;

// Initialise the stack
void rtStackInit(RtStack *sP);

bool rtStackIsEmpty(RtStack *sP);

bool rtStackIsFull(RtStack *sP);

void rtStackPush(RtStack *sP, RtSnapshot *element);

// Pop an element from the stack
void rtStackPop(RtStack *sP);

// Return a reference to the top value of the stack
RtSnapshot *rtStackTop(RtStack *sP);

#endif