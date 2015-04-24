#ifndef _RAYTRACE_STACK_H
#define _RAYTRACE_STACK_H

// Includes
#include "vec.h"
#include "ray.h"
#include "material.h"
#include "intersection.h"

#define RTSTACK_MAXSIZE 11

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
void rtStackInit(RtStack *sP) {
  // Empty
  sP->top = -1;
  sP->maxSize = RTSTACK_MAXSIZE;
}

bool rtStackIsEmpty(RtStack *sP) {
  return (sP->top < 0);
}

bool rtStackIsFull(RtStack *sP) {
  return (sP->top >= sP->maxSize - 1);
}

void rtStackPush(RtStack *sP, RtSnapshot *element) {
  // If the stack is not full
  if (!rtStackIsFull(sP)) {
    // Push element into the stack
    sP->elements[++sP->top] = *element;
  }
}

// Pop an element from the stack
void rtStackPop(RtStack *sP) {
  --sP->top;
}

// Return a reference to the top value of the stack
RtSnapshot *rtStackTop(RtStack *sP) {
  return &sP->elements[sP->top];
}

#endif