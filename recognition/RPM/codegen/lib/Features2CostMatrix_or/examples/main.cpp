//
// File: main.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 29-Nov-2021 18:46:19
//

//***********************************************************************
// This automatically generated example C main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************
// Include Files
#include "rt_nonfinite.h"
#include "Features2CostMatrix_or.h"
#include "main.h"
#include "Features2CostMatrix_or_terminate.h"
#include "Features2CostMatrix_or_emxAPI.h"
#include "Features2CostMatrix_or_initialize.h"

// Function Declarations
static double argInit_real_T();
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r();
static void main_Features2CostMatrix_or();

// Function Definitions

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : void
// Return Type  : emxArray_real_T *
//
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r()
{
  emxArray_real_T *result;
  int idx0;
  int idx1;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreate_real_T(2, 2);

  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_Features2CostMatrix_or()
{
  emxArray_real_T *C;
  emxArray_real_T *F1;
  emxArray_real_T *F2;
  emxArray_real_T *Points1;
  emxArray_real_T *Points2;
  emxInitArray_real_T(&C, 2);

  // Initialize function 'Features2CostMatrix_or' input arguments.
  // Initialize function input argument 'F1'.
  F1 = c_argInit_UnboundedxUnbounded_r();

  // Initialize function input argument 'F2'.
  F2 = c_argInit_UnboundedxUnbounded_r();

  // Initialize function input argument 'Points1'.
  Points1 = c_argInit_UnboundedxUnbounded_r();

  // Initialize function input argument 'Points2'.
  Points2 = c_argInit_UnboundedxUnbounded_r();

  // Call the entry-point 'Features2CostMatrix_or'.
  Features2CostMatrix_or(F1, F2, Points1, Points2, argInit_real_T(), C);
  emxDestroyArray_real_T(C);
  emxDestroyArray_real_T(Points2);
  emxDestroyArray_real_T(Points1);
  emxDestroyArray_real_T(F2);
  emxDestroyArray_real_T(F1);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // Initialize the application.
  // You do not need to do this more than one time.
  Features2CostMatrix_or_initialize();

  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_Features2CostMatrix_or();

  // Terminate the application.
  // You do not need to do this more than one time.
  Features2CostMatrix_or_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
