/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 29-Nov-2021 17:59:00
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "Features2CostMatrix.h"
#include "main.h"
#include "Features2CostMatrix_terminate.h"
#include "Features2CostMatrix_emxAPI.h"
#include "Features2CostMatrix_initialize.h"

/* Function Declarations */
static int argInit_int32_T(void);
static double argInit_real_T(void);
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void);
static void main_Features2CostMatrix(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : int
 */
static int argInit_int32_T(void)
{
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : emxArray_real_T *
 */
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void)
{
  emxArray_real_T *result;
  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreate_real_T(2, 2);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_Features2CostMatrix(void)
{
  emxArray_real_T *C;
  emxArray_real_T *F1;
  emxArray_real_T *F2;
  emxArray_real_T *Points1;
  emxArray_real_T *Points2;
  emxInitArray_real_T(&C, 2);

  /* Initialize function 'Features2CostMatrix' input arguments. */
  /* Initialize function input argument 'F1'. */
  F1 = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'F2'. */
  F2 = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'Points1'. */
  Points1 = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'Points2'. */
  Points2 = c_argInit_UnboundedxUnbounded_r();

  /* Call the entry-point 'Features2CostMatrix'. */
  Features2CostMatrix(F1, F2, Points1, Points2, argInit_int32_T(), C);
  emxDestroyArray_real_T(C);
  emxDestroyArray_real_T(Points2);
  emxDestroyArray_real_T(Points1);
  emxDestroyArray_real_T(F2);
  emxDestroyArray_real_T(F1);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  Features2CostMatrix_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_Features2CostMatrix();

  /* Terminate the application.
     You do not need to do this more than one time. */
  Features2CostMatrix_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
