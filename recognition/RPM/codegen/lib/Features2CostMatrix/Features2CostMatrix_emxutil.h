/*
 * File: Features2CostMatrix_emxutil.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 29-Nov-2021 17:59:00
 */

#ifndef FEATURES2COSTMATRIX_EMXUTIL_H
#define FEATURES2COSTMATRIX_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "Features2CostMatrix_types.h"

/* Function Declarations */
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for Features2CostMatrix_emxutil.h
 *
 * [EOF]
 */
