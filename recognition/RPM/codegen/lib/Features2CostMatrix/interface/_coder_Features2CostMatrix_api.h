/*
 * File: _coder_Features2CostMatrix_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 29-Nov-2021 17:59:00
 */

#ifndef _CODER_FEATURES2COSTMATRIX_API_H
#define _CODER_FEATURES2COSTMATRIX_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Features2CostMatrix_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void Features2CostMatrix(emxArray_real_T *F1, emxArray_real_T *F2,
  emxArray_real_T *Points1, emxArray_real_T *Points2, int32_T maxdist,
  emxArray_real_T *C);
extern void Features2CostMatrix_api(const mxArray * const prhs[5], int32_T nlhs,
  const mxArray *plhs[1]);
extern void Features2CostMatrix_atexit(void);
extern void Features2CostMatrix_initialize(void);
extern void Features2CostMatrix_terminate(void);
extern void Features2CostMatrix_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Features2CostMatrix_api.h
 *
 * [EOF]
 */
