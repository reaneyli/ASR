/*
 * File: _coder_Features2CostMatrix_or_mex.cpp
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 29-Nov-2021 18:46:19
 */

/* Include Files */
#include "_coder_Features2CostMatrix_or_api.h"
#include "_coder_Features2CostMatrix_or_mex.h"

/* Function Declarations */
static void c_Features2CostMatrix_or_mexFun(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[5]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[1]
 *                int32_T nrhs
 *                const mxArray *prhs[5]
 * Return Type  : void
 */
static void c_Features2CostMatrix_or_mexFun(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[5])
{
  const mxArray *outputs[1];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        22, "Features2CostMatrix_or");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "Features2CostMatrix_or");
  }

  /* Call the function. */
  Features2CostMatrix_or_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, plhs, outputs);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(Features2CostMatrix_or_atexit);

  /* Module initialization. */
  Features2CostMatrix_or_initialize();

  /* Dispatch the entry-point. */
  c_Features2CostMatrix_or_mexFun(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  Features2CostMatrix_or_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_Features2CostMatrix_or_mex.cpp
 *
 * [EOF]
 */
