/*
 * _coder_triangular_algorithm_TBM_mex.c
 *
 * Code generation for function 'triangular_algorithm_TBM'
 *
 */

/* Include files */
#include "mex.h"
#include "_coder_triangular_algorithm_TBM_api.h"
#include "triangular_algorithm_TBM_initialize.h"
#include "triangular_algorithm_TBM_terminate.h"

/* Function Declarations */
static void triangular_algorithm_TBM_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "triangular_algorithm_TBM", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
void *emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void triangular_algorithm_TBM_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *outputs[3];
  const mxArray *inputs[5];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  emlrtStack st = { NULL, NULL, NULL };
  /* Module initialization. */
  triangular_algorithm_TBM_initialize(&emlrtContextGlobal);
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 5, mxCHAR_CLASS, 24, "triangular_algorithm_TBM");
  } else if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 24, "triangular_algorithm_TBM");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = prhs[n];
  }
  /* Call the function. */
  triangular_algorithm_TBM_api(inputs, outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  triangular_algorithm_TBM_terminate();
}

void triangular_algorithm_TBM_atexit_wrapper(void)
{
   triangular_algorithm_TBM_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(triangular_algorithm_TBM_atexit_wrapper);
  /* Dispatch the entry-point. */
  triangular_algorithm_TBM_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (_coder_triangular_algorithm_TBM_mex.c) */
