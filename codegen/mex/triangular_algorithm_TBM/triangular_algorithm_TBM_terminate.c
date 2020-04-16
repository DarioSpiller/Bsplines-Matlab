/*
 * triangular_algorithm_TBM_terminate.c
 *
 * Code generation for function 'triangular_algorithm_TBM_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "triangular_algorithm_TBM.h"
#include "triangular_algorithm_TBM_terminate.h"

/* Function Definitions */
void triangular_algorithm_TBM_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void triangular_algorithm_TBM_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (triangular_algorithm_TBM_terminate.c) */
