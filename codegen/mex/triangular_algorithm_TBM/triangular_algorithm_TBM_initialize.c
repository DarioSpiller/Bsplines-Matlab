/*
 * triangular_algorithm_TBM_initialize.c
 *
 * Code generation for function 'triangular_algorithm_TBM_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "triangular_algorithm_TBM.h"
#include "triangular_algorithm_TBM_initialize.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar;

/* Function Definitions */
void triangular_algorithm_TBM_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (triangular_algorithm_TBM_initialize.c) */
