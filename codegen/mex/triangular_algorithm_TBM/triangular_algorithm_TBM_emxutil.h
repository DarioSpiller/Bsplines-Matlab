/*
 * triangular_algorithm_TBM_emxutil.h
 *
 * Code generation for function 'triangular_algorithm_TBM_emxutil'
 *
 */

#ifndef __TRIANGULAR_ALGORITHM_TBM_EMXUTIL_H__
#define __TRIANGULAR_ALGORITHM_TBM_EMXUTIL_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "triangular_algorithm_TBM_types.h"

/* Function Declarations */
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions,
  boolean_T doPush);

#endif

/* End of code generation (triangular_algorithm_TBM_emxutil.h) */
