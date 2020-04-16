/*
 * triangular_algorithm_TBM.h
 *
 * Code generation for function 'triangular_algorithm_TBM'
 *
 */

#ifndef __TRIANGULAR_ALGORITHM_TBM_H__
#define __TRIANGULAR_ALGORITHM_TBM_H__

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
extern void triangular_algorithm_TBM(int8_T n, int8_T p, real_T u, int8_T m,
  const emxArray_real_T *knot_points, real_T N_0_data[], int32_T N_0_size[2],
  real_T N_1_data[], int32_T N_1_size[2], real_T N_2_data[], int32_T N_2_size[2]);

#endif

/* End of code generation (triangular_algorithm_TBM.h) */
