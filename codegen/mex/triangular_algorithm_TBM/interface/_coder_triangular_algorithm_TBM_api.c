/*
 * _coder_triangular_algorithm_TBM_api.c
 *
 * Code generation for function '_coder_triangular_algorithm_TBM_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "triangular_algorithm_TBM.h"
#include "_coder_triangular_algorithm_TBM_api.h"
#include "triangular_algorithm_TBM_emxutil.h"

/* Function Declarations */
static int8_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static const mxArray *b_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2]);
static real_T c_emlrt_marshallIn(const mxArray *u, const char_T *identifier);
static const mxArray *c_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2]);
static real_T d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static void e_emlrt_marshallIn(const mxArray *knot_points, const char_T
  *identifier, emxArray_real_T *y);
static int8_T emlrt_marshallIn(const mxArray *n, const char_T *identifier);
static const mxArray *emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2]);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static int8_T g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static real_T h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);

/* Function Definitions */
static int8_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  int8_T y;
  y = g_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  static const int32_T iv1[2] = { 0, 0 };

  const mxArray *m1;
  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv1, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m1, (void *)u_data);
  emlrtSetDimensions((mxArray *)m1, u_size, 2);
  emlrtAssign(&y, m1);
  return y;
}

static real_T c_emlrt_marshallIn(const mxArray *u, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(emlrtAlias(u), &thisId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *c_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  static const int32_T iv2[2] = { 0, 0 };

  const mxArray *m2;
  y = NULL;
  m2 = emlrtCreateNumericArray(2, iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m2, (void *)u_data);
  emlrtSetDimensions((mxArray *)m2, u_size, 2);
  emlrtAssign(&y, m2);
  return y;
}

static real_T d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = h_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const mxArray *knot_points, const char_T
  *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(emlrtAlias(knot_points), &thisId, y);
  emlrtDestroyArray(&knot_points);
}

static int8_T emlrt_marshallIn(const mxArray *n, const char_T *identifier)
{
  int8_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(emlrtAlias(n), &thisId);
  emlrtDestroyArray(&n);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  static const int32_T iv0[2] = { 0, 0 };

  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u_data);
  emlrtSetDimensions((mxArray *)m0, u_size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  i_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static int8_T g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  int8_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "int8", false, 0U, 0);
  ret = *(int8_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv3[2];
  boolean_T bv0[2];
  int32_T i2;
  static const boolean_T bv1[2] = { false, true };

  int32_T iv4[2];
  for (i2 = 0; i2 < 2; i2++) {
    iv3[i2] = 1 + -2 * i2;
    bv0[i2] = bv1[i2];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 2U,
    iv3, bv0, iv4);
  ret->size[0] = iv4[0];
  ret->size[1] = iv4[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

void triangular_algorithm_TBM_api(const mxArray * const prhs[5], const mxArray
  *plhs[3])
{
  real_T (*N_0_data)[127];
  real_T (*N_1_data)[126];
  real_T (*N_2_data)[125];
  emxArray_real_T *knot_points;
  int8_T n;
  int8_T p;
  real_T u;
  int8_T m;
  int32_T N_2_size[2];
  int32_T N_1_size[2];
  int32_T N_0_size[2];
  N_0_data = (real_T (*)[127])mxMalloc(sizeof(real_T [127]));
  N_1_data = (real_T (*)[126])mxMalloc(sizeof(real_T [126]));
  N_2_data = (real_T (*)[125])mxMalloc(sizeof(real_T [125]));
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&knot_points, 2, true);

  /* Marshall function inputs */
  n = emlrt_marshallIn(emlrtAliasP(prhs[0]), "n");
  p = emlrt_marshallIn(emlrtAliasP(prhs[1]), "p");
  u = c_emlrt_marshallIn(emlrtAliasP(prhs[2]), "u");
  m = emlrt_marshallIn(emlrtAliasP(prhs[3]), "m");
  e_emlrt_marshallIn(emlrtAlias(prhs[4]), "knot_points", knot_points);

  /* Invoke the target function */
  triangular_algorithm_TBM(n, p, u, m, knot_points, *N_0_data, N_0_size,
    *N_1_data, N_1_size, *N_2_data, N_2_size);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*N_0_data, N_0_size);
  plhs[1] = b_emlrt_marshallOut(*N_1_data, N_1_size);
  plhs[2] = c_emlrt_marshallOut(*N_2_data, N_2_size);
  knot_points->canFreeData = false;
  emxFree_real_T(&knot_points);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (_coder_triangular_algorithm_TBM_api.c) */
