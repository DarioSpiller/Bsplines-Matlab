/*
 * triangular_algorithm_TBM_types.h
 *
 * Code generation for function 'triangular_algorithm_TBM'
 *
 */

#ifndef __TRIANGULAR_ALGORITHM_TBM_TYPES_H__
#define __TRIANGULAR_ALGORITHM_TBM_TYPES_H__

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_int8_T_1x125
#define struct_emxArray_int8_T_1x125
struct emxArray_int8_T_1x125
{
    int8_T data[125];
    int32_T size[2];
};
#endif /*struct_emxArray_int8_T_1x125*/
#ifndef typedef_emxArray_int8_T_1x125
#define typedef_emxArray_int8_T_1x125
typedef struct emxArray_int8_T_1x125 emxArray_int8_T_1x125;
#endif /*typedef_emxArray_int8_T_1x125*/
#ifndef struct_emxArray_int8_T_1x126
#define struct_emxArray_int8_T_1x126
struct emxArray_int8_T_1x126
{
    int8_T data[126];
    int32_T size[2];
};
#endif /*struct_emxArray_int8_T_1x126*/
#ifndef typedef_emxArray_int8_T_1x126
#define typedef_emxArray_int8_T_1x126
typedef struct emxArray_int8_T_1x126 emxArray_int8_T_1x126;
#endif /*typedef_emxArray_int8_T_1x126*/
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
#endif /*struct_emxArray_real_T*/
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /*typedef_emxArray_real_T*/

#endif
/* End of code generation (triangular_algorithm_TBM_types.h) */
