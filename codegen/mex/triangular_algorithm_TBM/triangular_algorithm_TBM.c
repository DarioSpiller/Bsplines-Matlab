/*
 * triangular_algorithm_TBM.c
 *
 * Code generation for function 'triangular_algorithm_TBM'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "triangular_algorithm_TBM.h"

/* Function Definitions */

/*
 * function [N_0,N_1,N_2] = triangular_algorithm_TBM(n,p,u,m,knot_points)
 */
void triangular_algorithm_TBM(int8_T n, int8_T p, real_T u, int8_T m, const
  emxArray_real_T *knot_points, real_T N_0_data[], int32_T N_0_size[2], real_T
  N_1_data[], int32_T N_1_size[2], real_T N_2_data[], int32_T N_2_size[2])
{
  int8_T k;
  int8_T max_length;
  int32_T b_n;
  int8_T L_N_1_data[126];
  int8_T yk;
  int32_T b_k;
  int32_T c_n;
  int8_T L_N_2_data[125];
  int32_T i0;
  real_T N_data[126];
  real_T num1;
  real_T num2;
  boolean_T zero_num1;
  boolean_T zero_num2;
  int8_T d;
  real_T den;
  real_T coeff;
  int8_T i1;
  int8_T i;
  real_T B;

  /*  Input: n, p, m, u, and m+1 clamped knots { u0, ..., um } */
  /*  Output: Coefficients N0,p(u), N1,p(u), ..., Nn,p(u) in N[0], N[1], ..., N[n] */
  /*  TO-BE-MEXED VERSION! */
  /*  Algorithm */
  /* % degree 0 coefficient */
  /*  rule out special cases */
  /* 'triangular_algorithm_TBM:12' if u == knot_points(1) */
  if (u == knot_points->data[0]) {
    /*  u corrisponde al primo knot point */
    /* 'triangular_algorithm_TBM:14' k = int8(1); */
    k = 1;
  } else if (u == knot_points->data[(int8_T)(m + 1) - 1]) {
    /* 'triangular_algorithm_TBM:15' elseif u == knot_points(m+1) */
    /*  u corrisponde all'ultimo knot point */
    /* 'triangular_algorithm_TBM:17' k = n+p-1; */
    k = (int8_T)((int8_T)(n + p) - 1);
  } else {
    /* 'triangular_algorithm_TBM:18' else */
    /*  now u is between the first and the last knot point */
    /*  Let u be in knot span [uk,uk+1); */
    /*  Original MATLAB code */
    /*  k = int8(sum(u > knot_points)); */
    /*  C-like code */
    /* 'triangular_algorithm_TBM:26' k = int8(0); */
    /* 'triangular_algorithm_TBM:27' while u>knot_points(k+1) */
    for (k = 0; u > knot_points->data[(int8_T)(k + 1) - 1]; k++) {
      /* 'triangular_algorithm_TBM:28' k = k + 1; */
    }
  }

  /* 'triangular_algorithm_TBM:32' max_length = n+p-1; */
  max_length = (int8_T)((int8_T)(n + p) - 1);

  /* 'triangular_algorithm_TBM:33' L_N_1 = 2:n; */
  if (n < 2) {
    b_n = 0;
  } else {
    b_n = n - 1;
  }

  if (b_n > 0) {
    L_N_1_data[0] = 2;
    yk = 2;
    for (b_k = 2; b_k <= b_n; b_k++) {
      yk++;
      L_N_1_data[b_k - 1] = yk;
    }
  }

  /* 'triangular_algorithm_TBM:34' L_N_2 = 3:n; */
  if (n < 3) {
    c_n = 0;
  } else {
    c_n = n - 2;
  }

  if (c_n > 0) {
    L_N_2_data[0] = 3;
    yk = 3;
    for (b_k = 2; b_k <= c_n; b_k++) {
      yk++;
      L_N_2_data[b_k - 1] = yk;
    }
  }

  /* 'triangular_algorithm_TBM:35' N = zeros(1,max_length); */
  b_k = max_length;
  for (i0 = 0; i0 < b_k; i0++) {
    N_data[i0] = 0.0;
  }

  /* 'triangular_algorithm_TBM:36' N_0 = zeros(1,n); */
  /* 'triangular_algorithm_TBM:37' N_1 = N_0(L_N_1); */
  N_1_size[0] = 1;
  N_1_size[1] = b_n;
  for (i0 = 0; i0 < b_n; i0++) {
    N_1_data[i0] = 0.0;
  }

  /*  derivata prima */
  /* 'triangular_algorithm_TBM:38' N_2 = N_0(L_N_2); */
  N_2_size[0] = 1;
  N_2_size[1] = c_n;
  for (i0 = 0; i0 < c_n; i0++) {
    N_2_data[i0] = 0.0;
  }

  /*  derivata seconda */
  /* 'triangular_algorithm_TBM:40' num = 0; */
  /* 'triangular_algorithm_TBM:41' den = 0; */
  /* 'triangular_algorithm_TBM:42' coeff = 0; */
  /* 'triangular_algorithm_TBM:43' coeff1 = 0; */
  /* 'triangular_algorithm_TBM:44' coeff2 = 0; */
  /* 'triangular_algorithm_TBM:46' N(k) = 1.0; */
  N_data[k - 1] = 1.0;

  /* % degree d goes from 1 to p */
  /*  globally used quantities */
  /* 'triangular_algorithm_TBM:51' KP_k = knot_points(k); */
  /* 'triangular_algorithm_TBM:52' KP_kplus1 = knot_points(k+1); */
  /* 'triangular_algorithm_TBM:53' num1 = (KP_kplus1 - u); */
  num1 = knot_points->data[k] - u;

  /* 'triangular_algorithm_TBM:54' num2 = (u - KP_k); */
  num2 = u - knot_points->data[k - 1];

  /* 'triangular_algorithm_TBM:55' zero_num1 = num1 == 0; */
  zero_num1 = (num1 == 0.0);

  /* 'triangular_algorithm_TBM:56' zero_num2 = num2 == 0; */
  zero_num2 = (num2 == 0.0);

  /* 'triangular_algorithm_TBM:58' for d = 1:p-1 */
  yk = (int8_T)(p - 1);
  for (d = 1; d <= yk; d++) {
    /* 'triangular_algorithm_TBM:60' if k - d > 0 */
    if ((int8_T)(k - d) > 0) {
      /*  right (south-west corner) term only */
      /* 'triangular_algorithm_TBM:62' den = (KP_kplus1 - knot_points(k-d+1)); */
      den = knot_points->data[k] - knot_points->data[(int8_T)((int8_T)(k - d) +
        1) - 1];

      /* 'triangular_algorithm_TBM:64' coeff = num1/den; */
      coeff = num1 / den;

      /*  to deal with singularities */
      /* 'triangular_algorithm_TBM:66' if zero_num1 && den == 0 */
      if (zero_num1 && (den == 0.0)) {
        /* 'triangular_algorithm_TBM:67' coeff = 1; */
        coeff = 1.0;
      }

      /* 'triangular_algorithm_TBM:70' N(k-d) = coeff * N((k-d)+1); */
      N_data[(int8_T)(k - d) - 1] = coeff * N_data[(int8_T)((int8_T)(k - d) + 1)
        - 1];
    }

    /* 'triangular_algorithm_TBM:75' if  d >= 2 && k>=2  && k < n + p - 1 */
    if ((d >= 2) && (k >= 2) && (k < (int8_T)((int8_T)(n + p) - 1))) {
      /*  compute internal terms */
      /* 'triangular_algorithm_TBM:77' for i = k-d+1:k-1 */
      i1 = (int8_T)(k - 1);
      for (i = (int8_T)((k - d) + 1); i <= i1; i++) {
        /* 'triangular_algorithm_TBM:78' if i > 0 && i + d < m + 1 */
        if ((i > 0) && ((int8_T)(i + d) < (int8_T)(m + 1))) {
          /* 'triangular_algorithm_TBM:79' coeff1 = (u - knot_points(i))/(knot_points(i+d) - knot_points(i)); */
          den = knot_points->data[(int8_T)(i + d) - 1] - knot_points->data[i - 1];

          /* 'triangular_algorithm_TBM:80' coeff2 = (knot_points(i+d+1) - u)/(knot_points(i+d+1) - knot_points(i+1)); */
          coeff = knot_points->data[(int8_T)((int8_T)(i + d) + 1) - 1] - u;
          B = knot_points->data[(int8_T)((int8_T)(i + d) + 1) - 1] -
            knot_points->data[i];

          /* 'triangular_algorithm_TBM:81' N(i) = coeff1 * N(i) + coeff2 * N(i+1); */
          N_data[i - 1] = (u - knot_points->data[i - 1]) / den * N_data[i - 1] +
            coeff / B * N_data[i];
        }
      }
    }

    /* 'triangular_algorithm_TBM:86' if k <= n + (p-1) - d */
    if (k <= (int8_T)((int8_T)(n + (int8_T)(p - 1)) - d)) {
      /*  let (north-west corner) term only */
      /* 'triangular_algorithm_TBM:88' den = (knot_points(k+d) - KP_k); */
      den = knot_points->data[(int8_T)(k + d) - 1] - knot_points->data[k - 1];

      /* 'triangular_algorithm_TBM:90' coeff = num2/den; */
      coeff = num2 / den;

      /*  to deal with singularities */
      /* 'triangular_algorithm_TBM:92' if zero_num2 && den == 0 */
      if (zero_num2 && (den == 0.0)) {
        /* 'triangular_algorithm_TBM:93' coeff = 1; */
        coeff = 1.0;
      }

      /* 'triangular_algorithm_TBM:96' N(k) = coeff * N(k); */
      N_data[k - 1] *= coeff;
    }

    /* 'triangular_algorithm_TBM:100' if d == p-2 */
    if (d == (int8_T)(p - 2)) {
      /*  First derivative */
      /* 'triangular_algorithm_TBM:103' if N(1) ~= 1 && N(end) ~= 1 */
      if ((N_data[0] != 1.0) && (N_data[max_length - 1] != 1.0)) {
        /*  Common case */
        /* 'triangular_algorithm_TBM:105' N_1 = N(L_N_1); */
        N_1_size[0] = 1;
        N_1_size[1] = b_n;
        for (i0 = 0; i0 < b_n; i0++) {
          N_1_data[i0] = N_data[L_N_1_data[i0] - 1];
        }
      } else if (N_data[0] == 1.0) {
        /* 'triangular_algorithm_TBM:106' elseif N(1) == 1 */
        /*  to deal with singularities */
        /* 'triangular_algorithm_TBM:108' N_1 = N(L_N_1 - 1); */
        N_1_size[0] = 1;
        N_1_size[1] = b_n;
        for (i0 = 0; i0 < b_n; i0++) {
          N_1_data[i0] = N_data[(int8_T)(L_N_1_data[i0] - 1) - 1];
        }
      } else {
        if (N_data[max_length - 1] == 1.0) {
          /* 'triangular_algorithm_TBM:109' elseif N(end) == 1 */
          /*  to deal with singularities */
          /* 'triangular_algorithm_TBM:111' N_1(end) = 1; */
          N_1_data[N_1_size[1] - 1] = 1.0;
        }
      }
    } else {
      if (d == (int8_T)(p - 3)) {
        /* 'triangular_algorithm_TBM:114' elseif d == p-3 */
        /*  Second derivative */
        /* 'triangular_algorithm_TBM:117' if N(1) ~= 1 && N(end) ~= 1 */
        if ((N_data[0] != 1.0) && (N_data[max_length - 1] != 1.0)) {
          /*  Common case */
          /* 'triangular_algorithm_TBM:119' N_2 = N(L_N_2); */
          N_2_size[0] = 1;
          N_2_size[1] = c_n;
          for (i0 = 0; i0 < c_n; i0++) {
            N_2_data[i0] = N_data[L_N_2_data[i0] - 1];
          }
        } else if (N_data[0] == 1.0) {
          /* 'triangular_algorithm_TBM:120' elseif N(1) == 1 */
          /*  to deal with singularities */
          /* 'triangular_algorithm_TBM:122' N_2 = N(L_N_2 - 2); */
          N_2_size[0] = 1;
          N_2_size[1] = c_n;
          for (i0 = 0; i0 < c_n; i0++) {
            N_2_data[i0] = N_data[(int8_T)(L_N_2_data[i0] - 2) - 1];
          }
        } else {
          if (N_data[max_length - 1] == 1.0) {
            /* 'triangular_algorithm_TBM:123' elseif N(end) == 1 */
            /*  to deal with singularities */
            /* 'triangular_algorithm_TBM:125' N_2(end) = 1; */
            N_2_data[N_2_size[1] - 1] = 1.0;
          }
        }
      }
    }
  }

  /*  array N[0..n] has the coefficients. */
  /* 'triangular_algorithm_TBM:134' N_0 = N(1:n); */
  if (1 > n) {
    b_k = 0;
  } else {
    b_k = n;
  }

  N_0_size[0] = 1;
  N_0_size[1] = b_k;
  for (i0 = 0; i0 < b_k; i0++) {
    N_0_data[i0] = N_data[i0];
  }
}

/* End of code generation (triangular_algorithm_TBM.c) */
