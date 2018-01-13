/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     13 March 2006                                       */
/*      Revision:    1.0.00                                              */
/*      Purpose:     Matrix solution operations                          */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include "bpl.h"

// Solve symmetric real linear system Ax=b
// Source: BChA NIVC MGU, refer to as: ash0d_c
int BPL_Matrix_SymmetricLinearSolve(
  double* a,    // IN:  matrix A (size N*N)
  double* b,    // IN:  vector b
  double* x,    // OUT: vector x
  double* s,    // TMP: buffer for N values
  int n,        // IN:  order of system
  int isFirst)  // IN:  is function called first time with this matrix
{
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
  // System generated locals
  int a_dim1, a_offset, i__1, i__2;
  // Local variables
  double f;
  int i__, j, k;
  double z__;
  int i1, j1, n1;

    // Parameter adjustments
    --s;
    --x;
    --b;
    a_dim1 = n;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    // Function Body
    if (isFirst)
    {
      s[1] = 1. / a_ref(1, 1);
      j1 = 1;
      j = 2;
      goto L4;
L1:
      j1 = j;
      ++j;
      i__1 = j1;
      for (i__ = 2; i__ <= i__1; ++i__)
      {
        i1 = i__ - 1;
        z__ = a_ref(i__, j);
        i__2 = i1;
        for (k = 1; k <= i__2; ++k)
          z__ -= a_ref(k, j) * a_ref(k, i__);
        a_ref(i__, j) = z__;
      }
L4:
      z__ = a_ref(j, j);
      i__1 = j1;
      for (k = 1; k <= i__1; ++k)
      {
        f = a_ref(k, j);
        z__ -= f * s[k] * f;
      }
      s[j] = 1. / z__;
      i__1 = j1;
      for (i__ = 1; i__ <= i__1; ++i__)
      {
        a_ref(i__, j) = a_ref(i__, j) * s[i__];
      }
      if (j != n)
        goto L1;
    }
    x[1] = b[1];
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__)
    {
      z__ = b[i__];
      i1 = i__ - 1;
      i__2 = i1;
      for (k = 1; k <= i__2; ++k)
        z__ -= a_ref(k, i__) * x[k];
      x[i__] = z__;
    }
    n1 = n + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      i__ = n1 - j;
      z__ = x[i__] * s[i__];
      i1 = i__ + 1;
      if (i1 <= n)
      {
        i__2 = n;
        for (k = i1; k <= i__2; ++k)
          z__ -= a_ref(i__, k) * x[k];
      }
      x[i__] = z__;
    }
    return 0;
#undef a_ref
}

// float version of SymmetricLinearSolve
int BPL_Matrix_SymmetricLinearSolveFloat(
  float* a,     // IN:  matrix A (size N*N)
  float* b,     // IN:  vector b
  float* x,     // OUT: vector x
  float* s,     // TMP: buffer for N values
  int n,        // IN:  order of system
  int isFirst)  // IN:  is function called first time with this matrix
{
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
  // System generated locals
  int a_dim1, a_offset, i__1, i__2;
  // Local variables
  float f;
  int i__, j, k;
  float z__;
  int i1, j1, n1;

    // Parameter adjustments
    --s;
    --x;
    --b;
    a_dim1 = n;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    // Function Body
    if (isFirst)
    {
      s[1] = 1.f/a_ref(1, 1);
      j1 = 1;
      j = 2;
      goto L4;
L1:
      j1 = j;
      ++j;
      i__1 = j1;
      for (i__ = 2; i__ <= i__1; ++i__)
      {
        i1 = i__ - 1;
        z__ = a_ref(i__, j);
        i__2 = i1;
        for (k = 1; k <= i__2; ++k)
          z__ -= a_ref(k, j) * a_ref(k, i__);
        a_ref(i__, j) = z__;
      }
L4:
      z__ = a_ref(j, j);
      i__1 = j1;
      for (k = 1; k <= i__1; ++k)
      {
        f = a_ref(k, j);
        z__ -= f * s[k] * f;
      }
      s[j] = 1.f/z__;
      i__1 = j1;
      for (i__ = 1; i__ <= i__1; ++i__)
      {
        a_ref(i__, j) = a_ref(i__, j) * s[i__];
      }
      if (j != n)
        goto L1;
    }
    x[1] = b[1];
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__)
    {
      z__ = b[i__];
      i1 = i__ - 1;
      i__2 = i1;
      for (k = 1; k <= i__2; ++k)
        z__ -= a_ref(k, i__) * x[k];
      x[i__] = z__;
    }
    n1 = n + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      i__ = n1 - j;
      z__ = x[i__] * s[i__];
      i1 = i__ + 1;
      if (i1 <= n)
      {
        i__2 = n;
        for (k = i1; k <= i__2; ++k)
          z__ -= a_ref(i__, k) * x[k];
      }
      x[i__] = z__;
    }
    return 0;
#undef a_ref
}
