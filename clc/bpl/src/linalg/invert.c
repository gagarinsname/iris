#include <math.h>
#include "stddefs.h"

staticFunc int ownBPL_afh1d_c(
  double *a,
  int n)
{
  double x,rn;
  int i__,j,k,ip,iq,ir;
  int ip1=0;

    // Parameter adjustments
    --a;
    // Function Body
    rn = 1. / ((double) n * 16.);
    ip = 1;
    for (i__ = 1; i__ <= n; ++i__)
    {
      iq = ip;
      ir = 1;
      for (j = 1; j <= i__; ++j)
      {
        x = a[ip];
        if (j != 1)
        {
          for (k = iq; k <= ip1; ++k)
          {
            x -= a[k] * a[ir];
            ++ir;
          }
        }
        if (i__ != j)
          a[ip] = x / a[ir];
        else
        {
          if (a[ip] + x * rn <= a[ip])
            //ÇAÄAHHAß MATPÈÖA HE ßBËßETCß ÏOËOÆÈTEËÜHO OÏPEÄEËEHHOÉ
            //ÏOÝTOMÓ TPEÓÃOËÜHOE PAÇËOÆEHÈE ÝTOÉ MATPÈÖÛ HEBOÇMOÆHO
            return 65;
          a[ip] = sqrt(x);
        }
        ip1 = ip;
        ++ip;
        ++ir;
      }
    }
    return 0;
}

staticFunc void ownBPL_ash1dt_c(
  double *a,
  double *b,
  double *x,
  int n)
{
  int i__, k;
  double t;
  int n1, ii, kk, ip, iq, is, iw, im1;

    // Parameter adjustments
    --x;
    --b;
    --a;
    // Function Body
    ip = 1;
    iw = 0;
    for (i__ = 1; i__ <= n; ++i__)
    {
      t = b[i__];
      im1 = i__ - 1;
      if (iw == 0)
      {
        if (t != 0.)
          iw = i__;
        ip += im1;
      }
      else
      {
        ip = ip + iw - 1;
        for (k = iw; k <= im1; ++k)
        {
          t -= a[ip] * x[k];
          ++ip;
        }
      }
      x[i__] = t / a[ip];
      ++ip;
    }
    n1 = n + 1;
    for (i__ = 1; i__ <= n; ++i__)
    {
      ii = n1 - i__;
      --ip;
      is = ip;
      iq = ii + 1;
      t = x[ii];
      if (n >= iq)
      {
        kk = n;
        for (k = iq; k <= n; ++k)
        {
          t -= a[is] * x[kk];
          --kk;
          is -= kk;
        }
      }
      x[ii] = t / a[is];
    }
}

// invert packed matrix
int BPL_Matrix_InvertPacked(
  double *ainv,   // OUT: inverted matrix in packed form
  double *a,      // IN:  source matrix in packed form (destroyed)
  int n)          // IN:  size of matrix
{
  int i__,j,l,i1,n1,ln,retval;

    // Parameter adjustments
    --ainv;
    --a;
    // Function Body
    l = 1;
    n1 = n - 1;
    ln = n;
    if ((retval = ownBPL_afh1d_c(&a[1], n))!=0)
      return retval;
    for (i__ = 1; i__ <= n; ++i__)
    {
      for (j = l; j <= ln; ++j)
        ainv[j] = 0.;
      i1 = l + i__ - 1;
      ainv[i1] = 1.;
      ownBPL_ash1dt_c(&a[1], &ainv[l], &ainv[l], n);
      l += i__;
      ln = l + n1;
    }
    return 0;
}

// pack symmetric matrix
void BPL_Matrix_Pack(
  double *b,
  const double *a,
  int n)
{
  int i__,j,k;

    // Parameter adjustments
    a -= n + 1;
    --b;
    // Function Body
    k = 1;
    for (i__ = 1; i__ <= n; ++i__)
      for (j = 1; j <= i__; ++j)
        b[k++] = a[i__*n+j];
}

// unpack to symmetric matrix
void BPL_Matrix_Unpack(
  double *b,
  const double *a,
  int n)
{
  int i__, j, k, i1;

    // Parameter adjustments
    --a;
    b -= n + 1;
    // Function Body
    i1 = n * (n + 1) / 2;
    j = n;
    do
    {
      for (k = 1; k <= j; ++k)
      {
        i__ = j + 1 - k;
        b[j*n + i__] = a[i1];
        --i1;
      }
      --j;
    } while (j >= 1);
    if (n < 2)
      return;
    for (i__ = 2; i__ <= n; ++i__)
    {
      i1 = i__ - 1;
      for (j = 1; j <= i1; ++j)
        b[j*n + i__] = b[i__*n + j];
    }
}

// unpack to symmetric matrix
void BPL_Matrix_UnpackFloat(
  float *b,
  const float *a,
  int n)
{
  int i__, j, k, i1;

    // Parameter adjustments
    --a;
    b -= n + 1;
    // Function Body
    i1 = n * (n + 1) / 2;
    j = n;
    do
    {
      for (k = 1; k <= j; ++k)
      {
        i__ = j + 1 - k;
        b[j*n + i__] = a[i1];
        --i1;
      }
      --j;
    } while (j >= 1);
    if (n < 2)
      return;
    for (i__ = 2; i__ <= n; ++i__)
    {
      i1 = i__ - 1;
      for (j = 1; j <= i1; ++j)
        b[j*n + i__] = b[i__*n + j];
    }
}
