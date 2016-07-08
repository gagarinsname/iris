#include <stdio.h>
#include <math.h>
#include "stddefs.h"

#define a_ref(a_1,a_2) a[(a_2)*m + (a_1)]
#define u_ref(a_1,a_2) u[(a_2)*m + (a_1)]
#define v_ref(a_1,a_2) v[(a_2)*n + (a_1)]
#define d_sign(a,b) (((b)>=0)?(((a)>=0)?(a):(-(a))):(((a)>= 0)?(-(a)):(a)))
#define max_(a,b) (((a)>(b))?(a):(b))

int BPL_Matrix_SVD(
  const double *a,  // вещественный двумерный массив размера m на n, в первых m стpоках которого задается исходная матрица
  double *w,  // вещественный вектоp длины n значений w(k) = dk вычисленных сингулярных чисел матрицы a
  double *u,  // вещественный двумерный массив размера m на n, в котоpом при matu = .true. содержатся вычисленные первые n столбцов матрицы u; при matu = .false. массив u используется как рабочий
  int matu,   // признак необходимости вычисления матрицы u
  double *v,  // вещественный двумерный массив размера m на n, в котоpом при matv = .true. содержится вычисленная матрица v; при matv = .false. массив v не используется
  int matv,   // признак необходимости вычисления матрицы v
  int m,      // число стpок матриц a и u
  int n,      // число столбцов матриц a и u и порядок матрицы v
  double *rv1)// вещественный рабочий вектоp длины n
{
  double d__1, d__2, d__3, d__4,ama999, c__, f, g, h__;
  double s, x, y, z__, scale,anorm;
  int i__,j,k,l = 0,i1,k1,l1 = 0,ii,kk,ll,mn,its;

    // Parameter adjustments
    --rv1;
    --w;
    v -= (n + 1);
    u -= (m + 1);
    a -= (m + 1);
    // Function Body
    for (i__ = 1; i__ <= m; ++i__)
    {
      for (j = 1; j <= n; ++j)
        u_ref(i__, j) = a_ref(i__, j);
    }
    i__ = 0;
    g = 0.;
    scale = 0.;
    anorm = 0.;
    for (i__ = 1; i__ <= n; ++i__)
    {
      l = i__ + 1;
      rv1[i__] = scale * g;
      g = 0.;
      s = 0.;
      scale = 0.;
      if (i__ <= m)
      {
        for (k = i__; k <= m; ++k)
        {
          d__1 = u_ref(k, i__);
          scale += MIA_abs(d__1);
        }
        if (scale != 0.)
        {
          for (k = i__; k <= m; ++k)
          {
            u_ref(k, i__) /= scale;
            // Computing 2nd power
            d__1 = u_ref(k, i__);
            s += d__1 * d__1;
          }
          f = u_ref(i__, i__);
          d__1 = sqrt(s);
          g = -d_sign(d__1,f);
          h__ = f * g - s;
          u_ref(i__, i__) = f - g;
          if (i__ != n)
          {
            for (j = l; j <= n; ++j)
            {
              s = 0.;
              for (k = i__; k <= m; ++k)
                s += u_ref(k, i__) * u_ref(k, j);
              f = s / h__;
              for (k = i__; k <= m; ++k)
                u_ref(k, j) = u_ref(k, j) + f * u_ref(k, i__);
            }
          }
          for (k = i__; k <= m; ++k)
            u_ref(k, i__) *= scale;
        }
      }
      w[i__] = scale * g;
      g = 0.;
      s = 0.;
      scale = 0.;
      if (!(i__ > m || i__ == n))
      {
        for (k = l; k <= n; ++k)
        {
          d__1 = u_ref(i__, k);
          scale += MIA_abs(d__1);
        }
        if (scale != 0.)
        {
          for (k = l; k <= n; ++k)
          {
            u_ref(i__, k) /= scale;
            // Computing 2nd power
            d__1 = u_ref(i__, k);
            s += d__1 * d__1;
          }
          f = u_ref(i__, l);
          d__1 = sqrt(s);
          g = -d_sign(d__1,f);
          h__ = f * g - s;
          u_ref(i__, l) = f - g;
          for (k = l; k <= n; ++k)
            rv1[k] = u_ref(i__, k) / h__;
          if (i__ != m)
          {
            for (j = l; j <= m; ++j)
            {
              s = 0.;
              if (l <= n)
              {
                for (k = l; k <= n; ++k)
                  s += u_ref(j, k) * u_ref(i__, k);
                for (k = l; k <= n; ++k)
                  u_ref(j, k) = u_ref(j, k) + s * rv1[k];
              }
            }
          }
          for (k = l; k <= n; ++k)
            u_ref(i__, k) = scale * u_ref(i__, k);
        }
      }
      // Computing MAX
      d__1 = w[i__];
      d__2 = rv1[i__];
      d__3 = anorm;
      d__4 = MIA_abs(d__1)+MIA_abs(d__2);
      anorm = max_(d__3,d__4);
    }
    if (matv)
    {
      for (ii = 1; ii <= n; ++ii)
      {
        i__ = n + 1 - ii;
        if (i__ != n)
        {
          if (g != 0.)
          {
            for (j = l; j <= n; ++j)
              v_ref(j, i__) = u_ref(i__, j) / u_ref(i__, l) / g;
            for (j = l; j <= n; ++j)
            {
              s = 0.;
              for (k = l; k <= n; ++k)
                s += u_ref(i__, k) * v_ref(k, j);
              for (k = l; k <= n; ++k)
                v_ref(k, j) = v_ref(k, j) + s * v_ref(k, i__);
            }
          }
          for (j = l; j <= n; ++j)
          {
            v_ref(i__, j) = 0.;
            v_ref(j, i__) = 0.;
          }
        }
        v_ref(i__, i__) = 1.;
        g = rv1[i__];
        l = i__;
      }
    }
    if (matu)
    {
      mn = n;
      if (m < n)
        mn = m;
      for (ii = 1; ii <= mn; ++ii)
      {
        i__ = mn + 1 - ii;
        l = i__ + 1;
        g = w[i__];
        if (i__ < n)
        {
          for (j = l; j <= n; ++j)
            u_ref(i__, j) = 0.;
        }
        if (g == 0.)
        {
          for (j = i__; j <= m; ++j)
            u_ref(j, i__) = 0.;
        }
        else
        {
          if (i__ != mn)
          {
            for (j = l; j <= n; ++j)
            {
              s = 0.;
              for (k = l; k <= m; ++k)
                s += u_ref(k, i__) * u_ref(k, j);
              f = s / u_ref(i__, i__) / g;
              for (k = i__; k <= m; ++k)
                u_ref(k, j) = u_ref(k, j) + f * u_ref(k, i__);
            }
          }
          for (j = i__; j <= m; ++j)
            u_ref(j, i__) /= g;
        }
        u_ref(i__, i__) += 1.;
      }
    }
    for (kk = 1; kk <= n; ++kk)
    {
      k1 = n - kk;
      k = k1 + 1;
      its = 0;
L37:
      for (ll = 1; ll <= k; ++ll)
      {
        l1 = k - ll;
        l = l1 + 1;
        d__1 = rv1[l];
        if (MIA_abs(d__1) + anorm == anorm)
          goto L42;
        d__1 = w[l1];
        if (MIA_abs(d__1) + anorm == anorm)
          break;
      }
      c__ = 0.;
      s = 1.;
      for (i__ = l; i__ <= k; ++i__)
      {
        f = s * rv1[i__];
        rv1[i__] = c__ * rv1[i__];
        if (MIA_abs(f) + anorm == anorm)
          break;
        g = w[i__];
        h__ = sqrt(f * f + g * g);
        w[i__] = h__;
        c__ = g / h__;
        s = -f / h__;
        if (matu)
        {
          for (j = 1; j <= m; ++j)
          {
            y = u_ref(j, l1);
            z__ = u_ref(j, i__);
            u_ref(j, l1) = y * c__ + z__ * s;
            u_ref(j, i__) = -y * s + z__ * c__;
          }
        }
      }
L42:
      z__ = w[k];
      if (l == k)
        goto L48;
      if (its == 30)
        //" БИБЛИOTEKA HИBЦ MГУ, ПOДПPOГPAMMA ACP1R(ACP1D): ФATAЛЬHAЯ OШИБKA\n\n"
        //" N  1 - PAЗЛOЖEHИE ПPEKPAЩEHO, T.K. ДЛЯ BЫЧИCЛEHИЯ\n"
        //"        CИHГУЛЯPHOГO ЗHAЧEHИЯ C ИHДEKCOM IERR TPEБУETCЯ\n"
        //"        БOЛEE 30 ИTEPAЦИЙ. ПPABИЛЬHO BЫЧИCЛEHЫ CИHГУЛЯPHЫE\n"
        //"        ЗHAЧEHИЯ C ИHДEKCAMИ IERR+1,IERR+2,... \n");
        return k;
      ++its;
      x = w[l];
      y = w[k1];
      g = rv1[k1];
      h__ = rv1[k];
      f = ((y - z__) * (y + z__) + (g - h__) * (g + h__)) / (h__ * 2. * y);
      g = sqrt(f * f + 1.);
      f = ((x-z__)*(x+z__)+h__ * (y / (f + d_sign(g,f)) - h__))/x;
      c__ = 1.;
      s = 1.;
      for (i1 = l; i1 <= k1; ++i1)
      {
        i__ = i1 + 1;
        g = rv1[i__];
        y = w[i__];
        h__ = s * g;
        g = c__ * g;
        // Computing MAX
        d__1 = MIA_abs(f);
        d__2 = MIA_abs(h__);
        ama999 = max_(d__1,d__2);
        // Computing 2nd power
        d__1 = f / ama999;
        // Computing 2nd power */
        d__2 = h__ / ama999;
        z__ = ama999 * sqrt(d__1 * d__1 + d__2 * d__2);
        rv1[i1] = z__;
        c__ = f / z__;
        s = h__ / z__;
        f = x * c__ + g * s;
        g = -x * s + g * c__;
        h__ = y * s;
        y *= c__;
        if (matv)
        {
          for (j = 1; j <= n; ++j)
          {
            x = v_ref(j, i1);
            z__ = v_ref(j, i__);
            v_ref(j, i1) = x * c__ + z__ * s;
            v_ref(j, i__) = -x * s + z__ * c__;
          }
        }
        z__ = sqrt(f * f + h__ * h__);
        w[i1] = z__;
        if (z__ != 0.)
        {
          c__ = f / z__;
          s = h__ / z__;
        }
        f = c__ * g + s * y;
        x = -s * g + c__ * y;
        if (matu)
        {
          for (j = 1; j <= m; ++j)
          {
            y = u_ref(j, i1);
            z__ = u_ref(j, i__);
            u_ref(j, i1) = y * c__ + z__ * s;
            u_ref(j, i__) = -y * s + z__ * c__;
          }
        }
      }
      rv1[l] = 0.;
      rv1[k] = f;
      w[k] = x;
      goto L37;
L48:
      if (z__ < 0.)
      {
        w[k] = -z__;
        if (matv)
        {
          for (j = 1; j <= n; ++j)
            v_ref(j, k) = -v_ref(j, k);
        }
      }
    }
    return 0;
}

// im = U*L*Vt
void BPL_Matrix_RestoreImageFromSVD(
  unsigned char* im,
  const double* U,
  const double* V,
  const double* L,
  int xs,
  int ys,
  int n)
{
  int i,j,k;
  double s;

    for (i=0;i<ys;i++)
      for (j=0;j<xs;j++)
      {
        s = 0.;
        for (k=n-1;k>=0;k--)
          s += L[k]*V[k*ys+i]*U[k*xs+j];
        k = (int)(s+.5);
        if (k>255)
          k = 255;
        else
          if (k<0)
            k = 0;
        im[i*xs+j] = (unsigned char)k;
      }
}
