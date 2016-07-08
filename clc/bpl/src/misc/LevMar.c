/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     2012                                                */
/*      Modified:    20 April 2012 - "BPL-ized" from BCha NIVC MGU format*/
/*      Revision:    2.0.00                                              */
/*      Purpose:     Levengerg-Marquardt algorithm                       */
/*      Authors:                                                         */
/*          NIVC MGU (http://www.srcc.msu.su)                            */
/*          http://num-anal.srcc.msu.ru/lib_na/cat/mn_htm_c/mnavr_c.htm  */
/*-----------------------------------------------------------------------*/
#include <math.h>
#include "bpl.h"

#define myabs(v) (((v)>0)?(v):(-(v)))
#define mymax(v,u) (((v)>(u))?(v):(u))
#define mymin(v,u) (((v)<(u))?(v):(u))

static int ownBPL_afh1r_c(
  double *a,
  int n)
{
  // System generated locals
  int i__1, i__2, i__3;
  // Local variables
  static int i__, j, k;
  static double x;
  static int ip, iq, ir;
  static double rn;
  static int ip1;

    // Parameter adjustments
    --a;
    // Function Body
    ip1 = 1;
    rn = 1. / (n * 16.);
    ip = 1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      iq = ip;
      ir = 1;
      i__2 = i__;
      for (j = 1; j <= i__2; ++j)
      {
        x = a[ip];
        if (j != 1)
        {
          i__3 = ip1;
          for (k = iq; k <= i__3; ++k)
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
            return 65;
          a[ip] = (double)sqrt(x);
        }
        ip1 = ip;
        ++ip;
        ++ir;
      }
    }
    return 0;
}

static void ownBPL_ash1rt_c(
  double *a,
  double *b,
  double *x,
  int n)
{
  // System generated locals
  int i__1, i__2;
  // Local variables
  static int i__, k;
  static double t;
  static int n1, ii, kk, ip, iq, is, iw, im1;

    // Parameter adjustments
    --x;
    --b;
    --a;
    // Function Body
    ip = 1;
    iw = 0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__)
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
        i__2 = im1;
        for (k = iw; k <= i__2; ++k)
        {
          t -= a[ip] * x[k];
          ++ip;
        }
      }
      x[i__] = t / a[ip];
      ++ip;
    }
    n1 = n + 1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      ii = n1 - i__;
      --ip;
      is = ip;
      iq = ii + 1;
      t = x[ii];
      if (n >= iq)
      {
        kk = n;
        i__2 = n;
        for (k = iq; k <= i__2; ++k)
        {
          t -= a[is] * x[kk];
          --kk;
          is -= kk;
        }
      }
      x[ii] = t / a[is];
    }
}

static int ownBPL_ash1r_c(
  double *a,
  double *b,
  double *x,
  int n)
{
  int ierr=0;

    // Parameter adjustments
    --a;
    --x;
    --b;
    // Function Body
    ierr = ownBPL_afh1r_c(&a[1], n);
    if (ierr != 0)
      return ierr;
    ownBPL_ash1rt_c(&a[1], &b[1], &x[1], n);
    return ierr;
}

// un-constraint minimization of differentiable function of multiple variables
// represented as sum of squares by Levenberg-Marquardt (LM) procedure
//  min F(x) , F(x)=f_1^2+...+f_m^2
int BPL_OPT_LevenbergMarquardt(
  BPL_OPT_LevenbergMarquardt_Callback func,      // function callback
  int m,          // number of sub-functions
  int n,          // dimension of variables space
  int nsig,       // 10^-nsig is a required precision of argument
  double eps,     // required precision of function
  double delta,   // required precision of gradient
  int maxfn,      // maximum allowed callbacks to 'func'
  int iopt,       // method of optimization: 0|1|2 = Brown scheme|LM scheme with 
                  //   LM parameter auto selection|LM scheme regulated by 'parm'
  double *parm,   // [0]-initial param value, [1]-scalar multiplier to change
                  //   [2]-upper limit for param, [3]-criterion
  double *x,      // IN/OUT: initial point/end point
  double *ssq,    // OUT: functional value at the end
  double *f,      // OUT: m values of f_i at the end
  double *xjac,   // m*n matrix - Jacobian of f
  int ixjac,      // number of strings in 'xjac'
  double *xjtj,   // service, size=(n+1)*n/2
  double *work,   // service, size = 5n + 2m + (n + 1) * n/2
  int *infer)     // OUT: criterion of stop: 0-no convergence, 1-argument precision reached
                  //  2-function precision reached, 3-gradient precision reached
{
  // Initialized data
  const double sig = 6.3;
  // System generated locals
  int i__1, i__2, i__3;
  double r__1, r__2, r__3, r__4, r__5, r__6;
  double d__1, d__2;
  // Local variables
  int ijac=0, ifml;
  double prec;
  int ifpl;
  double xdif;
  int ifpu, iter;
  double fosq=0.;
  int ifml1, ifpl1;
  double cons2, erl2x=0., g;
  int i__, j, k, l, ieval;
  double xdabs, sqdif, xhold, dnorm;
  int izero, ibad1j, ical1i, ical1j, igrad1, ixbad1, irad1i;
  double delta2;
  int iscal1, ifml1i, ifpl1i, idelx1, ielx1i, ielx1j, inew1j, ixnew1;
  double al, fosqs4=0., hh;
  int ki, li, im;
  double fo=0.;
  int ir, is, js, igradl;
  double up=0.;
  int iscall, igradu, iscalu, idelxl;
  double relcon, onesfo=0.;
  int idelxu, icount=0;
  double ssqold=0.;
  int ixnewl;
  double rel;
  int iri;
  double rhh, dsq;
  int isw;
  double sum, erl2;
  int ibad;
  // mia vars
  int ier=0;

    // Parameter adjustments
    --f;
    --x;
    --parm;
    --xjac;
    --xjtj;
    --work;
    // Function Body
    if (m <= 0 || n <= 0 || iopt < 0 || iopt > 2)
      return 130;
    if ((iopt == 2)&&(parm[2] <= 1. || parm[1] <= 0.))
      return 130;
    d__1 = (double) 10.;
    d__2 = (double) (-sig - 1.);
    prec = pow(d__1,d__2);
    d__1 = (double) 10.;
    d__2 = (double) (-sig * .5);
    rel = pow(d__1,d__2);
    i__1 = -nsig;
    relcon = pow(10., (double)i__1);
    igrad1 = (n + 1) * n / 2;
    igradl = igrad1 + 1;
    igradu = igrad1 + n;
    idelx1 = igradu;
    idelxl = idelx1 + 1;
    idelxu = idelx1 + n;
    iscal1 = idelxu;
    iscall = iscal1 + 1;
    iscalu = iscal1 + n;
    ixnew1 = iscalu;
    ixnewl = ixnew1 + 1;
    ixbad1 = ixnew1 + n;
    ifpl1 = ixbad1 + n;
    ifpl = ifpl1 + 1;
    ifpu = ifpl1 + m;
    ifml1 = ifpu;
    ifml = ifml1 + 1;
    al = 1.;
    cons2 = .1;
    if (iopt == 0)
      goto L40;
    if (iopt == 1)
      goto L20;
    al = parm[1];
    fo = parm[2];
    up = parm[3];
    cons2 = parm[4] * .5;
    goto L30;
L20:
    al = .01;
    fo = 2.;
    up = 120.;
L30:
    onesfo = 1. / fo;
    fosq = fo * fo;
/* Computing 4th power */
    r__1 = fosq, r__1 *= r__1;
    fosqs4 = r__1 * r__1;
L40:
    ieval = 0;
    delta2 = delta * .5;
    erl2 = 1e10;
    ibad = -99;
    isw = 1;
    iter = -1;
    *infer = 0;
    ier = 0;
    i__1 = idelxu;
    for (j = idelxl; j <= i__1; ++j)
      work[j] = 0.;
    goto L340;
L60:
    ssqold = *ssq;
    if (*infer > 0 || ijac >= n || iopt == 0 || icount > 0)
      goto L110;
    ++ijac;
    dsq = 0.;
    i__1 = idelxu;
    for (j = idelxl; j <= i__1; ++j)
      dsq += work[j] * work[j];
    if (dsq <= 0.)
      goto L110;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      ifml1i = ifml1 + i__;
      g = f[i__] - work[ifml1i];
      k = i__;
      i__2 = idelxu;
      for (j = idelxl; j <= i__2; ++j)
      {
        g += xjac[k] * work[j];
        k += m;
      }
      g /= dsq;
      k = i__;
      i__2 = idelxu;
      for (j = idelxl; j <= i__2; ++j)
      {
        xjac[k] -= g * work[j];
        k += m;
      }
    }
    goto L160;
L110:
    ijac = 0;
    k = 0;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      xdabs = (r__1 = x[j], myabs(r__1));
      hh = rel * mymax(xdabs,.1);
      xhold = x[j];
      x[j] += hh;
      (*func)(&x[1], m, n, &work[ifpl]);
      ++ieval;
      x[j] = xhold;
      if (isw == 1)
      {
        rhh = 1. / hh;
        i__2 = m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            ++k;
            ifpl1i = ifpl1 + i__;
            xjac[k] = (work[ifpl1i] - f[i__]) * rhh;
        }
      }
      else
      {
        x[j] = xhold - hh;
        (*func)(&x[1], m, n, &work[ifml]);
        ++ieval;
        x[j] = xhold;
        rhh = .5 / hh;
        i__2 = ifpu;
        for (i__ = ifpl; i__ <= i__2; ++i__)
        {
          ++k;
          im = i__ + m;
          xjac[k] = (work[i__] - work[im]) * rhh;
        }
      }
    }
L160:
    erl2x = erl2;
    erl2 = 0.;
    k = 0;
    i__1 = igradu;
    for (j = igradl; j <= i__1; ++j)
    {
      sum = 0.;
      i__2 = m;
      for (i__ = 1; i__ <= i__2; ++i__)
      {
        ++k;
        sum += xjac[k] * f[i__];
      }
      work[j] = sum;
      erl2 += sum * sum;
    }
    erl2 = (double)sqrt(erl2);
    if (ijac > 0)
      goto L190;
    if (erl2 <= delta2)
      *infer += 4;
    if (erl2 <= cons2)
      isw = 2;
L190:
    l = 0;
    is = -m;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      is += m;
      js = 0;
      i__2 = i__;
      for (j = 1; j <= i__2; ++j)
      {
        ++l;
        sum = 0.;
        i__3 = m;
        for (k = 1; k <= i__3; ++k)
        {
          li = is + k;
          ++js;
          sum += xjac[li] * xjac[js];
        }
        xjtj[l] = sum;
      }
    }
    if (*infer > 0)
      goto L640;
    if (ieval >= maxfn)
      goto L590;
    if (iopt == 0)
      goto L240;
    k = 0;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      k += j;
      ical1j = iscal1 + j;
      work[ical1j] = xjtj[k];
    }
    goto L270;
L240:
    dnorm = 0.;
    k = 0;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      k += j;
      ical1j = iscal1 + j;
      work[ical1j] = (double)sqrt(xjtj[k]);
      dnorm += xjtj[k] * xjtj[k];
    }
    dnorm = 1. / (double)sqrt(dnorm);
    i__1 = iscalu;
    for (j = iscall; j <= i__1; ++j)
        work[j] = work[j] * dnorm * erl2;
L270:
    icount = 0;
L280:
    k = 0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = i__;
      for (j = 1; j <= i__2; ++j)
      {
        ++k;
        work[k] = xjtj[k];
      }
      ical1i = iscal1 + i__;
      work[k] += work[ical1i] * al;
      ielx1i = idelx1 + i__;
      irad1i = igrad1 + i__;
      work[ielx1i] = work[irad1i];
    }
L310:
    ier = ownBPL_ash1r_c(&work[1], &work[idelxl], &work[ixnewl], n);
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      ielx1j = idelx1 + j;
      inew1j = ixnew1 + j;
      work[ielx1j] = work[inew1j];
    }
    if (ier == 0)
      goto L330;
    ier = 0;
    if (ijac > 0)
      goto L110;
    if (ibad <= 0)
      goto L490;
    if (ibad >= 2)
      goto L630;
    goto L390;
L330:
    if (ibad != -99)
      ibad = 0;
L340:
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      ielx1j = idelx1 + j;
      inew1j = ixnew1 + j;
      work[inew1j] = x[j] - work[ielx1j];
    }
    (*func)(&work[ixnewl], m, n, &work[ifpl]);
    ++ieval;
    *ssq = 0.;
    i__1 = ifpu;
    for (i__ = ifpl; i__ <= i__1; ++i__)
      *ssq += work[i__] * work[i__];
    if (iter >= 0)
      goto L380;
    iter = 0;
    ssqold = *ssq;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      ifpl1i = ifpl1 + i__;
      f[i__] = work[ifpl1i];
    }
    goto L110;
L380:
    if (iopt == 0)
      goto L440;
    if (*ssq <= ssqold)
      goto L420;
L390:
    ++icount;
    al *= fosq;
    if (ijac == 0)
      goto L400;
    if (icount >= 4 || al > up)
      goto L410;
L400:
    if (al <= up)
      goto L280;
    if (ibad == 1)
      goto L630;
    goto L610;
L410:
    al /= fosqs4;
    goto L110;
L420:
    if (icount == 0)
      al /= fo;
    if (erl2x <= 0.)
      goto L430;
    g = erl2 / erl2x;
    if (erl2 < erl2x)
      al *= mymax(onesfo,g);
    if (erl2 > erl2x)
      al *= mymin(fo,g);
L430:
    al = mymax(al,prec);
L440:
    ++iter;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      inew1j = ixnew1 + j;
      x[j] = work[inew1j];
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      ifml1i = ifml1 + i__;
      work[ifml1i] = f[i__];
      ifpl1i = ifpl1 + i__;
      f[i__] = work[ifpl1i];
    }
    if (icount > 0 || ijac > 0)
      goto L60;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      ielx1j = idelx1 + j;
/* Computing MAX */
      r__3 = (r__2 = x[j], myabs(r__2));
      xdif = (r__1 = work[ielx1j], myabs(r__1)) / mymax(r__3,.1);
      if (xdif > relcon)
        goto L480;
    }
    *infer = 1;
L480:
    sqdif = (r__1 = *ssq - ssqold, myabs(r__1)) / mymax(ssqold,.1);
    if (sqdif <= eps)
      *infer += 2;
    goto L60;
L490:
    if (ibad < 0)
      goto L520;
    else
      if (ibad == 0)
        goto L500;
      else
        goto L540;
L500:
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      ibad1j = ixbad1 + j;
      xhold = work[ibad1j];
/* Computing MAX */
      r__2 = .1, r__3 = myabs(xhold);
      if ((r__1 = x[j] - xhold, myabs(r__1)) > relcon * mymax(r__2,r__3))
        goto L520;
    }
    goto L600;
L520:
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      ibad1j = ixbad1 + j;
      work[ibad1j] = x[j];
    }
    ibad = 1;
L540:
    if (iopt != 0)
      goto L570;
    k = 0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = i__;
      for (j = 1; j <= i__2; ++j)
      {
        ++k;
        work[k] = xjtj[k];
      }
      ical1i = iscal1 + i__;
      work[k] = 1.5 * (xjtj[k] + al * erl2 * work[ical1i]) + rel;
    }
    ibad = 2;
    goto L310;
L570:
    izero = 0;
    i__1 = iscalu;
    for (j = iscall; j <= i__1; ++j)
      if (work[j] <= 0.)
      {
        ++izero;
        work[j] = 1.;
      }
    if (izero < n)
      goto L280;
    ier = 38;
    goto L640;
L590:
    ++(ier);
L600:
    ++(ier);
L610:
    ++(ier);
    ++(ier);
L630:
    ier += 129;
    if (ier == 130)
      return ier;
L640:
    work[1] = erl2 + erl2;
    work[2] = (double) ieval;
    *ssq = ssqold;
    g = sig;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      ielx1j = idelx1 + j;
      xhold = (r__1 = work[ielx1j], myabs(r__1));
      if (xhold <= 0.)
        goto L650;
// Computing MIN
// Computing MAX
      r__5 = .1, r__6 = (r__1 = x[j], myabs(r__1));
      r__4 = xhold / mymax(r__5,r__6);
      r__2 = g, r__3 = -log10(r__4);
      g = mymin(r__2,r__3);
L650:;
    }
    work[3] = g;
    work[4] = al;
    work[5] = (double) iter;
    if (ixjac == m)
      return ier;
    li = n - 1;
    k = li * m;
    ir = li * ixjac;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
      if (j == n)
        continue;
      i__2 = m;
      for (i__ = 1; i__ <= i__2; ++i__)
      {
        iri = ir + i__;
        ki = k + i__;
        xjac[iri] = xjac[ki];
      }
      k -= m;
      ir -= ixjac;
    }
    return ier;
}
