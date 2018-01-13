/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     1988                                                */
/*      Modified:    29 October 2004 - added QL eigen (v2)               */
/*      Revision:    2.0.00                                              */
/*      Purpose:     eigenvector calculation                             */
/*      Authors:                                                         */
/*        1.Numercal Recipes in C (www.nr.com) created it (1988)         */
/*        2.Matthew Turk typed this in (22-3-1991)                       */
/*        3.Thad Starner added error handling (6-10-1991)                */
/*        4.Victor Kuznetsov adapted for C++ CTDL project (15-5-1997)    */
/*        5.Ivan Matveev adapted this (23-3-2001)                        */
/*        6.NIVC MGU (http://www.srcc.msu.su)                            */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <malloc.h>
#include <memory.h>
#include "../ownbpl.h"

/*
   Use Jacobi Transformations of a symetric matrix to obtain the
   normalized eigenvectors and eigen values
   -- Taken from NUMERICAL RECIPES IN C (pp 364-366) --
*/

#define EPS 0.000000001

// single Jacobi rotation
#define ROTATE(a, i, j, k, l) \
  {\
    g = a[i*size+j];\
    h = a[k*size+l];\
    a[i*size+j] = g-s*(h+g*tau);\
    a[k*size+l] = h+s*(g-h*tau);\
  }

// Computes all eigenvalues and eigenvectors of a real symmetric matrix A
// On output, elements of A above the diagonal are destroyed,
// D returns the eigenvalues of A, and V is a matrix whose columns are the
// eigenvectors of A.
MIA_RESULT_CODE BPL_Matrix_EigenJacobiRotations(
  double* Matr, // IN: symmetric matrix which eigens are calculated
                //     elements above the diagonal are destroyed
  int size,     // IN: size of this matrix
  double* vals, // OUT: eigevalues of A
  double* Vect, // OUT: eigenvectors of A (in columns)
  int* sweeps,  // number of Jacobi rotation sweeps
                // IN: thresholding the stop. OUT: reached in fact
  double* thr,  // sum of the off-diagonal elements
                // IN: thresholding the stop. OUT: reached in fact
  int* pnrot)   // OUT: elementary rotations performed
{
  int j,iq,ip,ip_times_n,i_times_n1,i,SweepCountDown,nrot;
  double tresh,theta,tau,t,s,h,g,c,*b,*z,OffDiagSum;

    // allocate the temporarily used space
    b = (double*)malloc(size*sizeof(double));
    z = (double*)malloc(size*sizeof(double));
    // clear z
    memset(z,0,size*sizeof(z[0]));
    // clear Vect
    memset(Vect,0,size*size*sizeof(Vect[0]));
    for (i=0,i_times_n1=0;  i<size;  i++,i_times_n1+=(size+1))
    {
      // Initialize Vect as the identity matrix
      Vect[i_times_n1] = 1.;
      // Initailize b and vals to be diagonal of Matr
      b[i] = vals[i] = Matr[i_times_n1];
    }
    // no rotations made yet
    nrot = 0;
    // get the maximum number of rotations
    SweepCountDown = *sweeps;
    // while the given number of sweeps...
    while (SweepCountDown--)
    {
      // Summarize the off-diagonal elements
      OffDiagSum = 0.;
      for (ip_times_n=0,ip=0;  ip<size-1;  ip++,ip_times_n+=size)
        for (iq=ip+1;iq<size;iq++)
  	      OffDiagSum += fabs(Matr[ip_times_n+iq]);
      // if we have converged free the working vectors and return
      if (OffDiagSum<=*thr)
  	  {
  	    free(b);
  	    free(z);
        *thr = OffDiagSum;
        *sweeps -= SweepCountDown;
        *pnrot = nrot;
        return ERR_OK;
  	  }
      // on the first three iteration cut away substantially big elements
      // after this cut away each non-zero element
      tresh = (i<4) ? (.2*OffDiagSum/(size*size)) : 0.;
      for (ip_times_n=0,ip=0; ip<size-1; ip++,ip_times_n+=size)
  	  {
  	    for (iq=ip+1;iq<size;++iq)
  	    {
  	      g = 100.0*fabs(Matr[ip_times_n + iq]);
          // After four sweeps, skip the rotation if the off-diagonal
          // element is small
          // This test is taken directly from the text and looks
          // a little suspect to me...
          if (i > 3 && g < EPS)
            Matr[ip_times_n + iq] = 0.;
  	      else
            if (fabs(Matr[ip_times_n+iq]) <= tresh)
              continue;
            else
    	  	  {
    		      h = vals[iq]-vals[ip];
    		      if (g < EPS)
    		        t = (fabs(Matr[ip_times_n+iq]) > EPS) ?
                         (Matr[ip_times_n+iq]/h) : 0.;
    		      else
    		      {
    		        theta = (fabs(h) < EPS) ? 0.:.5*h/Matr[ip_times_n+iq];
    		        if (theta >= 0.)
    			        t = 1./(theta+sqrt(1.+theta*theta));
                else
    		          t = 1./(theta-sqrt(1.+theta*theta));
    		      }
    		      c = 1./sqrt(1.+t*t);
    		      s = t*c;
    		      tau = s/(1.+c);
    		      h = t*Matr[ip_times_n+iq];
    		      z[ip] -= h;
    		      z[iq] += h;
    		      vals[ip] -= h;
    		      vals[iq] += h;
    		      Matr[ip_times_n+iq] = 0.;
              // rotations of matrix
              // case of rotations 1<=j<p
    		      for(j=0;j<ip;j++)
    		        ROTATE(Matr,j,ip,j,iq);
              // case of rotations p<j<q
    		      for(j=ip+1;j<iq;j++)
    		        ROTATE(Matr,ip,j,j,iq);
              // case of rotations q<j<=n
    		      for(j=iq+1;j<size;j++)
    		        ROTATE(Matr,ip,j,iq,j);
              // rotations of eigen
    		      for(j=0;j<size;j++)
    		        ROTATE(Vect,j,ip,j,iq);
    		      nrot++;
    		    }
  	    }
  	  }
      for (ip=0;ip<size;++ip)
  	  {
  	    b[ip] += z[ip];
  	    vals[ip] = b[ip];
  	    z[ip] = 0.;
  	  }
    }
    // we come here if could not converge in a given number of rotations
    OffDiagSum = 0.;
    for(ip_times_n=0,ip=0;ip<size-1;ip++,ip_times_n+=size)
      for(iq=ip+1;iq<size;iq++)
        OffDiagSum += fabs(Matr[ip_times_n + iq]);
    free(b);
    free(z);
    *thr = OffDiagSum;
    *sweeps -= SweepCountDown;
    *pnrot = nrot;
  	return MIA_LIN_JACOBI_FAILED_CONVERGE;
}

// Eigen values by QL with shift - from http://www.srcc.msu.su
#define v_ref(a_1,a_2) a[(a_2)*n + (a_1)]
#define d_sign(a,b) (((b)>=0)?(((a)>=0)?(a):(-(a))):(((a)>= 0)?(-(a)):(a)))
#define sys051 1.1107652e-16

// solution of common eigen value problem Ax=lx by QL algorithm
// where A is real symmetric
// see www.srcc.msu.su/num_anal/lib_na/cat/ae_htm_c/aeh1r_c.htm
MIA_RESULT_CODE BPL_Matrix_Eigen_Symm(
  double *a,        // IN:  source matrix / OUT: eigenvector matrix
  double *ev,       // OUT: eigen values
  int n,            // IN:  size of matrix
  void* extbuftmp,  // IN:  temporary buffer
  int nTMP)         // IN:  size of temporary buffer
{
  int b,c__,d__,e,f,g,m,o,p,q,r__,s,t,u;
  double d__1,d__2,h__,i__,j,k,l,w,x,y,z__,ba,bb,bc = 0,bd,*rab1;

    // check args
    if ((a==NULL)||(ev==NULL)||(extbuftmp==NULL))
      return ERR_GEN_NULLPOINTER;
    if (n<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // allocate
    if (nTMP<n*(int)sizeof(a[0]))
      return ERR_GEN_NOMEMORY;
    rab1 = (double*)extbuftmp;
    //
    --rab1;
    --ev;
    a -= (n + 1);
    if (n != 1)
    {
      for (f = 2; f <= n; ++f)
      {
        b = n + 2 - f;
        e = b - 1;
        j = 0.;
        l = 0.;
        if (e < 2)
        {
          rab1[b] = v_ref(b, e);
          ev[b] = j;
          continue;
        }
        for (d__ = 1; d__ <= e; ++d__)
          l += (d__1 = v_ref(b, d__), MIA_abs(d__1));
        if (l == 0.)
        {
          rab1[b] = v_ref(b, e);
          ev[b] = j;
          continue;
        }
        for (d__ = 1; d__ <= e; ++d__)
        {
          v_ref(b, d__) = v_ref(b, d__) / l;
          j += v_ref(b, d__) * v_ref(b, d__);
        }
        h__ = v_ref(b, e);
        d__1 = sqrt(j);
        i__ = -d_sign(d__1,h__);
        rab1[b] = l * i__;
        j -= h__ * i__;
        v_ref(b, e) = h__ - i__;
        h__ = 0.;
        for (c__ = 1; c__ <= e; ++c__)
        {
          v_ref(c__, b) = v_ref(b, c__) / j;
          i__ = 0.;
          for (d__ = 1; d__ <= c__; ++d__)
            i__ += v_ref(c__, d__) * v_ref(b, d__);
          g = c__ + 1;
          if (e >= g)
          {
            for (d__ = g; d__ <= e; ++d__)
              i__ += v_ref(d__, c__) * v_ref(b, d__);
          }
          rab1[c__] = i__ / j;
          h__ += rab1[c__] * v_ref(b, c__);
        }
        k = h__ / (j + j);
        for (c__ = 1; c__ <= e; ++c__)
        {
          h__ = v_ref(b, c__);
          i__ = rab1[c__] - k * h__;
          rab1[c__] = i__;
          for (d__ = 1; d__ <= c__; ++d__)
            v_ref(c__, d__) = v_ref(c__, d__) - h__ * rab1[d__] - i__ * v_ref(b, d__);
        }
        ev[b] = j;
      }
    }
    ev[1] = 0.;
    rab1[1] = 0.;
    for (b = 1; b <= n; ++b)
    {
      e = b - 1;
      if (ev[b] != 0.)
      {
        for (c__ = 1; c__ <= e; ++c__)
        {
          i__ = 0.;
          for (d__ = 1; d__ <= e; ++d__)
            i__ += v_ref(b, d__) * v_ref(d__, c__);
          for (d__ = 1; d__ <= e; ++d__)
            v_ref(d__, c__) = v_ref(d__, c__) - i__ * v_ref(d__, b);
        }
      }
      ev[b] = v_ref(b, b);
      v_ref(b, b) = 1.;
      if (e < 1)
        continue;
      for (c__ = 1; c__ <= e; ++c__)
      {
        v_ref(b, c__) = 0.;
        v_ref(c__, b) = 0.;
      }
    }
    if (n == 1)
      return ERR_OK;
    for (m = 2; m <= n; ++m)
      rab1[m - 1] = rab1[m];
    y = 0.;
    w = 0.;
    rab1[n] = 0.;
    for (q = 1; q <= n; ++q)
    {
      o = 0;
      ba = sys051 * ((d__1 = ev[q], MIA_abs(d__1)) + (d__2 = rab1[q], MIA_abs(d__2)));
      if (w < ba)
        w = ba;
      for (r__ = q; r__ <= n; ++r__)
        if ((d__1 = rab1[r__], MIA_abs(d__1)) <= w)
          break;
      if (r__ == q)
      {
        ev[q] += y;
        continue;
      }
      do
      {
        if (o == 30)
          return MIA_LIN_JACOBI_FAILED_CONVERGE;// q is a number of vectors calculated OK;
        ++o;
        t = q + 1;
        z__ = ev[q];
        bb = (ev[t] - z__) / (rab1[q] * 2.);
        if (MIA_abs(bb) > 1.)
        { // Computing 2nd power
          d__1 = 1 / bb;
          bc = MIA_abs(bb) * sqrt(d__1 * d__1 + 1.);
        }
        if (MIA_abs(bb) <= 1.)
          bc = sqrt(bb * bb + 1.);
        ev[q] = rab1[q] / (bb + d_sign(bc,bb));
        ba = z__ - ev[q];
        for (m = t; m <= n; ++m)
          ev[m] -= ba;
        y += ba;
        bb = ev[r__];
        x = 1.;
        bd = 0.;
        u = r__ - q;
        for (s = 1; s <= u; ++s)
        {
          m = r__ - s;
          z__ = x * rab1[m];
          ba = x * bb;
          if (MIA_abs(bb) < (d__1 = rab1[m], MIA_abs(d__1)))
          {
            x = bb / rab1[m];
            bc = sqrt(x * x + 1.);
            rab1[m + 1] = bd * rab1[m] * bc;
            bd = 1. / bc;
            x *= bd;
          }
          else
          {
            x = rab1[m] / bb;
            bc = sqrt(x * x + 1.);
            rab1[m + 1] = bd * bb * bc;
            bd = x / bc;
            x = 1. / bc;
          }
          bb = x * ev[m] - bd * z__;
          ev[m + 1] = ba + bd * (x * z__ + bd * ev[m]);
          for (p = 1; p <= n; ++p)
          {
            ba = v_ref(p, m + 1);
            v_ref(p, m + 1) = bd * v_ref(p, m) + x * ba;
            v_ref(p, m) = x * v_ref(p, m) - bd * ba;
          }
        }
        rab1[q] = bd * bb;
        ev[q] = x * bb;
      }
      while ((d__1 = rab1[q], MIA_abs(d__1)) > w);
      ev[q] += y;
    }
    for (s = 2; s <= n; ++s)
    {
      m = s - 1;
      p = m;
      bb = ev[m];
      for (o = s; o <= n; ++o)
      {
        if (ev[o] >= bb)
          continue;
        p = o;
        bb = ev[o];
      }
      if (p == m)
        continue;
      ev[p] = ev[m];
      ev[m] = bb;
      for (o = 1; o <= n; ++o)
      {
        bb = v_ref(o, m);
        v_ref(o, m) = v_ref(o, p);
        v_ref(o, p) = bb;
      }
    }
    return ERR_OK;
}
#undef v_ref

#define r_sign(a,b) (((b)>=0)?(((a)>=0)?(a):(-(a))):(((a)>= 0)?(-(a)):(a)))
//#define sys051__ 1.192093e-7
#define sys051__ 1.0e-16
//#define MIA_abs(a) (((a)>0)?(a):(-(a)))

// solution of generalized eigen problem Ax=lBx
// where A,B - real symmetric, B is positively defined
// see www.srcc.msu.su/num_anal/lib_na/cat/ag_htm_c/agh0r_c.htm
MIA_RESULT_CODE BPL_Matrix_EigenGeneralized_Symm(
  double *a,        // IN:  left-size matrix (A)
  double *b,        // IN:  right-side matrix (B)
  double *v,        // OUT: eigen vectors
  double *ev,       // OUT: eigen values
  int n,            // IN:  dimensionality
  void* extbuftmp,  // IN:  temporary buffer
  int nTMP)         // IN:  size of temporary buffer
{
  int c__, d__, e, f, g, k, l, m, o, p, q,x, y, z__;
  int ba, bb, bc, bd, be, bo, bp, bq, br, bs;
  int nk, nl, nm, nx, nx1, nba, nbb;
  double r__1,r__2,bt,i__, j = 0,r__, s, t, u, w;
  double bf, bg, bh, bi, bj, bk, bl, bm;
  double *rab;

    // check args
    if ((a==NULL)||(b==NULL)||(v==NULL)||(ev==NULL)||(extbuftmp==NULL))
      return ERR_GEN_NULLPOINTER;
    if (n<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // allocate
    if (nTMP<2*n*(int)sizeof(a[0]))
      return ERR_GEN_NOMEMORY;
    rab = (double*)extbuftmp;
    // Parameter adjustments
    --rab;
    --ev;
    v -= n + 1;
    b -= n + 1;
    a -= n + 1;
    // Function Body
    for (c__ = 1; c__ <= n; ++c__)
    {
      f = c__ - 1;
      for (d__ = c__; d__ <= n; ++d__)
      {
        i__ = b[d__*n + c__];
        if (c__ != 1)
          for (e = 1; e <= f; ++e)
            i__ -= b[e*n+c__]*b[e*n+d__];
        if (d__ != c__)
          b[c__*n + d__] = i__ / j;
        else
        {
          if (i__ <= 0.)
            return n*7+1;
          j = sqrt(i__);
          rab[c__] = j;
        }
      }
    }
    for (c__ = 1; c__ <= n; ++c__)
    {
      f = c__ - 1;
      j = rab[c__];
      for (d__ = c__; d__ <= n; ++d__)
      {
        i__ = a[d__*n + c__];
        if (c__ != 1)
          for (e = 1; e <= f; ++e)
            i__ -= b[e*n + c__] * a[e*n + d__];
        a[c__*n + d__] = i__ / j;
      }
    }
    for (d__ = 1; d__ <= n; ++d__)
    {
      g = d__ - 1;
      for (c__ = d__; c__ <= n; ++c__)
      {
        i__ = a[d__*n + c__];
        if (c__ != d__)
        {
          f = c__ - 1;
          for (e = d__; e <= f; ++e)
            i__ -= a[d__*n + e] * b[e*n + c__];
        }
        if (d__ != 1)
          for (e = 1; e <= g; ++e)
            i__ -= a[e*n + d__] * b[e*n + c__];
        a[d__*n + c__] = i__ / rab[c__];
      }
    }
    for (k = 1; k <= n; ++k)
      for (l = 1; l <= k; ++l)
        v[l*n + k] = a[l*n + k];
    for (p = 2; p <= n; ++p)
    {
      k = n + 2 - p;
      nk = n + k;
      o = k - 1;
      t = 0.;
      w = 0.;
      if (o < 2)
      {
        rab[nk] = v[o*n + k];
        ev[k] = t;
        continue;
      }
      for (m = 1; m <= o; ++m)
        w += (r__1 = v[m*n + k], MIA_abs(r__1));
      if (w == 0.)
      {
        rab[nk] = v[o*n + k];
        ev[k] = t;
        continue;
      }
      for (m = 1; m <= o; ++m)
      {
        v[m*n+k] /= w;
        t += v[m*n+k]*v[m*n+k];
      }
      r__ = v[o*n+k];
      r__1 = sqrt(t);
      s = -r_sign(r__1,r__);
      rab[nk] = w * s;
      t -= r__ * s;
      v[o*n+k] = r__ - s;
      r__ = 0.;
      for (l = 1; l <= o; ++l)
      {
        nl = n + l;
        v[k*n+l] = v[l*n+k] / t;
        s = 0.;
        for (m = 1; m <= l; ++m)
          s += v[m*n +l] * v[m*n+k];
        q = l + 1;
        if (o >= q)
          for (m = q; m <= o; ++m)
            s += v[l*n+m] * v[m*n+k];
        rab[nl] = s / t;
        r__ += rab[nl] * v[l*n+k];
      }
      u = r__ / (t + t);
      for (l = 1; l <= o; ++l)
      {
        nl = n + l;
        r__ = v[l*n+k];
        s = rab[nl] - u * r__;
        rab[nl] = s;
        for (m = 1; m <= l; ++m)
        {
          nm = n + m;
          v[m*n+l] = v[m*n+l] - r__ * rab[nm] - s * v[m*n+k];
        }
      }
      ev[k] = t;
    }
    ev[1] = 0.;
    rab[n + 1] = 0.;
    for (k = 1; k <= n; ++k)
    {
      o = k - 1;
      if (ev[k] != 0.)
      {
        for (l = 1; l <= o; ++l)
        {
          s = 0.;
          for (m = 1; m <= o; ++m)
            s += v[m*n+k] * v[l*n+m];
          for (m = 1; m <= o; ++m)
            v[l*n+m] -= s * v[k*n+m];
        }
      }
      ev[k] = v[k*n+k];
      v[k*n+k] = 1.;
      if (o >= 1)
        for (l = 1; l <= o; ++l)
          v[k*n+l] = v[l*n+k] = 0.;
    }
    for (x = 2; x <= n; ++x)
    {
      nx1 = n + x - 1;
      nx = n + x;
      rab[nx1] = rab[nx];
    }
    bh = 0.;
    bf = 0.;
    rab[2*n] = 0.;
    for (ba = 1; ba <= n; ++ba)
    {
      nba = n + ba;
      y = 0;
      bj = sys051__*((r__1 = ev[ba], MIA_abs(r__1)) + (r__2 = rab[nba], MIA_abs(r__2)));
      if (bf < bj)
        bf = bj;
      for (bb = ba; bb <= n; ++bb)
      {
        nbb = n + bb;
        if ((r__1 = rab[nbb], MIA_abs(r__1)) <= bf)
          break;
      }
      if (bb == ba)
      {
        ev[ba] += bh;
        continue;
      }
      do
      {
        if (y == 30)
          return ba;
        ++y;
        bd = ba + 1;
        bi = ev[ba];
        bk = (ev[bd] - bi) / (rab[nba] * 2.);
        if (MIA_abs(bk) > 1.)
        { // Computing 2nd power
          r__1 = 1 / bk;
          bl = MIA_abs(bk) * sqrt(r__1 * r__1 + 1.);
        }
        else
          bl = sqrt(bk * bk + 1.);
        ev[ba] = rab[nba] / (bk + r_sign(bl,bk));
        bj = bi - ev[ba];
        for (x = bd; x <= n; ++x)
          ev[x] -= bj;
        bh += bj;
        bk = ev[bb];
        bg = 1.;
        bm = 0.;
        be = bb - ba;
        for (bc = 1; bc <= be; ++bc)
        {
          x = bb - bc;
          nx = n + x;
          nx1 = n + x + 1;
          bi = bg * rab[nx];
          bj = bg * bk;
          if (MIA_abs(bk) >= (r__1 = rab[nx], MIA_abs(r__1)))
          {
            bg = rab[nx] / bk;
            bl = sqrt(bg * bg + 1.);
            rab[nx1] = bm * bk * bl;
            bm = bg / bl;
            bg = 1. / bl;
          }
          else
          {
            bg = bk / rab[nx];
            bl = sqrt(bg * bg + 1.);
            rab[nx1] = bm * rab[nx] * bl;
            bm = 1. / bl;
            bg *= bm;
          }
          bk = bg * ev[x] - bm * bi;
          ev[x + 1] = bj + bm * (bg * bi + bm * ev[x]);
          for (z__ = 1; z__ <= n; ++z__)
          {
            bj = v[(x+1)*n +z__];
            v[(x+1)*n+z__] = bm * v[x*n +z__] + bg * bj;
            v[x*n+z__] = bg * v[x*n+z__] - bm * bj;
          }
        }
        rab[nba] = bm * bk;
        ev[ba] = bg * bk;
      }
      while ((r__1 = rab[nba], MIA_abs(r__1)) > bf);  // ->do
      ev[ba] += bh;
    }
    for (bc = 2; bc <= n; ++bc)
    {
      x = bc - 1;
      z__ = x;
      bk = ev[x];
      for (y = bc; y <= n; ++y)
        if (ev[y] < bk)
        {
          z__ = y;
          bk = ev[y];
        }
      if (z__ != x)
      {
        ev[z__] = ev[x];
        ev[x] = bk;
        for (y = 1; y <= n; ++y)
        {
          bk = v[x*n+y];
          v[x*n+y] = v[z__*n+y];
          v[z__*n+y] = bk;
        }
      }
    }
    for (bp = 1; bp <= n; ++bp)
    {
      for (bs = 1; bs <= n; ++bs)
      {
        bo = n + 1 - bs;
        br = bo + 1;
        bt = v[bp*n+bo];
        if (bo != n)
          for (bq = br; bq <= n; ++bq)
            bt -= b[bo*n + bq] * v[bp*n+bq];
        v[bp*n+bo] = bt / rab[bo];
      }
    }
    return 0;
}

#define ot116 1.1605

staticFunc int ownBPL_am12c_c(
  double *ajr,
  double *aji,
  double *ajp1r,
  double *ajp1i,
  double *c__,
  double *sr,
  double *si)
{
  double r__;

    if (*ajp1r == 0. && *ajp1i == 0.)
    {
      *c__ = 1.;
      *sr = 0.;
      *si = 0.;
      return 0;
    }
    if (*ajr == 0. && *aji == 0.)
    {
      *c__ = 0.;
      *sr = 1.;
      *si = 0.;
      return 0;
    }
    r__ = (double)sqrt(*ajr * *ajr + *aji * *aji);
    *c__ = r__;
    *sr = (*ajr * *ajp1r + *aji * *ajp1i) / r__;
    *si = (*ajr * *ajp1i - *aji * *ajp1r) / r__;
    r__ = 1. / (double)sqrt(*c__ * *c__ + *sr * *sr + *si * *si);
    *c__ *= r__;
    *sr *= r__;
    *si *= r__;
    return 0;
}

staticFunc int ownBPL_am12r_c(
  double *aj,
  double *ajp1,
  double *uj,
  double *ujp1,
  double *vj,
  double *vjp1)
{
  double r__, s;

    if (*ajp1 == 0.)
      *uj = 0.;
    else
    {
      s = MIA_abs(*aj) + MIA_abs(*ajp1);
      *uj = *aj / s;
      *ujp1 = *ajp1 / s;
      r__ = (double)sqrt(*uj * *uj + *ujp1 * *ujp1);
      if (*uj < 0.)
        r__ = -r__;
      *vj = -(*uj + r__) / r__;
      *vjp1 = -(*ujp1) / r__;
      *uj = 1.;
      *ujp1 = *vjp1 / *vj;
    }
    return 0;
}

staticFunc int ownBPL_am13r_c(
  double *aj,
  double *ajp1,
  double *ajp2,
  double *uj,
  double *ujp1,
  double *ujp2,
  double *vj,
  double *vjp1,
  double *vjp2)
{
  double r__, s, rd;

    if (*ajp1 == 0. && *ajp2 == 0.)
      *uj = 0.;
    else
    {
      s = MIA_abs(*aj) + MIA_abs(*ajp1) + MIA_abs(*ajp2);
      rd = 1. / s;
      *uj = *aj * rd;
      *ujp1 = *ajp1 * rd;
      *ujp2 = *ajp2 * rd;
      r__ = (double)sqrt(*uj * *uj + *ujp1 * *ujp1 + *ujp2 * *ujp2);
      if (*uj < 0.)
        r__ = -r__;
      rd = 1. / r__;
      *vj = -(*uj + r__) * rd;
      *vjp1 = -(*ujp1) * rd;
      *vjp2 = -(*ujp2) * rd;
      *uj = 1.;
      *ujp1 = *vjp1 / *vj;
      *ujp2 = *vjp2 / *vj;
    }
    return 0;
}

#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]
staticFunc int ownBPL_afg3r_c(
  double *a,
  double *b,
  double *z__,
  int *n,
  int *m)
{
  int a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2, i__3;
  int i__, j, k, l, k1, l1, wantx, lb, nk1, nm1, nm2;
  double r__1, r__, s, t, u1, u2, v1, v2, sd, rho;

    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    b_dim1 = *n;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    /* Function Body */
    wantx = 0;
    if (*m != 1)
    {
      wantx = 1;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
      {
        i__2 = *n;
        for (j = 1; j <= i__2; ++j)
          z___ref(i__, j) = 0.;
        z___ref(i__, i__) = 1.;
      }
    }
    nm1 = *n - 1;
    if (*n <= 1)
      return 0;
    i__1 = nm1;
    for (l = 1; l <= i__1; ++l)
    {
      l1 = l + 1;
      s = 0.;
      i__2 = *n;
      for (i__ = l1; i__ <= i__2; ++i__)
        s += (r__1 = b_ref(i__, l), MIA_abs(r__1));
      if (s == 0.)
        continue;
      s += (r__1 = b_ref(l, l), MIA_abs(r__1));
      r__ = 0.;
      sd = 1. / s;
      i__2 = *n;
      for (i__ = l; i__ <= i__2; ++i__)
      {
        b_ref(i__, l) = b_ref(i__, l) * sd;
        /* Computing 2nd power */
        r__1 = b_ref(i__, l);
        r__ += r__1 * r__1;
      }
      r__ = (double)sqrt(r__);
      if (b_ref(l, l) < 0.)
        r__ = -r__;
      b_ref(l, l) = b_ref(l, l) + r__;
      rho = r__ * b_ref(l, l);
      sd = 1. / rho;
      i__2 = *n;
      for (j = l1; j <= i__2; ++j)
      {
        t = 0.;
        i__3 = *n;
        for (i__ = l; i__ <= i__3; ++i__)
          t += b_ref(i__, l) * b_ref(i__, j);
        t = -t * sd;
        i__3 = *n;
        for (i__ = l; i__ <= i__3; ++i__)
          b_ref(i__, j) = b_ref(i__, j) + t * b_ref(i__, l);
      }
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        t = 0.;
        i__3 = *n;
        for (i__ = l; i__ <= i__3; ++i__)
          t += b_ref(i__, l) * a_ref(i__, j);
        t = -t * sd;
        i__3 = *n;
        for (i__ = l; i__ <= i__3; ++i__)
          a_ref(i__, j) = a_ref(i__, j) + t * b_ref(i__, l);
      }
      b_ref(l, l) = -s * r__;
      i__2 = *n;
      for (i__ = l1; i__ <= i__2; ++i__)
        b_ref(i__, l) = 0.;
    }
    if (*n <= 2)
      return 0;
    nm2 = *n - 2;
    i__1 = nm2;
    for (k = 1; k <= i__1; ++k)
    {
      k1 = k + 1;
      nk1 = *n - k1;
      i__2 = nk1;
      for (lb = 1; lb <= i__2; ++lb)
      {
        l = *n - lb;
        l1 = l + 1;
        ownBPL_am12r_c(&a_ref(l, k), &a_ref(l1, k), &u1, &u2, &v1, &v2);
        if (u1 != 1.)
          goto L16;
        i__3 = *n;
        for (j = k; j <= i__3; ++j)
        {
          t = a_ref(l, j) + u2 * a_ref(l1, j);
          a_ref(l, j) = a_ref(l, j) + t * v1;
          a_ref(l1, j) = a_ref(l1, j) + t * v2;
        }
        a_ref(l1, k) = 0.;
        i__3 = *n;
        for (j = l; j <= i__3; ++j)
        {
          t = b_ref(l, j) + u2 * b_ref(l1, j);
          b_ref(l, j) = b_ref(l, j) + t * v1;
          b_ref(l1, j) = b_ref(l1, j) + t * v2;
        }
L16:
        ownBPL_am12r_c(&b_ref(l1, l1), &b_ref(l1, l), &u1, &u2, &v1, &v2);
        if (u1 != 1.)
          continue;
        i__3 = l1;
        for (i__ = 1; i__ <= i__3; ++i__)
        {
          t = b_ref(i__, l1) + u2 * b_ref(i__, l);
          b_ref(i__, l1) = b_ref(i__, l1) + t * v1;
          b_ref(i__, l) = b_ref(i__, l) + t * v2;
        }
        b_ref(l1, l) = 0.;
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__)
        {
          t = a_ref(i__, l1) + u2 * a_ref(i__, l);
          a_ref(i__, l1) = a_ref(i__, l1) + t * v1;
          a_ref(i__, l) = a_ref(i__, l) + t * v2;
        }
        if (! wantx)
          continue;
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__)
        {
          t = z___ref(i__, l1) + u2 * z___ref(i__, l);
          z___ref(i__, l1) = z___ref(i__, l1) + t * v1;
          z___ref(i__, l) = z___ref(i__, l) + t * v2;
        }
      }
    }
    return 0;
}
#undef z___ref
#undef b_ref
#undef a_ref

#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]
staticFunc int ownBPL_aft0r_c(
  double *a,
  double *b,
  double *z__,
  int *n,
  const int *m,
  int *ierr)
{
  int a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2;
  int iter,i__, j, k,k1, k2, m1, m2, k3, l1,wantx,lb,km1,mid;
  int l=0,lor1=0,morn=0;
  double r__1,epsa, epsb,t, anorm, bnorm, const__;
  double a10, a11, a21, b11, a22, b12, a34, b34, a43;
  double a44, b22, a20, a33, b33,ani, bni,eps,b44, a12, a30;
  double old1=0.,old2=0.,v1=0.,v2=0.,u1=0.,u2=0.,v3=0.,u3=0.;

    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    b_dim1 = *n;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    /* Function Body */
    *ierr = 0;
    eps = sys051__;
    wantx = 0;
    if (*m == 0)
      wantx = 1;
    anorm = 0.;
    bnorm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      ani = 0.;
      if (i__ != 1)
        ani = (r__1 = a_ref(i__, i__ - 1), MIA_abs(r__1));
      bni = 0.;
      i__2 = *n;
      for (j = i__; j <= i__2; ++j)
      {
        ani += (r__1 = a_ref(i__, j), MIA_abs(r__1));
        bni += (r__1 = b_ref(i__, j), MIA_abs(r__1));
      }
      if (ani > anorm)
        anorm = ani;
      if (bni > bnorm)
        bnorm = bni;
    }
    if (anorm == 0.)
      anorm = eps;
    if (bnorm == 0.)
      bnorm = eps;
    epsa = eps * anorm;
    epsb = eps * bnorm;
    m1 = *n;
    iter = 0;
L3:
    if (m1 <= 2)
      return 0;
    i__1 = m1;
    for (lb = 1; lb <= i__1; ++lb)
    {
      l = m1 + 1 - lb;
      if (l == 1)
        goto L6;
      if ((r__1 = a_ref(l, l - 1), MIA_abs(r__1)) <= epsa)
        break;
    }
L5:
    a_ref(l, l - 1) = 0.;
    if (l < m1 - 1)
      goto L6;
    m1 = l - 1;
    iter = 0;
    goto L3;
L6:
    if ((r__1 = b_ref(l, l), MIA_abs(r__1)) > epsb)
      goto L9;
    b_ref(l, l) = 0.;
    l1 = l + 1;
    ownBPL_am12r_c(&a_ref(l, l), &a_ref(l1, l), &u1, &u2, &v1, &v2);
    if (u1 != 1.)
      goto L8;
    i__1 = *n;
    for (j = l; j <= i__1; ++j)
    {
      t = a_ref(l, j) + u2 * a_ref(l1, j);
      a_ref(l, j) = a_ref(l, j) + t * v1;
      a_ref(l1, j) = a_ref(l1, j) + t * v2;
      t = b_ref(l, j) + u2 * b_ref(l1, j);
      b_ref(l, j) = b_ref(l, j) + t * v1;
      b_ref(l1, j) = b_ref(l1, j) + t * v2;
    }
L8:
    l = l1;
    goto L5;
L9:
    m2 = m1 - 1;
    l1 = l + 1;
    const__ = .75f;
    ++iter;
    if (iter == 1)
      goto L10;
    if ((r__1 = a_ref(m1, m1 - 1), MIA_abs(r__1)) < const__ * old1)
      goto L10;
    if ((r__1 = a_ref(m1 - 1, m1 - 2), MIA_abs(r__1)) < const__ * old2)
      goto L10;
    if (iter == 10)
      goto L11;
    if (iter > 30)
    {
      *ierr = m1 + 128;
      return 0;
    }
L10:
    b11 = b_ref(l, l);
    b22 = b_ref(l1, l1);
    if (MIA_abs(b22) < epsb)
      b22 = epsb;
    b33 = b_ref(m2, m2);
    if (MIA_abs(b33) < epsb)
      b33 = epsb;
    b44 = b_ref(m1, m1);
    if (MIA_abs(b44) < epsb)
      b44 = epsb;
    a11 = a_ref(l, l) / b11;
    a12 = a_ref(l, l1) / b22;
    a21 = a_ref(l1, l) / b11;
    a22 = a_ref(l1, l1) / b22;
    a33 = a_ref(m2, m2) / b33;
    a34 = a_ref(m2, m1) / b44;
    a43 = a_ref(m1, m2) / b33;
    a44 = a_ref(m1, m1) / b44;
    b12 = b_ref(l, l1) / b22;
    b34 = b_ref(m2, m1) / b44;
    a10 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) / a21 + 
            a12 - a11 * b12;
    a20 = a22 - a11 - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34;
    a30 = a_ref(l + 2, l1) / b22;
    goto L12;
L11:
    a10 = 0.;
    a20 = 1.;
    a30 = ot116;
L12:
    old1 = (r__1 = a_ref(m1, m1 - 1), MIA_abs(r__1));
    old2 = (r__1 = a_ref(m1 - 1, m1 - 2), MIA_abs(r__1));
    if (! wantx)
      lor1 = l;
    if (wantx)
      lor1 = 1;
    if (! wantx)
      morn = m1;
    if (wantx)
      morn = *n;
    i__1 = m2;
    for (k = l; k <= i__1; ++k)
    {
      mid = k != m2;
      k1 = k + 1;
      k2 = k + 2;
      k3 = k + 3;
      if (k3 > m1)
        k3 = m1;
      km1 = k - 1;
      if (km1 < l)
        km1 = l;
      if (k == l)
        ownBPL_am13r_c(&a10, &a20, &a30, &u1, &u2, &u3, &v1, &v2, &v3);
      if (k > l && k < m2)
        ownBPL_am13r_c(&a_ref(k, km1), &a_ref(k1, km1), &a_ref(k2, km1), &u1, &u2,&u3, &v1, &v2, &v3);
      if (k == m2)
        ownBPL_am12r_c(&a_ref(k, km1), &a_ref(k1, km1), &u1, &u2, &v1, &v2);
      if (u1 != 1.)
        goto L14;
      i__2 = morn;
      for (j = km1; j <= i__2; ++j)
      {
        t = a_ref(k, j) + u2 * a_ref(k1, j);
        if (mid)
          t += u3 * a_ref(k2, j);
        a_ref(k, j) = a_ref(k, j) + t * v1;
        a_ref(k1, j) = a_ref(k1, j) + t * v2;
        if (mid)
          a_ref(k2, j) = a_ref(k2, j) + t * v3;
        t = b_ref(k, j) + u2 * b_ref(k1, j);
        if (mid)
          t += u3 * b_ref(k2, j);
        b_ref(k, j) = b_ref(k, j) + t * v1;
        b_ref(k1, j) = b_ref(k1, j) + t * v2;
        if (mid)
          b_ref(k2, j) = b_ref(k2, j) + t * v3;
      }
      if (k == l)
        goto L14;
      a_ref(k1, k - 1) = 0.;
      if (mid)
        a_ref(k2, k - 1) = 0.;
L14:
      if (k == m2)
        goto L17;
      ownBPL_am13r_c(&b_ref(k2, k2), &b_ref(k2, k1), &b_ref(k2, k), &u1, &u2, &u3, &v1, &v2, &v3);
      if (u1 != 1.)
        goto L17;
      i__2 = k3;
      for (i__ = lor1; i__ <= i__2; ++i__)
      {
        t = a_ref(i__, k2) + u2 * a_ref(i__, k1) + u3 * a_ref(i__, k);
        a_ref(i__, k2) = a_ref(i__, k2) + t * v1;
        a_ref(i__, k1) = a_ref(i__, k1) + t * v2;
        a_ref(i__, k) = a_ref(i__, k) + t * v3;
        t = b_ref(i__, k2) + u2 * b_ref(i__, k1) + u3 * b_ref(i__, k);
        b_ref(i__, k2) = b_ref(i__, k2) + t * v1;
        b_ref(i__, k1) = b_ref(i__, k1) + t * v2;
        b_ref(i__, k) = b_ref(i__, k) + t * v3;
      }
      b_ref(k2, k) = 0.;
      b_ref(k2, k1) = 0.;
      if (! wantx)
        goto L17;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
      {
        t = z___ref(i__, k2) + u2 * z___ref(i__, k1) + u3 * z___ref(i__,k);
        z___ref(i__, k2) = z___ref(i__, k2) + t * v1;
        z___ref(i__, k1) = z___ref(i__, k1) + t * v2;
        z___ref(i__, k) = z___ref(i__, k) + t * v3;
      }
L17:
      ownBPL_am12r_c(&b_ref(k1, k1), &b_ref(k1, k), &u1, &u2, &v1, &v2);
      if (u1 != 1.)
        continue;
      i__2 = k3;
      for (i__ = lor1; i__ <= i__2; ++i__)
      {
        t = a_ref(i__, k1) + u2 * a_ref(i__, k);
        a_ref(i__, k1) = a_ref(i__, k1) + t * v1;
        a_ref(i__, k) = a_ref(i__, k) + t * v2;
        t = b_ref(i__, k1) + u2 * b_ref(i__, k);
        b_ref(i__, k1) = b_ref(i__, k1) + t * v1;
        b_ref(i__, k) = b_ref(i__, k) + t * v2;
      }
      b_ref(k1, k) = 0.;
      if (! wantx)
        continue;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
      {
        t = z___ref(i__, k1) + u2 * z___ref(i__, k);
        z___ref(i__, k1) = z___ref(i__, k1) + t * v1;
        z___ref(i__, k) = z___ref(i__, k) + t * v2;
      }
    }
    goto L3;
}
#undef z___ref
#undef b_ref
#undef a_ref

#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]
staticFunc int ownBPL_agt0r_c(
  double *a,
  double *b,
  double *v,
  double *alfr,
  double *alfi,
  double *beta,
  int *n,
  int *ierr)
{
  const int c__0 = 0;
//  const int c__41 = 41;
  int a_dim1, a_offset, b_dim1, b_offset, v_dim1, v_offset, i__1, i__2;
  int iret,i__,j,l,m,l1;
  int flip=0,k=0,mi=0,mr=0;
  double r__1,r__2, r__3, r__4,betm, epsa, epsb,c__,f,h__,r__, s, t;
  double alfm,anorm, bnorm,a11, a21, b11, a22, b12, b22, a12,ei,bn,er,ti;
  double sk, sl, xi, yi, zi, tr, xr, yr, zr, a11i, a12i, a22i, bdi, 
         a21i, a11r, a21r, a12r, a22r, bdr, ani, bni,eps, tkk, 
         tlk, tll, tkl, ssi,ssr;
  double u1=0.,v1=0.,v2=0.,u2=0.,an=0.,e=0.,d__=0.,szi=0.,szr=0.,cz=0.;
  double sqi=0.,sqr=0.,cq=0.,almr=0.,almi=0.,slr=0.,sli=0.,dr=0.,di=0.;
  double skr=0.,ski=0.,tkkr=0.,tkki=0.,tklr=0.,tkli=0.;
  double tlki=0.,tlli=0.,tllr=0.,tlkr=0.;

    /* Parameter adjustments */
    --beta;
    --alfi;
    --alfr;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    b_dim1 = *n;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    /* Function Body */
    *ierr = 0;
    eps = sys051__;
    ownBPL_aft0r_c(&a[a_offset], &b[b_offset], &v[v_offset], n, &c__0, ierr);
    anorm = 0.;
    bnorm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      ani = 0.;
      if (i__ != 1)
        ani = (r__1 = a_ref(i__, i__ - 1), MIA_abs(r__1));
      bni = 0.;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        ani += (r__1 = a_ref(i__, j), MIA_abs(r__1));
        bni += (r__1 = b_ref(i__, j), MIA_abs(r__1));
      }
      if (ani > anorm)
        anorm = ani;
      if (bni > bnorm)
        bnorm = bni;
    }
    if (anorm == 0.)
      anorm = eps;
    if (bnorm == 0.)
      bnorm = eps;
    epsa = eps * anorm;
    epsb = eps * bnorm;
    m = *n;
L3:
    if (m == 1)
      goto L4;
    if (a_ref(m, m - 1) != 0.)
      goto L5;
L4:
    alfr[m] = a_ref(m, m);
    if (b_ref(m, m) < 0.)
      alfr[m] = -alfr[m];
    beta[m] = (r__1 = b_ref(m, m), MIA_abs(r__1));
    alfi[m] = 0.;
    --m;
    goto L16;
L5:
    l = m - 1;
    if ((r__1 = b_ref(l, l), MIA_abs(r__1)) > epsb)
      goto L6;
    b_ref(l, l) = 0.;
    ownBPL_am12r_c(&a_ref(l, l), &a_ref(m, l), &u1, &u2, &v1, &v2);
    goto L12;
L6:
    if ((r__1 = b_ref(m, m), MIA_abs(r__1)) > epsb)
      goto L7;
    b_ref(m, m) = 0.;
    ownBPL_am12r_c(&a_ref(m, m), &a_ref(m, l), &u1, &u2, &v1, &v2);
    bn = 0.;
    goto L8;
L7:
    an = (r__1 = a_ref(l, l), MIA_abs(r__1)) + (r__2 = a_ref(l, m), MIA_abs(r__2)) +
            (r__3 = a_ref(m, l), MIA_abs(r__3)) + (r__4 = a_ref(m, m), MIA_abs(r__4));
    bn = (r__1 = b_ref(l, l), MIA_abs(r__1)) + (r__2 = b_ref(l, m), MIA_abs(r__2)) +
            (r__3 = b_ref(m, m), MIA_abs(r__3));
    f = 1. / an;
    a11 = a_ref(l, l) * f;
    a12 = a_ref(l, m) * f;
    a21 = a_ref(m, l) * f;
    a22 = a_ref(m, m) * f;
    f = 1. / bn;
    b11 = b_ref(l, l) * f;
    b12 = b_ref(l, m) * f;
    b22 = b_ref(m, m) * f;
    e = a11 / b11;
    c__ = ((a22 - e * b22) / b22 - a21 * b12 / (b11 * b22)) * .5;
    d__ = c__ * c__ + a21 * (a12 - e * b12) / (b11 * b22);
    if (d__ < 0.)
      goto L15;
    if (c__ >= 0.)
      e += c__ + (double)sqrt(d__);
    if (c__ < 0.)
      e += c__ - (double)sqrt(d__);
    a11 -= e * b11;
    a12 -= e * b12;
    a22 -= e * b22;
    flip = MIA_abs(a11) + MIA_abs(a12) >= MIA_abs(a21) + MIA_abs(a22);
    if (flip)
      ownBPL_am12r_c(&a12, &a11, &u1, &u2, &v1, &v2);
    if (! flip)
      ownBPL_am12r_c(&a22, &a21, &u1, &u2, &v1, &v2);
L8:
    if (u1 != 1.)
      goto L11;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      t = a_ref(i__, m) + u2 * a_ref(i__, l);
      a_ref(i__, m) = a_ref(i__, m) + v1 * t;
      a_ref(i__, l) = a_ref(i__, l) + v2 * t;
      t = b_ref(i__, m) + u2 * b_ref(i__, l);
      b_ref(i__, m) = b_ref(i__, m) + v1 * t;
      b_ref(i__, l) = b_ref(i__, l) + v2 * t;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      t = v_ref(i__, m) + u2 * v_ref(i__, l);
      v_ref(i__, m) = v_ref(i__, m) + v1 * t;
      v_ref(i__, l) = v_ref(i__, l) + v2 * t;
    }
L11:
    if (bn == 0.)
      goto L14;
    flip = an >= MIA_abs(e) * bn;
    if (flip)
      ownBPL_am12r_c(&b_ref(l, l), &b_ref(m, l), &u1, &u2, &v1, &v2);
L12:
    if (u1 != 1.)
      goto L14;
    if (! flip)
      ownBPL_am12r_c(&a_ref(l, l), &a_ref(m, l), &u1, &u2, &v1, &v2);
    i__1 = *n;
    for (j = l; j <= i__1; ++j)
    {
      t = a_ref(l, j) + u2 * a_ref(m, j);
      a_ref(l, j) = a_ref(l, j) + v1 * t;
      a_ref(m, j) = a_ref(m, j) + v2 * t;
      t = b_ref(l, j) + u2 * b_ref(m, j);
      b_ref(l, j) = b_ref(l, j) + v1 * t;
      b_ref(m, j) = b_ref(m, j) + v2 * t;
    }
L14:
    a_ref(m, l) = 0.;
    b_ref(m, l) = 0.;
    alfr[l] = a_ref(l, l);
    alfr[m] = a_ref(m, m);
    if (b_ref(l, l) < 0.)
      alfr[l] = -alfr[l];
    if (b_ref(m, m) < 0.)
      alfr[m] = -alfr[m];
    beta[l] = (r__1 = b_ref(l, l), MIA_abs(r__1));
    beta[m] = (r__1 = b_ref(m, m), MIA_abs(r__1));
    alfi[m] = 0.;
    alfi[l] = 0.;
    m += -2;
    goto L16;
L15:
    er = e + c__;
    ei = (double)sqrt(-d__);
    a11r = a11 - er * b11;
    a11i = ei * b11;
    a12r = a12 - er * b12;
    a12i = ei * b12;
    a21r = a21;
    a21i = 0.;
    a22r = a22 - er * b22;
    a22i = ei * b22;
    flip = MIA_abs(a11r) + MIA_abs(a11i) + MIA_abs(a12r) + MIA_abs(a12i) >= MIA_abs(a21r) +
            MIA_abs(a22r) + MIA_abs(a22i);
    if (flip)
    {
      r__1 = -a11r;
      r__2 = -a11i;
      ownBPL_am12c_c(&a12r, &a12i, &r__1, &r__2, &cz, &szr, &szi);
    }
    if (! flip)
    {
      r__1 = -a21r;
      r__2 = -a21i;
      ownBPL_am12c_c(&a22r, &a22i, &r__1, &r__2, &cz, &szr, &szi);
    }
    flip = an >= (MIA_abs(er) + MIA_abs(ei)) * bn;
    if (flip)
    {
      r__1 = cz * b11 + szr * b12;
      r__2 = szi * b12;
      r__3 = szr * b22;
      r__4 = szi * b22;
      ownBPL_am12c_c(&r__1, &r__2, &r__3, &r__4, &cq, &sqr, &sqi);
    }
    if (! flip)
    {
      r__1 = cz * a11 + szr * a12;
      r__2 = szi * a12;
      r__3 = cz * a21 + szr * a22;
      r__4 = szi * a22;
      ownBPL_am12c_c(&r__1, &r__2, &r__3, &r__4, &cq, &sqr, &sqi);
    }
    ssr = sqr * szr + sqi * szi;
    ssi = sqr * szi - sqi * szr;
    tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
    ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
    bdr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
    bdi = cq * szi * b12 + ssi * b22;
    r__ = (double)sqrt(bdr * bdr + bdi * bdi);
    beta[l] = bn * r__;
    alfr[l] = an * (tr * bdr + ti * bdi) / r__;
    alfi[l] = an * (tr * bdi - ti * bdr) / r__;
    tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
    ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21;
    bdr = ssr * b11 - sqr * cz * b12 + cq * cz * b22;
    bdi = -ssi * b11 - sqi * cz * b12;
    r__ = (double)sqrt(bdr * bdr + bdi * bdi);
    beta[m] = bn * r__;
    alfr[m] = an * (tr * bdr + ti * bdi) / r__;
    alfi[m] = an * (tr * bdi - ti * bdr) / r__;
    m += -2;
L16:
    if (m > 0)
      goto L3;
    m = *n;
L17:
    if (alfi[m] != 0.)
      goto L24;
    alfm = alfr[m];
    betm = beta[m];
    b_ref(m, m) = 1.;
    l = m - 1;
    if (l == 0)
      goto L23;
L18:
    l1 = l + 1;
    sl = 0.;
    i__1 = m;
    for (j = l1; j <= i__1; ++j)
      sl += (betm * a_ref(l, j) - alfm * b_ref(l, j)) * b_ref(j, m);
    if (l == 1)
      goto L20;
    if (betm * a_ref(l, l - 1) != 0.)
      goto L21;
L20:
    d__ = betm * a_ref(l, l) - alfm * b_ref(l, l);
    if (d__ == 0.)
      d__ = (epsa + epsb) * .5;
    b_ref(l, m) = -sl / d__;
    --l;
    goto L23;
L21:
    k = l - 1;
    sk = 0.;
    i__1 = m;
    for (j = l1; j <= i__1; ++j)
      sk += (betm * a_ref(k, j) - alfm * b_ref(k, j)) * b_ref(j, m);
    tkk = betm * a_ref(k, k) - alfm * b_ref(k, k);
    tkl = betm * a_ref(k, l) - alfm * b_ref(k, l);
    tlk = betm * a_ref(l, k);
    tll = betm * a_ref(l, l) - alfm * b_ref(l, l);
    d__ = tkk * tll - tkl * tlk;
    if (d__ == 0.)
      d__ = (epsa + epsb) * .5;
    b_ref(l, m) = (tlk * sk - tkk * sl) / d__;
    flip = MIA_abs(tkk) >= MIA_abs(tlk);
    if (flip)
      b_ref(k, m) = -(sk + tkl * b_ref(l, m)) / tkk;
    if (! flip)
      b_ref(k, m) = -(sl + tll * b_ref(l, m)) / tlk;
    l += -2;
L23:
    if (l > 0)
      goto L18;
    --m;
    goto L35;
L24:
    almr = alfr[m - 1];
    almi = alfi[m - 1];
    betm = beta[m - 1];
    mr = m - 1;
    mi = m;
    b_ref(m - 1, mr) = almi * b_ref(m, m) / (betm * a_ref(m, m - 1));
    b_ref(m - 1, mi) = (betm * a_ref(m, m) - almr * b_ref(m, m)) / (betm * 
            a_ref(m, m - 1));
    b_ref(m, mr) = 0.;
    b_ref(m, mi) = -1.;
    l = m - 2;
    if (l == 0)
      goto L34;
L25:
    l1 = l + 1;
    slr = 0.;
    sli = 0.;
    i__1 = m;
    for (j = l1; j <= i__1; ++j)
    {
      tr = betm * a_ref(l, j) - almr * b_ref(l, j);
      ti = -almi * b_ref(l, j);
      slr = slr + tr * b_ref(j, mr) - ti * b_ref(j, mi);
      sli = sli + tr * b_ref(j, mi) + ti * b_ref(j, mr);
    }
    if (l == 1)
      goto L27;
    if (betm * a_ref(l, l - 1) != 0.)
      goto L29;
L27:
    dr = betm * a_ref(l, l) - almr * b_ref(l, l);
    di = -almi * b_ref(l, l);
    iret = 1;
    xr = -slr;
    xi = -sli;
    yr = dr;
    yi = di;
    goto L47;
L28:
    b_ref(l, mr) = zr;
    b_ref(l, mi) = zi;
    --l;
    goto L34;
L29:
    k = l - 1;
    skr = 0.;
    ski = 0.;
    i__1 = m;
    for (j = l1; j <= i__1; ++j)
    {
      tr = betm * a_ref(k, j) - almr * b_ref(k, j);
      ti = -almi * b_ref(k, j);
      skr = skr + tr * b_ref(j, mr) - ti * b_ref(j, mi);
      ski = ski + tr * b_ref(j, mi) + ti * b_ref(j, mr);
    }
    tkkr = betm * a_ref(k, k) - almr * b_ref(k, k);
    tkki = -almi * b_ref(k, k);
    tklr = betm * a_ref(k, l) - almr * b_ref(k, l);
    tkli = -almi * b_ref(k, l);
    tlkr = betm * a_ref(l, k);
    tlki = 0.;
    tllr = betm * a_ref(l, l) - almr * b_ref(l, l);
    tlli = -almi * b_ref(l, l);
    dr = tkkr * tllr - tkki * tlli - tklr * tlkr;
    di = tkkr * tlli + tkki * tllr - tkli * tlkr;
    if (dr == 0. && di == 0.)
      dr = (epsa + epsb) * .5;
    iret = 2;
    xr = tlkr * skr - tkkr * slr + tkki * sli;
    xi = tlkr * ski - tkkr * sli - tkki * slr;
    yr = dr;
    yi = di;
    goto L47;
L31:
    b_ref(l, mr) = zr;
    b_ref(l, mi) = zi;
    flip = MIA_abs(tkkr) + MIA_abs(tkki) >= MIA_abs(tlkr);
    iret = 3;
    if (flip)
      goto L32;
    xr = -slr - tllr * b_ref(l, mr) + tlli * b_ref(l, mi);
    xi = -sli - tllr * b_ref(l, mi) - tlli * b_ref(l, mr);
    yr = tlkr;
    yi = tlki;
    goto L47;
L32:
    xr = -skr - tklr * b_ref(l, mr) + tkli * b_ref(l, mi);
    xi = -ski - tklr * b_ref(l, mi) - tkli * b_ref(l, mr);
    yr = tkkr;
    yi = tkki;
    goto L47;
L33:
    b_ref(k, mr) = zr;
    b_ref(k, mi) = zi;
    l += -2;
L34:
    if (l > 0)
      goto L25;
    m += -2;
L35:
    if (m > 0)
      goto L17;
    m = *n;
L36:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      s = 0.;
      i__2 = m;
      for (j = 1; j <= i__2; ++j)
        s += v_ref(i__, j) * b_ref(j, m);
      v_ref(i__, m) = s;
    }
    --m;
    if (m > 0)
      goto L36;
    m = *n;
L39:
    s = 0.;
    if (alfi[m] != 0.)
      goto L42;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      r__ = (r__1 = v_ref(i__, m), MIA_abs(r__1));
      if (r__ < s)
        continue;
      s = r__;
      d__ = v_ref(i__, m);
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
      v_ref(i__, m) = v_ref(i__, m) / d__;
    --m;
    goto L46;
L42:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      r__ = (r__1 = v_ref(i__, m - 1), MIA_abs(r__1)) + (r__2 = v_ref(i__, m), 
              MIA_abs(r__2));
      /* Computing 2nd power */
      r__1 = v_ref(i__, m - 1) / r__;
      /* Computing 2nd power */
      r__2 = v_ref(i__, m) / r__;
      r__ *= (double)sqrt(r__1 * r__1 + r__2 * r__2);
      if (r__ < s)
        continue;
      s = r__;
      dr = v_ref(i__, m - 1);
      di = v_ref(i__, m);
    }
    iret = 4;
    i__ = 0;
L44:
    ++i__;
    xr = v_ref(i__, m - 1);
    xi = v_ref(i__, m);
    yr = dr;
    yi = di;
    goto L47;
L45:
    v_ref(i__, m - 1) = zr;
    v_ref(i__, m) = zi;
    if (i__ < *n)
      goto L44;
    m += -2;
L46:
    if (m > 0)
      goto L39;
    goto L50;
L47:
    if (MIA_abs(yr) < MIA_abs(yi))
      goto L48;
    h__ = yi / yr;
    f = yr + h__ * yi;
    zr = (xr + h__ * xi) / f;
    zi = (xi - h__ * xr) / f;
    goto L49;
L48:
    h__ = yr / yi;
    f = yi + h__ * yr;
    zr = (h__ * xr + xi) / f;
    zi = (h__ * xi - xr) / f;
L49:
    switch (iret)
    {
      case 1:  goto L28;
      case 2:  goto L31;
      case 3:  goto L33;
      case 4:  goto L45;
    }
L50:
//error report to stdio - commented
    //if (*ierr != 0) {
  //      utag10_c(ierr, &c__41);
    //}
    return 0;
}
#undef v_ref
#undef b_ref
#undef a_ref

int BPL_Matrix_EigenGeneralized(
  double *a,
  double *b,
  double *v,
  double *alfa,
  double *beta,
  double *wk,
  int *n,
  int *ierr)
{
  int c__0 = 0;
  int a_dim1, a_offset, b_dim1, b_offset, i__1;
  int i__, j, n2, ja, ig, is, iz2, npi, igz;

    /* Parameter adjustments */
    --v;
    --alfa;
    --wk;
    --beta;
    b_dim1 = *n;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    /* Function Body */
    ownBPL_afg3r_c(&a[a_offset], &b[b_offset], &v[1], n, &c__0);
    ownBPL_agt0r_c(&a[a_offset], &b[b_offset], &v[1], &alfa[1], &alfa[*n + 1], &beta[1], n, ierr);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      npi = *n + i__;
      wk[i__] = alfa[npi];
    }
    ja = *n + *n;
    j = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      alfa[ja - 1] = alfa[j];
      alfa[ja] = wk[j];
      ja += -2;
      --j;
    }
    iz2 = *n + *n;
    n2 = *n + *n;
    j = *n;
L3:
    if (j < 1)
      return 0;
    if (alfa[j * 2] == 0.)
      goto L6;
    is = iz2 * (j - 1) + 1;
    ig = *n * (j - 2) + 1;
    igz = ig + *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      v[is] = v[ig];
      v[is + 1] = -v[igz];
      is += 2;
      ++ig;
      ++igz;
    }
    is = iz2 * (j - 2) + 1;
    ig = is + iz2;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      v[is] = v[ig];
      v[is + 1] = -v[ig + 1];
      is += 2;
      ig += 2;
    }
    j += -2;
    goto L3;
L6:
    is = iz2 * (j - 1) + n2;
    ig = *n * j;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      v[is - 1] = v[ig];
      v[is] = 0.;
      is += -2;
      --ig;
    }
    --j;
    goto L3;
}
