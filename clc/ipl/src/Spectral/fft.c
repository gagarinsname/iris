/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  February 9, 1989 - Original fft_nd code.               */
/*      Modified: April 17, 1989 - Added fft_1d and fft_2d.              */
/*      Modified: Febuary 15, 1991 - Added fft_3d.                       */
/*      Revision: 1.2.00                                                 */
/*      Purpose:  Fast Fourier transform                                 */
/*      Authors:                                                         */
/*        John Gauch                                                     */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <malloc.h>
#include "ipl.h"

#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define TWOPI 6.28318530717959

//double cos();
//double sin();

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its one-dimensional discrete      */
/*           transform if ISIGN=1 or replaces DATA by its inverse transform  */
/*           if ISIGN=-1.  DATA is a complex array of length NN which is     */
/*           input as a real array of length 2*NN.  NN must be an integer    */
/*           power of 2 or the routine will abort and return FALSE.          */
/*                                                                           */
/* Note:     Because this code was adapted from a FORTRAN library, the       */
/*           data array is 1-indexed.  In other words, the first element     */
/*           of the array is assumed to be in data[1].  Because C is zero    */
/*           indexed, the first element of the array is in data[0].  Hence,  */
/*           we must subtract 1 from the data address at the start of this   */
/*           routine so references to data[1] will really access data[0].    */
/*---------------------------------------------------------------------------*/
int IPL_FFT_1D /* fft_1d */(
  double *data,
  int nn,
  int isign)
{
  int n, mmax, m, j, istep, i;
  double wtemp, wr, wi, wpr, wpi, theta;
  double tempr, tempi;

    /* Fix indexing problems (see above) */
    data = data - 1;
    /* Error checking section */
/* what is it !?
    i = 1;
    j = nn;
    while (j > 1)
    {
      i = i << 1;
      j = j >> 1;
    }
    if (nn != i)
      return 0; */
    /* Bit reversal section */
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2)
    {
      if (j > i)
      {
        SWAP(data[j], data[i]);
        SWAP(data[j + 1], data[i + 1]);
      }
      m = n >> 1;
      while (m >= 2 && j > m)
      {
        j -= m;
        m = m >> 1;
      }
      j += m;
    }
    /* Danielson-Lanczos section */
    mmax = 2;
    while (n > mmax)
    {
      istep = 2 * mmax;
      theta = TWOPI / (isign * mmax);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (m = 1; m < mmax; m += 2)
      {
        for (i = m; i <= n; i += istep)
        {
          j = i + mmax;
          tempr = (double) (wr * data[j] - wi * data[j + 1]);
          tempi = (double) (wr * data[j + 1] + wi * data[j]);
          data[j] = data[i] - tempr;
          data[j + 1] = data[i + 1] - tempi;
          data[i] += tempr;
          data[i + 1] += tempi;
        }
        wtemp = wr;
        wr += wr * wpr - wi * wpi;
        wi += wi * wpr + wtemp * wpi;
      }
      mmax = istep;
    }
    /* Normalizing section */
    if (isign == 1)
    {
      n = nn << 1;
      for (i = 1; i <= n; i++)
      data[i] = data[i] / nn;
      /* printf("divide by %d\n", nn); */
    }
    return 1;
}

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its two-dimensional discrete      */
/*           transform if ISIGN=1 or replaces DATA by its inverse transform  */
/*           if ISIGN=-1.  DATA is a complex array with NN columns and MM    */
/*           rows.  NN and MM must both be integer powers of 2 or the        */
/*           routine will abort and return FALSE.                          */
/*---------------------------------------------------------------------------*/
int IPL_FFT_2D(
  double data[],
  int nn,
  int mm,
  int isign)
{
  int i, j, index1, index2;
  double *copy;

    // Error checking section
    i = 1;
    j = nn;
    while (j > 1)
    {
      i = i << 1;
      j = j >> 1;
    }
    if (nn != i)
      return 0;
    i = 1;
    j = mm;
    while (j > 1)
    {
      i = i << 1;
      j = j >> 1;
    }
    if (mm != i)
      return 0;
    // Transform by ROWS for forward transform
    if (isign == 1)
    {
      index1 = 0;
      for (i = 0; i < mm; i++)
      {
        IPL_FFT_1D(&(data[index1]), nn, isign);
        index1 += (nn << 1);
      }
    }
    // Allocate space for temporary array
    copy = (double *) malloc((unsigned) (sizeof(double) * mm * 2));
    if (copy == NULL)
      return 0;
    // Transform by COLUMNS
    for (j = 0; j < nn; j++)
    {
      // Copy pixels into temp array
      index1 = (j << 1);
      index2 = 0;
      for (i = 0; i < mm; i++)
      {
        copy[index2++] = data[index1];
        copy[index2++] = data[index1 + 1];
        index1 += (nn << 1);
      }
      // Perform transform
      IPL_FFT_1D(copy, mm, isign);
      // Copy pixels back into data array
      index1 = (j << 1);
      index2 = 0;
      for (i = 0; i < mm; i++)
      {
        data[index1] = copy[index2++];
        data[index1 + 1] = copy[index2++];
        index1 += (nn << 1);
      }
    }
    // Free temporary space
    free((char *) copy);
    // Transform by ROWS for inverse transform
    if (isign == -1)
    {
      index1 = 0;
      for (i = 0; i < mm; i++)
      {
        IPL_FFT_1D(&(data[index1]), nn, isign);
        index1 += (nn << 1);
      }
    }
    return 1;
}

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its three-dimensional discrete    */
/*           transform if ISIGN=1 or replaces DATA by its inverse transform  */
/*           if ISIGN=-1.  DATA is a complex array with NN columns and MM    */
/*           rows and OO slices.  NN, MM and OO must be integer powers of    */
/*           2 or the routine will abort and return FALSE.                 */
/*---------------------------------------------------------------------------*/
int IPL_FFT_3D( /* fft_3d */
  double data[],
  int nn,
  int mm,
  int oo,
  int isign)
{
   int i, j, index1, index2, nnmm;
   double *copy;

   /* Error checking section */
   i = 1;
   j = nn;
   while (j > 1)
   {
      i = i << 1;
      j = j >> 1;
   }
   if (nn != i)
      return 0;
   i = 1;
   j = mm;
   while (j > 1)
   {
      i = i << 1;
      j = j >> 1;
   }
   if (mm != i)
      return 0;
   i = 1;
   j = oo;
   while (j > 1)
   {
      i = i << 1;
      j = j >> 1;
   }
   if (oo != i)
      return 0;

   /* Transform by ROWS and COLUMNS for forward transform */
   nnmm = nn * mm;
   if (isign == 1)
   {
      index1 = 0;
      for (i = 0; i < oo; i++)
      {
	 IPL_FFT_2D(&(data[index1]), nn, mm, isign);
	 index1 += (nnmm << 1);
      }
   }

   /* Allocate space for temporary array */
   copy = (double *) malloc((unsigned) (sizeof(double) * 2 * oo));
   if (copy == NULL)
      return 0;

   /* Transform by SLICES */
   for (j = 0; j < nnmm; j++)
   {
      /* Copy pixels into temp array */
      index1 = (j << 1);
      index2 = 0;
      for (i = 0; i < oo; i++)
      {
	 copy[index2++] = data[index1];
	 copy[index2++] = data[index1 + 1];
	 index1 += (nnmm << 1);
      }

      /* Perform transform */
      IPL_FFT_1D(copy, oo, isign);

      /* Copy pixels back into data array */
      index1 = (j << 1);
      index2 = 0;
      for (i = 0; i < oo; i++)
      {
	 data[index1] = copy[index2++];
	 data[index1 + 1] = copy[index2++];
	 index1 += (nnmm << 1);
      }
   }

   /* Free temporary space */
   free((char *) copy);

   /* Transform by COLUMNS and ROWS for inverse transform */
   if (isign == -1)
   {
      index1 = 0;
      for (i = 0; i < oo; i++)
      {
	 IPL_FFT_2D(&(data[index1]), nn, mm, isign);
	 index1 += (nnmm << 1);
      }
   }

   return 1;
}

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its DIMC-dimensional discrete     */
/*           Fourier transform, if SIGN is input as 1.  DIMV[0..DIMC-1]      */
/*           is an integer array containing the lengths of each dimension    */
/*           (number of complex values) , which MUST be a power of 2.        */
/*           DATA is a real array of length twice the product of these       */
/*           lengths, in which the data is stored in a multidimensional      */
/*           complex array: real and imaginary parts of this array are in    */
/*           consecutive locations, and the rightmost index of the array     */
/*           increases most rapidly as one proceeds along DATA.  For a two-  */
/*           dimensional array, this is equivalent to storing the array      */
/*           by rows (the C norm).  If SIGN is -1, DATA is replaced by its   */
/*           inverse transform times the product of the lengths of all       */
/*           dimensions (in other words, divide by the number of pixels      */
/*           in the image to get the inverse transform).                     */
/*                                                                           */
/* Note:     Because this code was adapted from a FORTRAN library, the       */
/*           data array is 1-indexed.  In other words, the first element     */
/*           of the array is assumed to be in data[1].  Because C is zero    */
/*           indexed, the first element of the array is in data[0].  Hence,  */
/*           the address of data[-1] must be passed to this routine to       */
/*           ensure that references to data[1] will really access data[0].   */
/*---------------------------------------------------------------------------*/
int IPL_FFT_ND( /* fft_nd */
  double data[],
  int nn[],
  int ndim,
  int isign)
{
   /* Local variables */
   int i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
   int ibit, idim, k1, k2, n, nprev, nrem, ntot;
   double tempi, tempr;
   double theta, wi, wpi, wpr, wr, wtemp;

   /* Compute total number of complex values */
   ntot = 1;
   for (idim = 0; idim < ndim; idim++)
      ntot *= nn[idim];

   /* Main loop over the dimensions of the image */
   nprev = 1;
   for (idim = ndim - 1; idim >= 0; idim--)
   {
      n = nn[idim];
      nrem = ntot / (n * nprev);
      ip1 = nprev << 1;
      ip2 = ip1 * n;
      ip3 = ip2 * nrem;
      i2rev = 1;

      /* This is the bit reversal section */
      for (i2 = 1; i2 <= ip2; i2 += ip1)
      {
	 if (i2 < i2rev)
	 {
	    for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2)
	    {
	       for (i3 = i1; i3 <= ip3; i3 += ip2)
	       {
		  i3rev = i2rev + i3 - i2;
		  SWAP(data[i3], data[i3rev]);
		  SWAP(data[i3 + 1], data[i3rev + 1]);
	       }
	    }
	 }
	 ibit = ip2 >> 1;
	 while ((ibit >= ip1) && (i2rev > ibit))
	 {
	    i2rev -= ibit;
	    ibit >>= 1;
	 }
	 i2rev += ibit;
      }

      /* Here is the Danielson-Lanczos section */
      ifp1 = ip1;
      while (ifp1 < ip2)
      {
	 ifp2 = ifp1 << 1;
	 theta = isign * TWOPI / (ifp2 / ip1);
	 wtemp = sin(0.5 * theta);
	 wpr = -2.0 * wtemp * wtemp;
	 wpi = sin(theta);
	 wr = 1.0;
	 wi = 0.0;
	 for (i3 = 1; i3 <= ifp1; i3 += ip1);
	 {
	    for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2)
	    {
	       for (i2 = i1; i2 <= ip3; i2 += ip2)
	       {
		  k1 = i2;
		  k2 = k1 + ifp1;
		  tempr = (double) (wr * data[k2] - wi * data[k2 + 1]);
		  tempi = (double) (wr * data[k2 + 1] - wi * data[k2]);
		  data[k2] = data[k1] - tempr;
		  data[k2 + 1] = data[k1 + 1] - tempi;
		  data[k1] = tempr;
		  data[k1 + 1] = tempi;
	       }
	    }
	    wtemp = wr;
	    wr += wr * wpr - wi * wpi;
	    wi += wi * wpr + wtemp * wpi;
	 }
	 ifp1 = ifp2;
      }

      /* Get ready for next image dimension */
      nprev *= n;
   }

   return 1;
}

/*
void im_free2D(char **Data)
{
   free(Data[0]);
   free(Data);
}
*/
