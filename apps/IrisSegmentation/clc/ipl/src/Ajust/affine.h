/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  13 June 2002                                           */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Image affine mapping by means of least squares method  */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include "bpl.h"

#ifndef __CIMAGEAFFINEMAPPING_H__
#define __CIMAGEAFFINEMAPPING_H__

typedef struct
{
  double S_xb_xb;
  double S_xb_yb;
  double S_yb_yb;
  double S_xb;
  double S_yb;
  double S_1;
  double S_xt_xb;
  double S_xt_yb;
  double S_xt;
  double S_yt_xb;
  double S_yt_yb;
  double S_yt;
  SLinearSystem cls_x;
  SLinearSystem cls_y;
  unsigned int isSolved;
} SImageAffineMapping;

void SImageAffineMapping_SImageAffineMapping(
  SImageAffineMapping* pThis);

void SImageAffineMapping_UpdateMapping(
  SImageAffineMapping* pThis,
  double xb,
  double yb,
  double xt,
  double yt);

void SImageAffineMapping_UnupdateMapping(
  SImageAffineMapping* pThis,
  double xb,
  double yb,
  double xt,
  double yt);

MIA_RESULT_CODE SImageAffineMapping_Map(
  SImageAffineMapping* pThis,
  double* x,
  double* y);

MIA_RESULT_CODE SImageAffineMapping_GetCoeffs(
  SImageAffineMapping* pThis,
  double *Coeffs);

#endif // __CIMAGEAFFINEMAPPING_H__
