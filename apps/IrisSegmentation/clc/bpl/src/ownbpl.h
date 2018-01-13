/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  8 February 2010                                        */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  basic procedures library                               */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __OWNBPL_H__
#define __OWNBPL_H__

#include "bpl.h"

/* linear algebra errors ----------------------------------------------*/
/* invalid object size */
#define MIA_LIN_INVALID_SIZE ENHANCED_ERROR_CODE(0x101)
/* Jacobi algorithm could not converge in a given number of rotations */
#define MIA_LIN_JACOBI_FAILED_CONVERGE ENHANCED_ERROR_CODE(0x102)
/* matrix is not square */
#define MIA_LIN_MATRIX_NOT_SQUARE ENHANCED_ERROR_CODE(0x103)
/* matrix is degenerate */
#define MIA_LIN_MATRIX_DEGENERATE ENHANCED_ERROR_CODE(0x104)

/* PCA & LDA errors ---------------------------------------------------------*/
/* not all declared samples are present */
#define MIA_PCA_NOSAMPLES ENHANCED_ERROR_CODE(0x141)
/* mean is not calculated */
#define MIA_PCA_NOMEAN ENHANCED_ERROR_CODE(0x142)
/* number of necessary vectors required is more than vector where rounding error is yet small */
#define MIA_PCA_INSUFFICIENT_RELIABLE_EIGENS_ROUNDING ENHANCED_ERROR_CODE(0x143)
/* eigen dicrepancy too big to be used in further calculations */
#define MIA_PCA_TOO_BIG_DISCREPANCY ENHANCED_ERROR_CODE(0x144)
/* L[0]/L[n]>M*M condtionality restriction does not allow sufficient eigens */
#define MIA_PCA_INSUFFICIENT_RELIABLE_EIGENS_CONDITIONALITY ENHANCED_ERROR_CODE(0x145)
/* LDA should have at least two classes */
#define MIA_LDA_TOO_FEW_CLASSES ENHANCED_ERROR_CODE(0x161)
/* class is not adequately populated */
#define MIA_LDA_TOO_SMALL_CLASS ENHANCED_ERROR_CODE(0x162)
/* class indexes should go with unit increment, but they do not */
#define MIA_LDA_UNCONTIGUIOS_INDEXES ENHANCED_ERROR_CODE(0x163)

/* Optimisation errors ----------------------------------------------*/
/* optimisation have not converged in a given number of iterations */
#define MIA_OPT_NOT_CONVERGED ENHANCED_ERROR_CODE(0x181)
/* not a convex function - unable to find global minimum */
#define MIA_OPT_NOT_CONVEX ENHANCED_ERROR_CODE(0x182)
/* function has one sign on segment edges - bissection rejected */
#define MIA_OPT_ONE_SIGN ENHANCED_ERROR_CODE(0x183)

#endif // __OWNBPL_H__
