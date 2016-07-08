/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  19 March 2002                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Fast Walsh and Haar transforms                         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include "errorcodes.h"

#define DATATYPE double
#include "fht.h"
#undef DATATYPE

/*
// TODO: add Haar transform normalisation
// ref. code:
staticFunc inline int _HLP_CalcEHCoef(
  unsigned int n,    // number for which coef is estimated
  unsigned int p)    // log2(size)
{
    return 1<<( (p-MIA_IntLog2Lo(n>>p))+
                (p-MIA_IntLog2Lo(n&((1<<p)-1))) );
}
*/
