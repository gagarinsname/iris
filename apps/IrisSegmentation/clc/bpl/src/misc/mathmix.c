/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  6 March 2004                                           */
/*      Modified: 22 July 2004 (v2)                                      */
/*      Modified: 27 August 2004 - MIA clarified it                      */
/*      Modified: 1 September 2004 - MIA changed atan to atan2           */
/*      Revision: 2.2.00                                                 */
/*      Purpose:  math funcrions                                         */
/*        sqrt                                                           */
/*        atan2                                                          */
/*        cos (sin is derived)                                           */
/*        log2                                                           */
/*        pow2 (exp, pow are derived)                                    */
/*      Authors:                                                         */
/*        Gevorg Hmayakyan (HGA)                                         */
/*        Konstantin Gankin (GKA)                                        */
/*        Ivan Matveev (MIA)                                             */
/*        Vladimir Novik (NVP)                                           */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#include <math.h>
#include "stddefs.h"
#include "bpl.h"

#define PI_FIX20    3294199   // representation of PI in fix20 arithmetic
#define PIx2_FIX20  6588397   // representation of 2*PI in fix20 arithmetic
#define PId2_FIX20  1647099   // representation of PI/2 in fix20 arithmetic

// include tables
#define __INCLUDE_TABLES_C__
#include "tables.h"

// square root
#define INNER_MBGSQRT(s)                      \
  temp = (g << (s)) + (1 << ((s) * 2 - 2));   \
  if (val >= temp)                            \
  {                                           \
    g += 1 << ((s)-1);                        \
    val -= temp;                              \
  }

unsigned int MIA_INT_sqrt(
  unsigned int val)
{
  unsigned int temp, g=0;

    if (val >= 0x40000000)
    {
      g = 0x8000; 
      val -= 0x40000000;
    }
    INNER_MBGSQRT (15)
    INNER_MBGSQRT (14)
    INNER_MBGSQRT (13)
    INNER_MBGSQRT (12)
    INNER_MBGSQRT (11)
    INNER_MBGSQRT (10)
    INNER_MBGSQRT ( 9)
    INNER_MBGSQRT ( 8)
    INNER_MBGSQRT ( 7)
    INNER_MBGSQRT ( 6)
    INNER_MBGSQRT ( 5)
    INNER_MBGSQRT ( 4)
    INNER_MBGSQRT ( 3)
    INNER_MBGSQRT ( 2)
    if (val >= g+g+1)
      g++;
    return g;
}

#undef INNER_MBGSQRT

// cosine with fix20 integer arithmetic
int MIA_FIX20_cos(
  int x)
{
  int neg,xn,tMin;

    // cos is an even function
    if (x<0)
      x = -x;
    // cos is cyclic
    if (x>=(int)(((float)PI)*(1<<21)))
      x %= (int)(((float)PI)*(1<<21));
    // cos(x+pi) = -cos(x)
    if (x>=(int)(((float)PI)*(1<<20)))
    {
      x -= (int)(((float)PI)*(1<<20));
      neg = -1;
    }
    else
      neg = 1;
    // cos(x+pi/2) = -cos(pi/2-x)
    if (x>=(int)(((float)PI)*(1<<19)))
    {
      neg = -neg;
      x = (int)(((float)PI)*(1<<20))-x;
    }
    // table
    xn = x>>10;
    tMin = table_cos[xn];
    tMin += ((table_cos[xn+1]-tMin)*(x-(xn<<10)))>>10;
    return (neg>0)?tMin:(-tMin);
}

#define INT_ORDER 32
#define LOG2_TAB_ORDER 10
#define POW2_TAB_ORDER 10
#define FIX_ORDER 20
#define FIX_ORDER_TBLPOW 31

// log2 of unsigned fix20 value to signed fix20 value
int MIA_FIX20_log2_of_FIX20(
  unsigned int x)
{
  int order,mask;

    if (!x)
      return (int)(0x80000000); // infinity
    // search for first unit from big end
    order = (INT_ORDER-1);
    mask = 1<<(INT_ORDER-1);
    while (!(x&mask))
    {
      order--;
      mask >>= 1;
    }
    // Shift the number to fit it's first unit to bit 
    // position INT_ORDER+LOG2_TAB_ORDER-FIX_ORDER. (22 here)
    // LOG2_TAB_ORDER bits (21-12) will give table index, 
    // and other INT_ORDER-FIX_ORDER bits (11-0) will be used for interpolation. 
    // Unit bit in position (22) is useless from now on. 
    if (order>INT_ORDER+LOG2_TAB_ORDER-FIX_ORDER)
      x >>= order-(LOG2_TAB_ORDER+INT_ORDER-FIX_ORDER);
    else
      x <<= (LOG2_TAB_ORDER+INT_ORDER-FIX_ORDER)-order;
    x &= (1<<(LOG2_TAB_ORDER+INT_ORDER-FIX_ORDER))-1;
    // calculate base index
    mask = x>>(INT_ORDER-FIX_ORDER);
    // 
    return ((order-FIX_ORDER)<<FIX_ORDER)+            // integer part
           table_log2[mask]+                          // table value
           (((table_log2[mask+1]-table_log2[mask])*
             (x&((1<<(INT_ORDER-FIX_ORDER))-1)))>>
             (INT_ORDER-FIX_ORDER));                  // interpolation
}

// log2 of unsigned integer to unsigned fix20
int MIA_FIX20_log2_of_int(
  unsigned int x)
{
  int order,mask;

    if (!x)
      return (int)(0x80000000); // infinity
    // search for first unit from big end
    order = (INT_ORDER-1);
    mask = 1<<(INT_ORDER-1);
    while (!(x&mask))
    {
      order--;
      mask >>= 1;
    }
    // align highest '1' to 33-th bit (omit it in fact)
    x <<= INT_ORDER-order;
    // get elder 10 bits as index 
    mask = x>>(INT_ORDER-LOG2_TAB_ORDER);
    // align subsequent bits to start from fix_order
    x = (x>>(INT_ORDER-LOG2_TAB_ORDER-FIX_ORDER))&((1<<FIX_ORDER)-1);
    // 
    return (order<<FIX_ORDER)+            // integer part
      table_log2[mask]+                          // table value
      (((table_log2[mask+1]-table_log2[mask])*x)>>FIX_ORDER); // interpolation
}

// pow2 of signed fix20 to unsigned fix20
unsigned int MIA_FIX20_pow2_of_FIX20(
  int x)
{
  int shift,index,interpol;

    shift = x>>FIX_ORDER;
    if (shift>=INT_ORDER-FIX_ORDER)
      return 0xffffffff; // infinity
    index = (x>>(FIX_ORDER-LOG2_TAB_ORDER))&((1<<LOG2_TAB_ORDER)-1);
    interpol = x&((1<<(FIX_ORDER-LOG2_TAB_ORDER))-1);
    // 
    if (shift>=0)
      return (table_pow2[index]+                          // table value
        (((table_pow2[index+1]-table_pow2[index])*interpol)>>10)) // interpolation
        <<shift;    // integer
    else
      return (table_pow2[index]+                          // table value
        (((table_pow2[index+1]-table_pow2[index])*interpol)>>10)) // interpolation
        >>shift;    // integer
}

// atan2 function with integer arguments
// trigonometric circle zones are denoted here: 
//    |   |
//  6 | 2 | 5
// ---+---+---
//  3 | 0 | 1
// ---+---+---
//  7 | 4 | 8
//    |   |
// 0 - center, 1,2,3,4 - axes, 5,6,7,8 - non-zero quadrants

int MIA_FIX20_atan2(
  int dy,
  int dx)
{
  int neg,subpi,inv,x,xn,tMin;

    if (dy==0)
      // zones 1,0,3
      return (dx>=0)?0:PI_FIX20;  // 0 | PI
    if (dx==0)
      // zones 2,4
      return  (dy>0)?PId2_FIX20:    //  PI/2
                    -PId2_FIX20;    // -PI/2
    // zones 5,6,7,8
    // map 7->6, 8->5
    if (dy<0)
    {
      neg = 1;
      dy = -dy;
    }
    else
      neg = 0;
    // zones 5,6
    // map 6->5
    if (dx<0)
    {
      subpi = 1;
      dx = -dx;
    }
    else
      subpi = 0;
    // zone 5
    // map higher half of 5 to lower half
    if (dy>=dx)
    {
      inv = 1;
      x = ((1<<20)*dx)/dy;
    }
    else
    {
      inv = 0;
      x = ((1<<20)*dy)/dx;
    }
    // table with interpolation
    xn = x>>10;
    tMin = table_atan[xn];
    tMin += ((table_atan[xn+1]-tMin)*(x-(xn<<10)))>>10;
    // count inversion
    if (inv)
      tMin = PId2_FIX20-tMin;
    // count pi
    if (subpi)
      tMin = PI_FIX20-tMin;
    // count sign
    return (neg)?(-tMin):tMin;
}

// return log2 of a number rounded down
unsigned int MIA_IntLog2Lo(
  unsigned int val)
{
  staticConst unsigned char ByteLog2LUT[256] = 
  { // starting value is never touched exept if val=0
    0,0,1,1,2,2,2,2, 3,3,3,3,3,3,3,3,   // 00-0f
    4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,   // 10-1f
    5,5,5,5,5,5,5,5, 5,5,5,5,5,5,5,5,   // 20-2f
    5,5,5,5,5,5,5,5, 5,5,5,5,5,5,5,5,   // 30-3f
    6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,   // 40-4f
    6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,   // 50-5f
    6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,   // 60-6f
    6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,   // 70-7f
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // 80-8f
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // 90-9f
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // a0-af
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // b0-bf
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // c0-cf
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // d0-df
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,   // e0-ef
    7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7    // f0-ff
  };
    return 
      (!(val&0xffff0000))?
        ( (!(val&0x0000ff00))? (ByteLog2LUT[val    ]   ) : (ByteLog2LUT[val>> 8]+ 8) ):
        ( (!(val&0xff000000))? (ByteLog2LUT[val>>16]+16) : (ByteLog2LUT[val>>24]+24) );
}

// return log2 of a number rounded up
unsigned int MIA_IntLog2Hi(
  unsigned int val)
{
    if (val)
      return MIA_IntLog2Lo(val-1)+1;
    return 0;
}

// search local minimum of a function of one variable
int BPL_OPT_LocateMinimum(
  BPL_OPT_LocateMinimum_Callback funct,
  const void* context,
  double *a,
  double *b,
  double *c__,
  double *eps,
  double *xmin,
  double *fmin)
{
  /* Initialized data */
  const double r__ = .61803399f;
  const double one = 1.f;
  /* System generated locals */
  double r__1, r__2;
  /* Local variables */
  double g, f0, f1, f2, f3, x0, x1, x2, x3;

    g = one - r__;
    x0 = *a;
    x3 = *c__;
    if ((r__1 = *c__ - *b, fabs(r__1)) <= (r__2 = *b - *a, fabs(r__2)))
      goto L10;
    x1 = *b;
    x2 = *b + g * (*c__ - *b);
    goto L20;
L10:
    x2 = *b;
    x1 = *b - g * (*b - *a);
L20:
    f1 = (*funct)(x1,context);
    f2 = (*funct)(x2,context);
L30:
    if ((r__1 = x3 - x0, fabs(r__1)) <= *eps * (fabs(x1) + fabs(x2)))
      goto L50;
    if (f2 > f1)
      goto L40;
    x0 = x1;
    x1 = x2;
    x2 = r__ * x1 + g * x3;
    f0 = f1;
    f1 = f2;
    f2 = (*funct)(x2,context);
    goto L30;
L40:
    x3 = x2;
    x2 = x1;
    x1 = r__ * x2 + g * x0;
    f3 = f2;
    f2 = f1;
    f1 = (*funct)(x1,context);
    goto L30;
L50:
    if (f1 >= f2)
      goto L60;
    *fmin = f1;
    *xmin = x1;
    goto L70;
L60:
    *fmin = f2;
    *xmin = x2;
L70:
    return 0;
}
