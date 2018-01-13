/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  8 February 2010                                        */
/*      Modified: 19 July 2010 - revised for CodCon                      */
/*      Revision: 2.0.00                                                 */
/*      Purpose:  basic procedures library                               */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __bpl_h__
#define __bpl_h__

// common headers
#include <math.h>
#include "stddefs.h"
#include "errorcodes.h"

EXTERNC MIA_RESULT_CODE BPL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize);      // IN:  string length

// uint32 random number generation
#define ES_rand(n)  (1664525L*(n)+1013904223L)

// bound the window to frame
#define FRAMEBOUND(x,w,xs) \
  if (x<0) \
  { \
    w += x; \
    x = 0;\
  } \
  if (x+w>xs) \
    w = xs-x;

//== heap sort algorithm ================================================
#define IMPLEMENT_HEAPSORT( func_name, T, less_than )                   \
void func_name( T* array, int length )                                  \
{                                                                       \
    const int bubble_level = 8;                                         \
                                                                        \
    struct                                                              \
    {                                                                   \
        int lb, ub;                                                     \
    }                                                                   \
    stack[48];                                                          \
                                                                        \
    int sp = 0;                                                         \
                                                                        \
    T   temp;                                                           \
    T   lb_val;                                                         \
                                                                        \
    stack[0].lb = 0;                                                    \
    stack[0].ub = length - 1;                                           \
                                                                        \
    while( sp >= 0 )                                                    \
    {                                                                   \
        int lb = stack[sp].lb;                                          \
        int ub = stack[sp--].ub;                                        \
                                                                        \
        for(;;)                                                         \
        {                                                               \
            int diff = ub - lb;                                         \
            if( diff < bubble_level )                                   \
            {                                                           \
                int i, j;                                               \
                T* arr = array + lb;                                    \
                                                                        \
                for( i = diff; i > 0; i-- )                             \
                {                                                       \
                    int f = 0;                                          \
                    for( j = 0; j < i; j++ )                            \
                        if( less_than( arr[j+1], arr[j] ))              \
                        {                                               \
                            temp = arr[j];                              \
                            arr[j] = arr[j+1];                          \
                            arr[j+1] = temp;                            \
                            f = 1;                                      \
                        }                                               \
                    if( !f ) break;                                     \
                }                                                       \
                break;                                                  \
            }                                                           \
            else                                                        \
            {                                                           \
                /* select pivot and exchange with 1st element */        \
                int  m = lb + (diff >> 1);                              \
                int  i = lb + 1, j = ub;                                \
                                                                        \
                lb_val = array[m];                                      \
                                                                        \
                array[m]  = array[lb];                                  \
                array[lb] = lb_val;                                     \
                                                                        \
                /* partition into two segments */                       \
                for(;;)                                                 \
                {                                                       \
                    for( ;i < j && less_than(array[i], lb_val); i++ );  \
                    for( ;j >= i && less_than(lb_val, array[j]); j-- ); \
                                                                        \
                    if( i >= j ) break;                                 \
                    temp = array[i];                                    \
                    array[i++] = array[j];                              \
                    array[j--] = temp;                                  \
                }                                                       \
                                                                        \
                /* pivot belongs in A[j] */                             \
                array[lb] = array[j];                                   \
                array[j]  = lb_val;                                     \
                                                                        \
                /* keep processing smallest segment,and stack largest*/ \
                if( j - lb <= ub - j )                                  \
                {                                                       \
                    if( j + 1 < ub )                                    \
                    {                                                   \
                        stack[++sp].lb   = j + 1;                       \
                        stack[sp].ub = ub;                              \
                    }                                                   \
                    ub = j - 1;                                         \
                }                                                       \
                else                                                    \
                {                                                       \
                    if( j - 1 > lb)                                     \
                    {                                                   \
                        stack[++sp].lb = lb;                            \
                        stack[sp].ub = j - 1;                           \
                    }                                                   \
                    lb = j + 1;                                         \
                }                                                       \
            }                                                           \
        }                                                               \
    }                                                                   \
}
#define DECLARE_HEAPSORT(func_name,T) void func_name(T*,int)

MIA_RESULT_CODE BPL_InterpolateSequence_SimpleIntDM(
  int* pOutVal,       // OUT: buffer for output sequence
  int nOutLen,        // IN:  output sequence length
  const int* pInVal,  // IN:  input sequence values
  int nInLen,         // IN:  input sequence length
  int nBegDup);       // IN:  number of beginning duplicates

EXTERNC int32 BPL_ReadInt32Bigendian(
  const void* ptr);

EXTERNC void BPL_InvertEndianness_uint32(
  void* ptr,
  int num);

//== interpolating points by elliptic figures ===================================
// calculate inertia-equivalent ellipse parameters of a cluster of points
EXTERNC int BPL_ELL_EllipseParamsOfPointGroup(
  double* ell_x,      // OUT: center position
  double* ell_y,      //
  double* ell_a,      // OUT: axis
  double* ell_b,      //
  double* ell_t,      // OUT: tilt
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M);             // number of points

// interpolate group of point by circle using MSD method
EXTERNC MIA_RESULT_CODE BPL_ELL_InterpolateByCircle(
  double* cir_x,      // OUT: center position
  double* cir_y,      //
  double* cir_r,      // OUT: axis
  double* pDiscr,     // OUT: discrepancy (set NULL to skip its calculation)
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M);             // number of points

// interpolate group of point by ellipse using MSD method
EXTERNC MIA_RESULT_CODE BPL_ELL_InterpolateByEllipse(
  double* ell_x,      // OUT: center position
  double* ell_y,      //
  double* ell_a,      // OUT: axes
  double* ell_b,      //
  double* ell_t,      // OUT: tilt
  double* pDiscr,     // OUT: discrepancy
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M);             // number of points

// special version for ES: no doubles, no int64
EXTERNC MIA_RESULT_CODE BPL_ELL_InterpolateByEllipseFloat(
  float* ell_x,       // OUT: center position
  float* ell_y,       //
  float* ell_a,       // OUT: axes
  float* ell_b,       //
  float* ell_t,       // OUT: tilt
  float* pDiscr,      // OUT: discrepancy
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M);             // number of points

// calculate inertia-equivalent ellipse parameters from simple moments
EXTERNC MIA_RESULT_CODE BPL_ELL_GetEquivalentEllipse(
  double* ellpar,         // OUT: xc,yc,a,b,tilt
  const double* moments); // IN:  M,Mx,My,Mxx,Mxy,Myy

//== fixed-point calculations ===========================================
EXTERNC int MIA_FIX20_atan2(int dy,int dx);
EXTERNC int MIA_FIX20_cos(int x);
EXTERNC int MIA_FIX20_log2_of_FIX20(unsigned int x);
EXTERNC int MIA_FIX20_log2_of_int(unsigned int x);
EXTERNC unsigned int MIA_FIX20_pow2_of_FIX20(int x);
#define MIA_FIX20_sin(x) MIA_FIX20_cos((x)-(int)(((float)PI)*(1<<19)))

//== CRC claculations ===================================================
EXTERNC unsigned int MIA_CalculateCRC32_M1(
  const void* pvBuf,
  int nLen,
  unsigned int nKey);

EXTERNC unsigned int MIA_CalculateCRC32_M2(
  unsigned int crc,
  const unsigned char* buf,
  unsigned int len);

EXTERNC unsigned int MIA_INT_sqrt(
  unsigned int val);

//== number of bits =====================================================
// return log2 of a number rounded down
EXTERNC unsigned int MIA_IntLog2Lo(
  unsigned int val);

// return log2 of a number rounded up
EXTERNC unsigned int MIA_IntLog2Hi(
  unsigned int val);

// return log2 of a number rounded down
#define MIA_IntLog2Lo_Def(logval,val) \
{ \
  int n; \
  uint64 v = (uint64)(val); \
  for (n=0;v;v>>=1,n++); \
  logval = n; \
}

// return log2 of a number rounded up
#define MIA_IntLog2Hi_Def(logval,val) \
{ \
  int n; \
  uint64 v = ((uint64)(val))-1; \
  for (n=0;v;v>>=1,n++); \
  logval = n+1; \
}

//== matrix operations ==================================================
// in-place transpose of square matrix
EXTERNC void BPL_Matrix_Transpose(
  double* matr,
  int size);

// pack symmetric matrix to 'triangular form'
EXTERNC void BPL_Matrix_Pack(
  double *b,        // OUT: packed buffer ('triangular form')
  const double *a,  // IN:  source matrix 
  int n);           // IN:  size of matrix (number of rows=number of columns)

// unpack 'triangular form' to symmetric matrix
EXTERNC void BPL_Matrix_Unpack(
  double *b,        // OUT: symmetric matrix
  const double *a,  // IN:  'triangular form'
  int n);           // IN:  matrix size

// float version of Unpack
EXTERNC void BPL_Matrix_UnpackFloat(
  float *b,         // OUT: symmetric matrix
  const float *a,   // IN:  'triangular form'
  int n);           // IN:  matrix size

// invert packed matrix
EXTERNC int BPL_Matrix_InvertPacked(
  double *ainv,   // OUT: inverted matrix in packed form
  double *a,      // IN:  source matrix in packed form (destroyed)
  int n);         // IN:  size of matrix

// Solve linear system Ax=b
// Source: BChA NIVC MGU, refer to as: ash0d_c
EXTERNC int BPL_Matrix_SymmetricLinearSolve(
  double* a,    // IN:  matrix A (size N*N)
  double* b,    // IN:  vector b
  double* x,    // OUT: vector x
  double* s,    // TMP: buffer for N values
  int n,        // IN:  order of system
  int isFirst); // IN:  is function called first time with this matrix

// float version of SymmetricLinearSolve
EXTERNC int BPL_Matrix_SymmetricLinearSolveFloat(
  float* a,     // IN:  matrix A (size N*N)
  float* b,     // IN:  vector b
  float* x,     // OUT: vector x
  float* s,     // TMP: buffer for N values
  int n,        // IN:  order of system
  int isFirst); // IN:  is function called first time with this matrix

EXTERNC int BPL_Matrix_SVD(
  const double *a,  // вещественный двумерны?массив размер?nm на n, ?первых m стpоках которого задает? исходн? матриц?
  double *w,  // вещественный вект? длин?n значений w(k) = dk вычисленны?сингулярных чисе?матриц?a
  double *u,  // вещественный двумерны?массив размер?nm на n, ?котоpом пр?matu = .true. содержат? вычисленны?первые n столбцов матриц?u; пр?matu = .false. массив u использует? ка?рабочи?
  int matu,   // призна?необходимост?вычислен? матриц?u: пр?matu = .true. матриц?u вычисляет?, пр?matu = .false. - не?
  double *v,  // вещественный двумерны?массив размер?nm на n, ?котоpом пр?matv = .true. содержит? вычисленная матриц?v; пр?matv = .false. массив v не использует?
  int matv,   // призна?необходимост?вычислен? матриц?v: пр?matv = .true. матриц?v вычисляет?, пр?matv = .false. - не?
  int m,      // числ?стpок матриц a ?u
  int n,      // числ?столбцов матриц a ?u ?по?до?матриц?v
  double *r); // IN: temporary memory, size = n*sizeof(double)

// solution of common eigen value problem Ax=lx by QL algorithm
// where A is real symmetric
// see www.srcc.msu.su/num_anal/lib_na/cat/ae_htm_c/aeh1r_c.htm
EXTERNC MIA_RESULT_CODE BPL_Matrix_Eigen_Symm(
  double *a,        // IN:  source matrix / OUT: eigenvector matrix
  double *ev,       // OUT: eigen values
  int n,            // IN:  size of matrix
  void* extbuftmp,  // IN:  temporary buffer
  int nTMP);        // IN:  size of temporary buffer

// solution of generalized eigen problem Ax=lBx
// where A,B - real symmetric, B is positively defined
// see www.srcc.msu.su/num_anal/lib_na/cat/ag_htm_c/agh0r_c.htm
EXTERNC MIA_RESULT_CODE BPL_Matrix_EigenGeneralized_Symm(
  double *a,        // IN:  left-size matrix (A)
  double *b,        // IN:  right-side matrix (B)
  double *v,        // OUT: eigen vectors
  double *ev,       // OUT: eigen values
  int n,            // IN:  dimensionality
  void* extbuftmp,  // IN:  temporary buffer
  int nTMP);        // IN:  size of temporary buffer

// solution of generalized eigen problem Ax=lBx
// no limitations put on A,B
// see www.srcc.msu.su/num_anal/lib_na/cat/ag_htm_c/agg0r_c.htm
EXTERNC int BPL_Matrix_EigenGeneralized(
  double *a,
  double *b,
  double *v,
  double *alfa,
  double *beta,
  double *wk,
  int *n,
  int *ierr);

// calculation by Jacobi rotations, see Numerical recipies
// preferrable for small and diagonally expressed matrices
// for large matrices QL algorithm is better
EXTERNC MIA_RESULT_CODE BPL_Matrix_EigenJacobiRotations(
  double* Matr, // IN: symmetric matrix which eigens are calculated
                //     elements above the diagonal are destroyed
  int size,     // IN: size of this matrix
  double* vals, // OUT: eigevalues of A
  double* Vect, // OUT: eigenvectors of A (in columns)
  int* sweeps,  // number of Jacobi rotation sweeps
                // IN: thresholding the stop. OUT: reached in fact
  double* thr,  // sum of the off-diagonal elements
                // IN: thresholding the stop. OUT: reached in fact
  int* pnrot);  // OUT: elementary rotations performed

enum
{
  IPL_MATRIX_MULTIPLY_AB   = 0,
  IPL_MATRIX_MULTIPLY_ABt  = 1,
  IPL_MATRIX_MULTIPLY_AtB  = 2,
  IPL_MATRIX_MULTIPLY_AtBt = 3
};

EXTERNC MIA_RESULT_CODE BPL_Matrix_Multiply(
  double *R,        // OUT: result matrix       dimensions depend on mode
  const double *A,  // IN:  first multiplier    height:width = l:m
  const double *B,  // IN:  second multiplier   height:width = m:n
  int aw,           // IN:  matrix A width
  int ah,           // IN:  matrix A height
  int bw,           // IN:  matrix B width
  int bh,           // IN:  matrix B height
  int mode);        // IN:  multiplication mode

EXTERNC MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonUnitar(
  const double* A,
  int n,
  double* pdErrMaxAbsDevDiag,   // OUT: maximum absolute deviation of diagonal elements from 1
  double* pdErrMaxAbsDevNdiag,  // OUT: maximum absolute deviation of off-diagonal elements from 0
  double* pdErrAvgAbsDevDiag,   // OUT: average absolute deviation of diagonal elements from 1
  double* pdErrAvgAbsDevNdiag); // OUT: average absolute deviation of off-diagonal elements from 0

// test AQAt=E
EXTERNC MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonUnitarGeneral(
  const double* A,
  const double* Q,
  int n,
  double* pdErrMaxAbsDevDiag,   // OUT: maximum absolute deviation of diagonal elements from 1
  double* pdErrMaxAbsDevNdiag,  // OUT: maximum absolute deviation of off-diagonal elements from 0
  double* pdErrAvgAbsDevDiag,   // OUT: average absolute deviation of diagonal elements from 1
  double* pdErrAvgAbsDevNdiag); // OUT: average absolute deviation of off-diagonal elements from 0

EXTERNC MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonEigen(
  const double* A,    // IN:  source matrix
  const double* Vec,  // IN:  eigen vectors
  const double* Val,  // IN:  eigen values
  int n,              // IN:  size of source matrix
  int m,              // IN:  number of eigen vectors to be tested
  double* maxdiff);

// Ax=LBx
EXTERNC MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonEigenGeneral(
  const double* A,    // IN:  source matrix at left
  const double* B,    // IN:  source matrix at right
  const double* Vec,  // IN:  eigen vectors
  const double* Val,  // IN:  eigen values
  int n,              // IN:  size of source matrix
  int m,              // IN:  number of eigen vectors to be tested
  double* maxdiff);

// im = U*L*Vt
EXTERNC void BPL_Matrix_RestoreImageFromSVD(
  unsigned char* im,
  const double* U,
  const double* V,
  const double* L,
  int xs,
  int ys,
  int n);

//== Error Calculation ========================================================
DECL_HANDLE(HErrCalc);

typedef struct
{
  double dThreshold;  // threshold in point
  double dFAR;        // far value in point
  double dFRR;        // frr value in point
  int nSkipped;       // number of skipped items, which did no change
} SDETCurvePoint;

// create Error Calculation object
EXTERNC MIA_RESULT_CODE BPL_ERCL_Create(
  HErrCalc* phEC,                 // OUT: handle to object
  int nBiggerMeansMoreLikely);    // IN:  type of comparison value

// destroy error calculation object
EXTERNC MIA_RESULT_CODE BPL_ERCL_Free(
  HErrCalc* phEC);  // IN/OUT: handle to object / is set to NULL

// add a comparison act item to object
EXTERNC MIA_RESULT_CODE BPL_ERCL_AddItem(
  HErrCalc hEC,     // IN:  handle to object
  double dValue,    // IN:  value of comparison function
  int nIsSame);     // IN:  is it comparison of same person?

// create DET curve points
EXTERNC MIA_RESULT_CODE BPL_ERCL_GenerateDETCurve(
  HErrCalc hEC,             // IN:  handle to object
  SDETCurvePoint** ppsDCP,  // OUT: pointer to array of curve points
  int* pnPoints);           // OUT: number of curve points

// free curve memory resource
EXTERNC MIA_RESULT_CODE BPL_ERCL_FreeDETCurve(
  SDETCurvePoint** ppsDCP); // IN/OUT: handle to array / is set to NULL

//== 1D statistics calculation ==========================================
typedef struct
{
  int nUsed;
  int nAlloced;
  int nSize;
  void *pvStore;
} SCollection;

// create collection
EXTERNC MIA_RESULT_CODE BPL_HIST_CollectionCreate(
  SCollection** ppsCol,
  int nSize);

// destroy collection
EXTERNC MIA_RESULT_CODE BPL_HIST_CollectionFree(
  SCollection** ppsCol);

// add item to collection
EXTERNC MIA_RESULT_CODE BPL_HIST_CollectionAddItem(
  SCollection* psCol, // IN:  handle to object
  const void* item);  // IN:  value

typedef struct
{
  int* pnHist;
  int nDivisor;
  int nItems;
  int nLength;
} SFlexibleHistogram;

EXTERNC MIA_RESULT_CODE BPL_HIST_FlexibleCreate(
  SFlexibleHistogram** ppsFH,
  int nLen);

EXTERNC MIA_RESULT_CODE BPL_HIST_FlexibleFree(
  SFlexibleHistogram** ppsFH);

EXTERNC MIA_RESULT_CODE BPL_HIST_FlexibleAdd(
  SFlexibleHistogram* psFH,
  int nItem);

typedef struct
{
  int mass;
  int minval;
  int maxval;
  int* hist;
  int* ptr_alloc;  // pointer to alocated
} SSimpleHist;

EXTERNC MIA_RESULT_CODE BPL_SHIST_Create(
  SSimpleHist** ppSH,
  int minval,
  int maxval);

EXTERNC MIA_RESULT_CODE BPL_SHIST_Free(
  SSimpleHist** ppSH);

EXTERNC MIA_RESULT_CODE BPL_SHIST_AddItem(
  SSimpleHist* pSH,
  int val);

//== linear systems =====================================================
typedef struct
{
  int nStr;
  int nCol;
  int isSolved;
  double** A;
  double*  b;
  double*  x;
} SLinearSystem;

// maintenance
EXTERNC MIA_RESULT_CODE BPL_LIN_AllocateSystem(
  SLinearSystem* pThis,
  int nStr,
  int nCol);

EXTERNC void BPL_LIN_FreeSystem(
  SLinearSystem* pThis);

// calculations
EXTERNC MIA_RESULT_CODE BPL_LIN_SolveSystem(
  SLinearSystem* pThis);

//====== Principal Components Analysing =========================================
// ID of PC analyser
#define MEMID_SPCANALYSER GENERATE_FOURCC_ID('S','P','C','A')
//
#define CHARSIZE 8                        // sizeof of data types used
#define SHORTSIZE 16                      //
#define INTSIZE 32                        //
#define VALSIZE (CHARSIZE+1)              // size of image item (uchar-uchar gives [-255;+255])
#define VECSIZE SHORTSIZE                 // size of vector
#define SC_POW (INTSIZE-VECSIZE-VALSIZE)  // log2(number of items)
#define SUMCHUNK (1<<SC_POW)              // number of items to be summarized in one heap
#define PC_NVEC_ALIGN 4                   // number of vectors in PC should be multiple of this

typedef struct
{
  int IDim;                 // dimension of original space (should be multiple of SUMCHUNK)
  int PDim;                 // dimension of PC space (should be multiple of 16)
  unsigned char* Mean;      // mean of source cluster
  short* Vecs;              // vector table: (idim(*%)SUMCHUNK)*pdim elements
  unsigned char* Shifts;    // normalizing coefficients
} SPCAnalyser;

// set structure according to block
EXTERNC MIA_RESULT_CODE BPL_PCA_ArrangePCStruct(
  SPCAnalyser* pSPCA,       // IN:  structure to be set
  const void* inbuf,        // IN:  buffer storing the structure
  int* nBytes);             // OUT: number of bytes used in buffer

// expand the incoming data
EXTERNC MIA_RESULT_CODE BPL_PCA_Expand(
  const SPCAnalyser* pSPCA, // IN:  PCA structure
  short* coefs,             // OUT: expansion coefficients
  const unsigned char* im,  // IN:  source image
  void* extbuftmp);

//== classifier ============================================================
DECL_HANDLE(HNDClassifier);

// initialize classifier structure from data set by ndc_build function
EXTERNC MIA_RESULT_CODE BPL_NDC_Initialise(
  HNDClassifier* phNDC, // OUT: handle to classifier
  int* pnValNumber,     // OUT: number of variables
  const void* pvPacked, // IN:  packed classifier data (can be cleaned after this function exits)
  void* pvMembuf,       // IN:  memory buffer for unpacking classifier
  int* pnBuflen);       // IN/OUT: number of bytes allocated/used in buffer

// free classifier instance
EXTERNC MIA_RESULT_CODE BPL_NDC_Free(
  void** ppvMembuf,     // OUT: membuf used for classifier
  HNDClassifier* phNDC);// IN:  handle to classifier

// classify the sample
EXTERNC MIA_RESULT_CODE BPL_NDC_Classify(
  int* pnClass,         // OUT: class number
  float* pdRisks,       // OUT: risk value for each class
  HNDClassifier hNDC,   // IN:  classifier handle
  const float* pdVec);  // OUT: vector to be classified

DECL_HANDLE(HNDClassBuilder);

// initialise classifier builder
EXTERNC MIA_RESULT_CODE BPL_NDC_Create(
  HNDClassBuilder* phNDB, // OUT: handle to classifier builder
  int nClasses,           // IN:  number of classes
  int nDimensions);       // IN:  number of dimensions

// free builder instance
EXTERNC MIA_RESULT_CODE BPL_NDC_FreeBuilder(
  HNDClassBuilder* phNDB);// IN:  handle to classifier

// add sample to teaching set
EXTERNC MIA_RESULT_CODE BPL_NDC_AddSample(
  HNDClassBuilder hNDB,   // IN:  classifier builder handle
  const double* pdVec,    // IN:  sample
  int nClass);            // IN:  class number

// build the classifier
EXTERNC MIA_RESULT_CODE BPL_NDC_Build(
  HNDClassBuilder hNDB,   // IN:  classifier builder handle
  void* pvBuf,            // OUT: buffer for classifier packed data
  int* pnLen,             // IN/OUT: bytes allocated/used for storing
  const double* pdPens);  // IN:  penalty matrix

typedef int (*BPL_OPT_LevenbergMarquardt_Callback)(double *x, int m, int n, double *f);

// un-constraint minimization of differentiable function of multiple variables
// represented as sum of squares by Levenberg-Marquardt (LM) procedure
//  min F(x) , F(x)=f_1^2+...+f_m^2
EXTERNC int BPL_OPT_LevenbergMarquardt(
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
  int *infer);    // OUT: criterion of stop: 0-no convergence, 1-argument precision reached
                  //  2-function precision reached, 3-gradient precision reached

typedef double (*BPL_OPT_LocateMinimum_Callback)(double value, const void* context);

// search local minimum of a function of one variable
EXTERNC int BPL_OPT_LocateMinimum(
  BPL_OPT_LocateMinimum_Callback funct,
  const void* context,
  double *a,
  double *b,
  double *c__,
  double *eps,
  double *xmin,
  double *fmin);

#endif // __bpl_h__
