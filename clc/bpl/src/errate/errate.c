/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  27 September 2006                                      */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  error calculation (FAR/FRR/DET)                        */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 67 // __FILENUM__TAG67

#include <math.h>
#include <malloc.h>
#include "bpl.h"

typedef struct
{
  double dValue;
  int nIsSame;
} SErClItem;

typedef struct
{
  int nCurrentItems;
  int nAllocedItems;
  int nBiggerMeansMoreLikely;
  int nSortedToPosition;
  SErClItem* psECI;
} SErrCalc;

// sorting for correlation-like measure (for instance similar=1, different=0)
// i.e. nBiggerMeansMoreLikely = 1
#define SErClItem_IsLess_Corr(s1,s2) (\
  ((s1).dValue!=(s2).dValue)?\
  ((s1).dValue >(s2).dValue):\
  ((s1).nIsSame<(s2).nIsSame)    )
staticFunc IMPLEMENT_HEAPSORT(ownBPL_heapSortSErClItem_Corr,SErClItem,SErClItem_IsLess_Corr)

// sorting for distance-like measure (for instance similar=0, different=1)
// i.e. nBiggerMeansMoreLikely = 0
#define SErClItem_IsLess_Dist(s1,s2) (\
  ((s1).dValue!=(s2).dValue)?\
  ((s1).dValue <(s2).dValue):\
  ((s1).nIsSame<(s2).nIsSame)    )
staticFunc IMPLEMENT_HEAPSORT(ownBPL_heapSortSErClItem_Dist,SErClItem,SErClItem_IsLess_Dist)

#define BPL_ERCL_ALLOC_CHUNK 1024

// create Error Calculation object
MIA_RESULT_CODE BPL_ERCL_Create(
  HErrCalc* phEC,                 // OUT: handle to object
  int nBiggerMeansMoreLikely)     // IN:  type of comparison value
{
  SErrCalc* psEC;

    if (phEC==NULL)
      return ERR_GEN_NULLPOINTER;
    *phEC = NULL;
    if ((psEC = (SErrCalc*)malloc(sizeof(psEC[0])))==NULL)
      return ERR_GEN_NOMEMORY;
    if ((psEC->psECI = (SErClItem*)malloc(sizeof(psEC->psECI[0])*BPL_ERCL_ALLOC_CHUNK))==NULL)
    {
      free(psEC);
      return ERR_GEN_NOMEMORY;
    }
    psEC->nCurrentItems = 0;
    psEC->nAllocedItems = BPL_ERCL_ALLOC_CHUNK;
    psEC->nBiggerMeansMoreLikely = nBiggerMeansMoreLikely;
    psEC->nSortedToPosition = 0;
    *phEC = (HErrCalc)psEC;
    return ERR_OK;
}

// destroy error calculation object
MIA_RESULT_CODE BPL_ERCL_Free(
  HErrCalc* phEC)   // IN/OUT: handle to object / is set to NULL
{
    if (phEC==NULL)
      return ERR_OK;
    if (*phEC==NULL)
      return ERR_OK;
    if (((SErrCalc*)(*phEC))->psECI)
      free(((SErrCalc*)(*phEC))->psECI);
    free(*phEC);
    *phEC = NULL;
    return ERR_OK;
}

// add a comparison act item to object
MIA_RESULT_CODE BPL_ERCL_AddItem(
  HErrCalc hEC,     // IN:  handle to object
  double dValue,    // IN:  value of comparison function
  int nIsSame)      // IN:  is it comparison of same person?
{
  SErClItem* psECI_temp;
  SErrCalc* psEC = (SErrCalc*)hEC;

    if (psEC==NULL)
      return ERR_GEN_NULLPOINTER;
    if (psEC->nCurrentItems==psEC->nAllocedItems)
    { // space exhausted, realloc
      if ((psEC->psECI = (SErClItem*)realloc(psECI_temp = psEC->psECI,(psEC->nAllocedItems)*2*sizeof(SErClItem)))==NULL)
      { // realloc failed to inflat memblock and leaves it as it was
        psEC->psECI = psECI_temp;
        // try to reserve less memory
        if ((psEC->psECI = (SErClItem*)realloc(psECI_temp = psEC->psECI,(size_t)((psEC->nAllocedItems)*1.2*sizeof(SErClItem))))==NULL)
        { // realloc failed to inflate memblock and leaves it as it was
          psEC->psECI = psECI_temp;
          // try to reserve less memory
          return ERR_GEN_NOMEMORY;
        }
        else
          psEC->nAllocedItems = (int)(psEC->nAllocedItems*1.2);
       }
      else
        psEC->nAllocedItems *= 2;
    }
    psEC->psECI[psEC->nCurrentItems].dValue = dValue;
    psEC->psECI[psEC->nCurrentItems].nIsSame = ((nIsSame)?1:0);
    psEC->nCurrentItems++;
    return ERR_OK;
}

// create DET curve points
MIA_RESULT_CODE BPL_ERCL_GenerateDETCurve(
  HErrCalc hEC,             // IN:  handle to object
  SDETCurvePoint** ppsDCP,  // OUT: pointer to array of curve points
  int* pnPoints)            // OUT: number of curve points
{
  int i,k_p,nPointsAlloced,nPointsSet;
  int     nalien,calien,calien_p,calien_pp;
  int     nsame ,csame ,csame_p ,csame_pp;
  double         thr   ,thr_p   ,thr_pp;
  SDETCurvePoint *psDCP,*psDCP_temp;
  SErrCalc* psEC = (SErrCalc*)hEC;

    if ((psEC==NULL)||(ppsDCP==NULL)||(pnPoints==NULL))
      return ERR_GEN_NULLPOINTER;
    if (psEC->nCurrentItems<2)
      return ERR_GEN_NOT_INITIALISED;
    *ppsDCP = NULL;
    // allocate curve points
    if ((psDCP = (SDETCurvePoint*)malloc((nPointsAlloced = BPL_ERCL_ALLOC_CHUNK)*sizeof(psDCP[0])))==NULL)
      return ERR_GEN_NOMEMORY;
    nPointsSet = 0;
    if (psEC->nSortedToPosition<psEC->nCurrentItems)
    { // sorting is needed
      if (psEC->nBiggerMeansMoreLikely)
        ownBPL_heapSortSErClItem_Corr(psEC->psECI,psEC->nCurrentItems);
      else
        ownBPL_heapSortSErClItem_Dist(psEC->psECI,psEC->nCurrentItems);
      psEC->nSortedToPosition = psEC->nCurrentItems;
    }
    // calculate number of aliens and sames
    for (nsame=0,i=psEC->nCurrentItems-1;i>=0;i--)
      nsame += psEC->psECI[i].nIsSame;
    nalien = psEC->nCurrentItems-nsame;
    if ((!nalien)||(!nsame))
      return ERR_GEN_NO_DATA;
    // 
//    csame = csame_p = csame_pp = 0;             // number of false accepts
  //  calien = calien_p = calien_pp = nalien;     // number of false rejects
    csame = csame_p = csame_pp = nsame;             // number of false accepts
    calien = calien_p = calien_pp = 0;     // number of false rejects
    thr = thr_p = thr_pp = psEC->psECI[0].dValue;  // sliding threshold
    k_p = 0;
    // main cycle for all comparisons
    for (i=0;i<psEC->nCurrentItems;i++)
    {
      // update queue
      thr_pp = thr_p;
      thr_p = thr;
      csame_pp = csame_p;
      csame_p = csame;
      calien_pp = calien_p;
      calien_p = calien;
      // update counters
      thr = psEC->psECI[i].dValue;
      if (psEC->psECI[i].nIsSame)
//        csame++;
        csame--;
      else
//        calien--;
        calien++;
      // is there a difference?
//      if ((csame_pp==csame_p)&&(csame_p==csame))
      if (csame_pp==csame)
        continue;
//      if ((calien_pp==calien_p)&&(calien_p==calien))
      if (calien_pp==calien)
        continue;
      // add to list
      if (nPointsSet==nPointsAlloced)
        // reallocate
        if ((psDCP = (SDETCurvePoint*)realloc(psDCP_temp = psDCP,(nPointsAlloced += BPL_ERCL_ALLOC_CHUNK)*sizeof(psDCP[0])))==NULL)
        {
          free(psDCP_temp);
          return ERR_GEN_NOMEMORY;
        }
      psDCP[nPointsSet].dThreshold = thr_p;
      psDCP[nPointsSet].dFRR = (double)csame_p/nsame;
      psDCP[nPointsSet].dFAR = (double)calien_p/nalien;
      psDCP[nPointsSet].nSkipped = i-k_p;
      nPointsSet++;
      // remember new position to calc skipped
      k_p = i;
    }
    *ppsDCP = psDCP;
    *pnPoints = nPointsSet;
    return ERR_OK;
}

// free curve memory resource
MIA_RESULT_CODE BPL_ERCL_FreeDETCurve(
  SDETCurvePoint** ppsDCP)  // IN/OUT: handle to array / is set to NULL
{
    if (ppsDCP==NULL)
      return ERR_OK;
    if (*ppsDCP==NULL)
      return ERR_OK;
    free(*ppsDCP);
    *ppsDCP = NULL;
    return ERR_OK;
}
