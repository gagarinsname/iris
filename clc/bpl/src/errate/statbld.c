/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  29 October 2010                                        */
/*      Revision: 3.0.00                                                 */
/*      Purpose:  1D statistics building                                 */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 68 // __FILENUM__TAG68

#include <malloc.h>
#include "bpl.h"

/*
// sorting for correlation-like measure (for instance similar=1, different=0)
// i.e. nBiggerMeansMoreLikely = 1
#define SErClItem_IsLess_Corr(s1,s2) (\
  ((s1).dValue!=(s2).dValue)?\
  ((s1).dValue >(s2).dValue):\
  ((s1).nIsSame<(s2).nIsSame)    )
staticFunc IMPLEMENT_HEAPSORT(ownBPL_heapSortSErClItem_Corr,SErClItem,SErClItem_IsLess_Corr)
*/

// create collection
MIA_RESULT_CODE BPL_HIST_CollectionCreate(
  SCollection** ppsCol,
  int nSize)
{
  SCollection* psCol=NULL;

    // check arguments
    if (ppsCol==NULL)
      return ERR_GEN_NULLPOINTER;
    if (*ppsCol!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    if (nSize<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // allocate
    if ((psCol = (SCollection*)malloc(sizeof(psCol[0])))==NULL)
      return ERR_GEN_NOMEMORY;
    memset(psCol,0,sizeof(psCol[0]));
    // init
    psCol->nSize = nSize;
    // output
    *ppsCol = psCol;
    return ERR_OK;
}

// destroy collection
MIA_RESULT_CODE BPL_HIST_CollectionFree(
  SCollection** ppsCol)
{
  SCollection* psCol;

    if (ppsCol==NULL)
      return ERR_OK;
    psCol = *ppsCol;
    if (psCol==NULL)
      return ERR_OK;
    if (psCol->pvStore)
      free(psCol->pvStore);
    free(psCol);
    *ppsCol = NULL;
    return ERR_OK;
}

// add item to collection
MIA_RESULT_CODE BPL_HIST_CollectionAddItem(
  SCollection* psCol, // IN:  handle to object
  const void* item)   // IN:  value
{
  int tmpint;
  void* pvTemp;

    if ((psCol==NULL)||(item==NULL))
      return ERR_GEN_NULLPOINTER;
    if (psCol->nUsed==psCol->nAlloced)
    { // space exhausted, realloc
      if (psCol->nAlloced)
        tmpint = psCol->nAlloced*2;
      else
        if (psCol->nSize<(1<<12))
          tmpint = (1<<12)/psCol->nSize;
        else
          tmpint = 1;
      if ((psCol->pvStore = (void*)realloc(pvTemp = psCol->pvStore,tmpint*psCol->nSize))==NULL)
      { // realloc failed to inflate memblock and leaves it as it was
        psCol->pvStore = pvTemp;
        return ERR_GEN_NOMEMORY;
      }
      psCol->nAlloced = tmpint;
    }
    memcpy(((char*)psCol->pvStore)+psCol->nUsed*psCol->nSize,item,psCol->nSize);
    psCol->nUsed++;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_HIST_FlexibleCreate(
  SFlexibleHistogram** ppsFH,
  int nLen)
{
  SFlexibleHistogram* psFH;

    if (ppsFH==NULL)
      return ERR_GEN_NULLPOINTER;
    *ppsFH = NULL;
    if (nLen<=0)
      return ERR_GEN_INVALID_PARAMETER;
    psFH = (SFlexibleHistogram*)malloc(sizeof(psFH[0])+sizeof(psFH->pnHist[0])*nLen);
    if (psFH==NULL)
      return ERR_GEN_NOMEMORY;
    psFH->pnHist = (int*)(psFH+1);
    memset(psFH->pnHist,0,sizeof(psFH->pnHist[0])*nLen);
    psFH->nDivisor = 1;
    psFH->nItems = 0;
    psFH->nLength = nLen;
    *ppsFH = psFH;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_HIST_FlexibleFree(
  SFlexibleHistogram** ppsFH)
{
    if (ppsFH==NULL)
      return ERR_GEN_NULLPOINTER;
    if (*ppsFH==NULL)
      return ERR_GEN_NULLPOINTER;
    free(*ppsFH);
    *ppsFH = NULL;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_HIST_FlexibleAdd(
  SFlexibleHistogram* psFH,
  int nItem)
{
  int i;

    if (psFH==NULL)
      return ERR_GEN_NULLPOINTER;
    if (nItem<0)
      return ERR_GEN_INVALID_PARAMETER;
    nItem /= psFH->nDivisor;
    while (nItem>=psFH->nLength)
    {
      for (i=0;i<psFH->nLength/2;i++)
        psFH->pnHist[i] = psFH->pnHist[i*2]+psFH->pnHist[i*2+1];
      // process odd residue
      if (psFH->nLength&1)
        psFH->pnHist[psFH->nLength/2+1] = psFH->pnHist[psFH->nLength-1];
      psFH->nDivisor *= 2;
      nItem /= 2;
    }
    psFH->pnHist[nItem]++;
    psFH->nItems++;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_SHIST_Create(
  SSimpleHist** ppSH,
  int minval,
  int maxval)
{
  SSimpleHist* pSH;

    if (ppSH==NULL)
      return ERR_GEN_NULLPOINTER;
    if (*ppSH!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    if (minval+1>=maxval)
      return ERR_GEN_INVALID_PARAMETER;
    if ((pSH = (SSimpleHist*)malloc(sizeof(*pSH)))==NULL)
      return ERR_GEN_NOMEMORY;
    if ((pSH->ptr_alloc = (int*)malloc(sizeof(pSH->hist[0])*(maxval+1-minval)))==NULL)
    {
      free(pSH);
      return ERR_GEN_NOMEMORY;
    }
    memset(pSH->ptr_alloc,0,sizeof(pSH->hist[0])*(maxval+1-minval));
    pSH->minval = minval;
    pSH->maxval = maxval;
    pSH->mass = 0;
    pSH->hist = pSH->ptr_alloc-minval;
    *ppSH = pSH;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_SHIST_Free(
  SSimpleHist** ppSH)
{
  SSimpleHist* pSH;

    if (ppSH==NULL)
      return ERR_GEN_NULLPOINTER;
    pSH = *ppSH;
    if (pSH==NULL)
      return ERR_OK;
    if (pSH->ptr_alloc)
      free(pSH->ptr_alloc);
    free(pSH);
    *ppSH = NULL;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_SHIST_AddItem(
  SSimpleHist* pSH,
  int val)
{
    if (pSH==NULL)
      return ERR_GEN_NULLPOINTER;
    if (val<pSH->minval)
      pSH->hist[pSH->minval]++;
    else
      if (val>pSH->maxval)
        pSH->hist[pSH->maxval]++;
      else
        pSH->hist[val]++;
    pSH->mass++;
    return ERR_OK;
}
