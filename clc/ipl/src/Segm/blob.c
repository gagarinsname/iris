/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  28 November 2012                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  blob detection                                         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 117 // __FILENUM__TAG117

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "bpl.h"
#include "dbg.h"
#include "cgl.h"
#include "../ownipl.h"

staticFunc int BlobDistance(
  const SBlobInfo *psBI,
  const SBlobData1 *psBD,
  int cBlob,
  int cNeib)
{
  int v;

    v  = (psBD[cBlob].nI+psBD[cBlob].nS/2)/psBD[cBlob].nS-
         (psBD[cNeib].nI+psBD[cNeib].nS/2)/psBD[cNeib].nS;
    return (v>0)?v:(-v);
}

MIA_RESULT_CODE IPL_CONOB_Enumerate_v1(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData1 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg)       // IN:  for debug, set to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  SBlobInfo *psBI = NULL;   // blob infos
  SBlobData1 *psBD = NULL;  // blob data
  int *pnNeibList = NULL;   // neighbor indices
  int x,y,cIdx,nNeibList_u,nNeibList_a;
//  char nama[FILENAME_MAX];
  int nBlob=0,cBlob,cNeib,cIdx1,cIdx2,nNeib,nIdx1,nIdx2,cIdx3;
  int nGeneration=0,nNextGenerationIdx,nMaxApproveVal,nMaxApproveIdx,tmpint;

    // check arguments
    if ((ppsBI==NULL)||(ppsBD==NULL)||(ppnNeibList==NULL)||(pnBlob==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if ((*ppsBI!=NULL)||(*ppsBD!=NULL)||(*ppnNeibList!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    // allocate
    if ((psBI = (SBlobInfo*)malloc(sizeof(psBI[0])*xs*ys*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_v1_exit;
    }
    if ((psBD = (SBlobData1*)malloc(sizeof(psBD[0])*xs*ys*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_v1_exit;
    }
    nNeibList_u = 0;
    if ((pnNeibList = (int*)malloc((nNeibList_a = xs*ys*4)*sizeof(pnNeibList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_v1_exit;
    }
    // set initial blobs from pixels
    nBlob = 0;
    for (y=0;y<ys;y++)
      for (x=0;x<xs;x++)
      {
        cIdx = y*xs+x;
        psBI[nBlob].nParent = -1;
        psBI[nBlob].nChild1 = -1;
        psBI[nBlob].nChild2 = -1;
        psBI[nBlob].nNeibStart = nNeibList_u;
        psBI[nBlob].nNeibCount = 0;
        if (y>0)
        { // upper
          pnNeibList[nNeibList_u++] = (y-1)*xs+x;
          psBI[nBlob].nNeibCount++;
        }
        if (x>0)
        { // left
          pnNeibList[nNeibList_u++] = y*xs+(x-1);
          psBI[nBlob].nNeibCount++;
        }
        if (x<xs-1)
        { // right
          pnNeibList[nNeibList_u++] = y*xs+(x+1);
          psBI[nBlob].nNeibCount++;
        }
        if (y<ys-1)
        { // upper
          pnNeibList[nNeibList_u++] = (y+1)*xs+x;
          psBI[nBlob].nNeibCount++;
        }
        psBD[nBlob].nS = 1;
        psBD[nBlob].nI = im[cIdx];
        nBlob++;
      }
    nNextGenerationIdx = nBlob;
    // run blob counter to the end of list
    for (cBlob=0;cBlob<nBlob;cBlob++)
    {
      if (cBlob==nNextGenerationIdx)
      {
        nGeneration++;
        nNextGenerationIdx = nBlob;
      }
      if (psBI[cBlob].nParent>=0)
        continue; // the object is a part of bigger one
      // test neibs for merging
      nMaxApproveVal = 0;
      nMaxApproveIdx = -1;
      for (cNeib=0;cNeib<psBI[cBlob].nNeibCount;cNeib++)
      {
        int nnei = pnNeibList[psBI[cBlob].nNeibStart+cNeib];
        if (cBlob>=nnei)
          continue;
        tmpint = 10-BlobDistance(psBI,psBD,cBlob,nnei);
        if (nMaxApproveVal<tmpint)
        {
          nMaxApproveVal = tmpint;
          nMaxApproveIdx = cNeib;
        }
      }
      if (nMaxApproveIdx==-1)
        continue; // cannot merge now, next object
      // set 'cneib' to index in global object list
      cNeib = pnNeibList[psBI[cBlob].nNeibStart+nMaxApproveIdx];
      // merge 'cblob' and 'cneib'
      // fill new blob info
      psBI[nBlob].nParent = -1;
      psBI[nBlob].nChild1 = cBlob;
      psBI[nBlob].nChild2 = cNeib;
      psBI[nBlob].nNeibStart = nNeibList_u;
      psBI[nBlob].nNeibCount = 0;
      // update parents in children of new blob
      psBI[cBlob].nParent = nBlob;
      psBI[cNeib].nParent = nBlob;
      // copy neibs from two child blobs to new one
      cIdx1 = psBI[cBlob].nNeibStart;
      nIdx1 = cIdx1+psBI[cBlob].nNeibCount;
      cIdx2 = psBI[cNeib].nNeibStart;
      nIdx2 = cIdx2+psBI[cNeib].nNeibCount;
      while ((cIdx1!=nIdx1)||(cIdx2!=nIdx2))
      {
        if (cIdx1==nIdx1)
          nNeib = pnNeibList[cIdx2++];
        else
          if (cIdx2==nIdx2)
            nNeib = pnNeibList[cIdx1++];
          else
            if (pnNeibList[cIdx1]==pnNeibList[cIdx2])
            {
              nNeib = -1;
              cIdx1++;
            }
            else
              if (pnNeibList[cIdx1]<pnNeibList[cIdx2])
                nNeib = pnNeibList[cIdx1++];
              else
                nNeib = pnNeibList[cIdx2++];
        if ((nNeib<0)||(nNeib==cBlob)||(nNeib==cNeib))
          continue;
        if (nNeibList_u==nNeibList_a)
          pnNeibList = (int*)realloc(pnNeibList,
                                (nNeibList_a+=xs*ys)*sizeof(pnNeibList[0]));
        pnNeibList[nNeibList_u++] = nNeib;
      }
      psBI[nBlob].nNeibCount = nNeibList_u-psBI[nBlob].nNeibStart;
      // fill new blob data
      psBD[nBlob].nS = psBD[cBlob].nS+psBD[cNeib].nS;
      psBD[nBlob].nI = psBD[cBlob].nI+psBD[cNeib].nI;
      // notify neibs about number changing
      for (cIdx1=psBI[nBlob].nNeibStart,nIdx1=cIdx1+psBI[nBlob].nNeibCount;cIdx1<nIdx1;cIdx1++)
      { // enumerating neibs
        nNeib = pnNeibList[cIdx1];  // number of the neib
        for (cIdx2=psBI[nNeib].nNeibStart,nIdx2=cIdx2+psBI[nNeib].nNeibCount;cIdx2<nIdx2;cIdx2++)
        { // enumerating neibs of the neibs
          if ((pnNeibList[cIdx2]==cBlob)||(pnNeibList[cIdx2]==cNeib))
          {
            for (cIdx3=cIdx2+1;cIdx3<nIdx2;cIdx3++)
              pnNeibList[cIdx3-1] = pnNeibList[cIdx3];
            nIdx2--;
            cIdx2--;
            psBI[nNeib].nNeibCount--;
          }
        }
        pnNeibList[psBI[nNeib].nNeibStart+psBI[nNeib].nNeibCount] = nBlob;
        psBI[nNeib].nNeibCount++;
      }
      nBlob++;
    }
/*    if (nambeg)
    {
      uint8* dbgim;
      strcpy(nama,nambeg);
      strcat(nama,"_0pnt.bmp");
      dbgim = (uint8*)malloc(xs*ys);
      for (cIdx=xs*ys-1;cIdx>=0;cIdx--)
      {
        for (cIdx2=cIdx;psBI[cIdx2].nParent>=0;cIdx2=psBI[cIdx2].nParent);
        dbgim[cIdx] = (uint8)((psBD[cIdx2].nI+psBD[cIdx2].nS/2)/psBD[cIdx2].nS);
      }
      DBGL_FILE_SaveUint8Image(dbgim,xs,xs,ys,nama);
      free(dbgim);
    }*/
IPL_CONOB_Enumerate_v1_exit:;
    if (res!=ERR_OK)
    {
      if (psBI)
        free(psBI);
      psBI = NULL;
      if (psBD)
        free(psBD);
      psBD = NULL;
      if (pnNeibList)
        free(pnNeibList);
      pnNeibList = NULL;
    }
    else
    {
      *ppsBI = psBI;
      *ppsBD = psBD;
      *ppnNeibList = pnNeibList;
      *pnBlob = nBlob;
    }
    return res;
}

#if 0
MIA_RESULT_CODE IPL_CONOB_Enumerate_v2(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData1 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg)       // IN:  for debug, set to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  SBlobInfo *psBI = NULL;   // blob infos
  SBlobData1 *psBD = NULL;  // blob data
  int *pnNeibList = NULL;   // neighbor indices
  int x,y,cIdx,nNeibList_u,nNeibList_a;
//  char nama[FILENAME_MAX];
  int nBlob=0,cBlob,cNeib,cIdx1,cIdx2,nNeib,nIdx1,nIdx2,cIdx3,nEra,nEraStartingBlob;
  int nGeneration=0,nNextGenerationIdx,nMaxApproveVal,nMaxApproveIdx,tmpint;

    // check arguments
    if ((ppsBI==NULL)||(ppsBD==NULL)||(ppnNeibList==NULL)||(pnBlob==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if ((*ppsBI!=NULL)||(*ppsBD!=NULL)||(*ppnNeibList!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    // allocate
    if ((psBI = (SBlobInfo*)malloc(sizeof(psBI[0])*xs*ys*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_v2_exit;
    }
    if ((psBD = (SBlobData1*)malloc(sizeof(psBD[0])*xs*ys*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_v2_exit;
    }
    nNeibList_u = 0;
    if ((pnNeibList = (int*)malloc((nNeibList_a = xs*ys*4)*sizeof(pnNeibList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_v2_exit;
    }
    // set initial blobs from pixels
    nBlob = 0;
    for (y=0;y<ys;y++)
      for (x=0;x<xs;x++)
      {
        cIdx = y*xs+x;
        psBI[nBlob].nParent = -1;
        psBI[nBlob].nChild1 = -1;
        psBI[nBlob].nChild2 = -1;
        psBI[nBlob].nNeibStart = nNeibList_u;
        psBI[nBlob].nNeibCount = 0;
        if (y>0)
        { // upper
          pnNeibList[nNeibList_u++] = (y-1)*xs+x;
          psBI[nBlob].nNeibCount++;
        }
        if (x>0)
        { // left
          pnNeibList[nNeibList_u++] = y*xs+(x-1);
          psBI[nBlob].nNeibCount++;
        }
        if (x<xs-1)
        { // right
          pnNeibList[nNeibList_u++] = y*xs+(x+1);
          psBI[nBlob].nNeibCount++;
        }
        if (y<ys-1)
        { // upper
          pnNeibList[nNeibList_u++] = (y+1)*xs+x;
          psBI[nBlob].nNeibCount++;
        }
        psBD[nBlob].nS = 1;
        psBD[nBlob].nR = im[3*cIdx+2];
        psBD[nBlob].nG = im[3*cIdx+1];
        psBD[nBlob].nB = im[3*cIdx+0];
        nBlob++;
      }
    nNextGenerationIdx = nBlob;
    nEra = 0;
    do
    {
      nEraStartingBlob = nBlob;
      // run blob counter to the end of list
      for (cBlob=0;cBlob<nBlob;cBlob++)
      {
        if (cBlob==nNextGenerationIdx)
        {
          nGeneration++;
          nNextGenerationIdx = nBlob;
        }
        if (psBI[cBlob].nParent>=0)
          continue; // the object is a part of bigger one
        // test neibs for merging
        if ( nEra&&(psBI[cBlob].nNeibCount==1)&&
            (psBD[pnNeibList[psBI[cBlob].nNeibStart]].nS<4000) )
          // only one neibor => surrounder, and with small size
          nMaxApproveIdx = 0;
        else
        {
          nMaxApproveVal = 0;
          nMaxApproveIdx = -1;
          for (cNeib=0;cNeib<psBI[cBlob].nNeibCount;cNeib++)
          {
            int nnei = pnNeibList[psBI[cBlob].nNeibStart+cNeib];
            if (cBlob>=nnei)
              continue;
            tmpint = 10-BlobDistance(psBI,psBD,cBlob,nnei);
            if (nMaxApproveVal<tmpint)
            {
              nMaxApproveVal = tmpint;
              nMaxApproveIdx = cNeib;
            }
          }
          if (nMaxApproveIdx==-1)
          {
            int nMassMax=0,nDistToMax=0,nIdxMax=-1;
            if (!nEra)
              continue; // cannot merge now, next object
            // find biggest neib
            for (cNeib=0;cNeib<psBI[cBlob].nNeibCount;cNeib++)
              if (nMassMax<psBD[pnNeibList[psBI[cBlob].nNeibStart+cNeib]].nS)
                nMassMax = psBD[nIdxMax = pnNeibList[psBI[cBlob].nNeibStart+cNeib]].nS;
            if (nIdxMax<0)
              continue;
            nDistToMax = BlobDistance(psBI,psBD,cBlob,nIdxMax);
            // check other neibs
            for (cNeib=0;cNeib<psBI[cBlob].nNeibCount;cNeib++)
            {
              if (cNeib==nIdxMax)
                continue;
              tmpint = BlobDistance(psBI,psBD,cBlob,pnNeibList[psBI[cBlob].nNeibStart+cNeib]);
              if ((tmpint+5<nDistToMax)&&(tmpint<30)&&
                  (psBD[pnNeibList[psBI[cBlob].nNeibStart+cNeib]].nS+
                   psBD[cBlob].nS<4000))
                break;
            }
            if (cNeib==psBI[cBlob].nNeibCount)
              continue;
            nMaxApproveIdx = cNeib;
          }
        }
        // set 'cneib' to index in global object list
        cNeib = pnNeibList[psBI[cBlob].nNeibStart+nMaxApproveIdx];
        // merge 'cblob' and 'cneib'
        // fill new blob info
        psBI[nBlob].nParent = -1;
        psBI[nBlob].nChild1 = cBlob;
        psBI[nBlob].nChild2 = cNeib;
        psBI[nBlob].nNeibStart = nNeibList_u;
        psBI[nBlob].nNeibCount = 0;
        // update parents in children of new blob
        psBI[cBlob].nParent = nBlob;
        psBI[cNeib].nParent = nBlob;
        // copy neibs from two child blobs to new one
        cIdx1 = psBI[cBlob].nNeibStart;
        nIdx1 = cIdx1+psBI[cBlob].nNeibCount;
        cIdx2 = psBI[cNeib].nNeibStart;
        nIdx2 = cIdx2+psBI[cNeib].nNeibCount;
        while ((cIdx1!=nIdx1)||(cIdx2!=nIdx2))
        {
          if (cIdx1==nIdx1)
            nNeib = pnNeibList[cIdx2++];
          else
            if (cIdx2==nIdx2)
              nNeib = pnNeibList[cIdx1++];
            else
              if (pnNeibList[cIdx1]==pnNeibList[cIdx2])
              {
                nNeib = -1;
                cIdx1++;
              }
              else
                if (pnNeibList[cIdx1]<pnNeibList[cIdx2])
                  nNeib = pnNeibList[cIdx1++];
                else
                  nNeib = pnNeibList[cIdx2++];
          if ((nNeib<0)||(nNeib==cBlob)||(nNeib==cNeib))
            continue;
          if (nNeibList_u==nNeibList_a)
            pnNeibList = (int*)realloc(pnNeibList,
                                  (nNeibList_a+=xs*ys)*sizeof(pnNeibList[0]));
          pnNeibList[nNeibList_u++] = nNeib;
        }
        psBI[nBlob].nNeibCount = nNeibList_u-psBI[nBlob].nNeibStart;
        // fill new blob data
        psBD[nBlob].nS = psBD[cBlob].nS+psBD[cNeib].nS;
        psBD[nBlob].nR = psBD[cBlob].nR+psBD[cNeib].nR;
        psBD[nBlob].nG = psBD[cBlob].nG+psBD[cNeib].nG;
        psBD[nBlob].nB = psBD[cBlob].nB+psBD[cNeib].nB;
        // notify neibs about number changing
        for (cIdx1=psBI[nBlob].nNeibStart,nIdx1=cIdx1+psBI[nBlob].nNeibCount;cIdx1<nIdx1;cIdx1++)
        { // enumerating neibs
          nNeib = pnNeibList[cIdx1];  // number of the neib
          for (cIdx2=psBI[nNeib].nNeibStart,nIdx2=cIdx2+psBI[nNeib].nNeibCount;cIdx2<nIdx2;cIdx2++)
          { // enumerating neibs of the neibs
            if ((pnNeibList[cIdx2]==cBlob)||(pnNeibList[cIdx2]==cNeib))
            {
              for (cIdx3=cIdx2+1;cIdx3<nIdx2;cIdx3++)
                pnNeibList[cIdx3-1] = pnNeibList[cIdx3];
              nIdx2--;
              cIdx2--;
              psBI[nNeib].nNeibCount--;
            }
          }
          pnNeibList[psBI[nNeib].nNeibStart+psBI[nNeib].nNeibCount] = nBlob;
          psBI[nNeib].nNeibCount++;
        }
        nBlob++;
      }
/*      if (nambeg)
      {
        uint8* dbgim;
        sprintf(nama,"%s_%02d_%02d.bmp",nambeg,nEra,nGeneration);
        dbgim = (uint8*)malloc(xs*ys*3);
        for (cIdx=xs*ys-1;cIdx>=0;cIdx--)
        {
          for (cIdx2=cIdx;psBI[cIdx2].nParent>=0;cIdx2=psBI[cIdx2].nParent);
          {
            dbgim[cIdx*3+2] = psBD[cIdx2].nR/psBD[cIdx2].nS;
            dbgim[cIdx*3+1] = psBD[cIdx2].nG/psBD[cIdx2].nS;
            dbgim[cIdx*3+0] = psBD[cIdx2].nB/psBD[cIdx2].nS;
          }
        }
        DBGL_FILE_SaveRGB8Image(dbgim,xs,xs,ys,nama);
        free(dbgim);
      }*/
      nEra++;
    }
    while (nEraStartingBlob!=nBlob);//(nNextGenerationIdx!=nBlob);
    // final
/*    if (nambeg)
    {
      uint8* dbgim;
      sprintf(nama,"%s_%02d_%02d_final.bmp",nambeg,nEra,nGeneration);
      dbgim = (uint8*)malloc(xs*ys*3);
      for (cIdx=xs*ys-1;cIdx>=0;cIdx--)
      {
        for (cIdx2=cIdx;psBI[cIdx2].nParent>=0;cIdx2=psBI[cIdx2].nParent);
        if ((psBD[cIdx2].nS<100)||(psBD[cIdx2].nS>5000))
        {
          dbgim[cIdx*3+2] = 0;
          dbgim[cIdx*3+1] = 0;
          dbgim[cIdx*3+0] = 0;
        }
        else
        {
          dbgim[cIdx*3+2] = psBD[cIdx2].nR/psBD[cIdx2].nS;
          dbgim[cIdx*3+1] = psBD[cIdx2].nG/psBD[cIdx2].nS;
          dbgim[cIdx*3+0] = psBD[cIdx2].nB/psBD[cIdx2].nS;
        }
      }
      DBGL_FILE_SaveRGB8Image(dbgim,xs,xs,ys,nama);
      free(dbgim);
    }*/
IPL_CONOB_Enumerate_v2_exit:;
    if (res!=ERR_OK)
    {
      if (psBI)
        free(psBI);
      psBI = NULL;
      if (psBD)
        free(psBD);
      psBD = NULL;
      if (pnNeibList)
        free(pnNeibList);
      pnNeibList = NULL;
    }
    else
    {
      *ppsBI = psBI;
      *ppsBD = psBD;
      *ppnNeibList = pnNeibList;
      *pnBlob = nBlob;
    }
    return res;
}
#endif //0
