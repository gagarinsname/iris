/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  9 April 2013                                           */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  linear chunks, merged by hierarchical clustering       */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 119 // __FILENUM__TAG119

#include <stdio.h>
#include <malloc.h>
#include <string.h>
//#include "dbg.h"
#include "bpl.h"
#include "../ownipl.h"

double BlobDistance_StraightChunks_v1(
  const SBlobData3* bd1,
  const SBlobData3* bd2);

// former IPL_CONOB_Enumerate_StraightChunks_v1
MIA_RESULT_CODE IPL_HIERCLUST_EnumerateLinear_Coords(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData3 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  const short* pgx,
  const short* pgy,
  int xs,                   // IN:  image size
  int ys)                   //
{
  MIA_RESULT_CODE res = ERR_OK;
  SBlobInfo *psBI = NULL;   // blob infos
  SBlobData3 *psBD = NULL;  // blob data
  int *pnNeibList = NULL;   // neighbor indices
  int *pnHelperIdxs = NULL; // neighbor indices
  int x,y,cIdx,nNeibList_u,nNeibList_a,nPts;
  int nBlob=0,cBlob,cNeib,cIdx1,cIdx2,nNeib,nIdx1,nIdx2,cIdx3,nnei;
  int nGeneration=0,nNextGenerationIdx,nMinDistIdx;
  double dDist,dMinDistVal;

    // check arguments
    if ((ppsBI==NULL)||(ppsBD==NULL)||(ppnNeibList==NULL)||(pnBlob==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if ((*ppsBI!=NULL)||(*ppsBD!=NULL)||(*ppnNeibList!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    // prepare helper array of indexes
    if ((pnHelperIdxs = (int*)malloc(xs*ys*sizeof(pnHelperIdxs[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_EnumerateLinear_Coords_exit;
    }
    memset(pnHelperIdxs,-1,xs*ys*sizeof(pnHelperIdxs[0]));
    // count number of existing points
    nPts = 0;
    for (cIdx=0;cIdx<xs*ys;cIdx++)
      if (im[cIdx])
        pnHelperIdxs[cIdx] = nPts++;
    if (!nPts)
      // no points, nothing to do
      goto HIERCLUST_EnumerateLinear_Coords_exit;
    // allocate
    if ((psBI = (SBlobInfo*)malloc(sizeof(psBI[0])*nPts*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_EnumerateLinear_Coords_exit;
    }
    if ((psBD = (SBlobData3*)malloc(sizeof(psBD[0])*nPts*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_EnumerateLinear_Coords_exit;
    }
    nNeibList_u = 0;
    if ((pnNeibList = (int*)malloc((nNeibList_a = nPts*8)*sizeof(pnNeibList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_EnumerateLinear_Coords_exit;
    }
    // set initial blobs from pixels
    nBlob = 0;
    for (y=0;y<ys;y++)
      for (x=0;x<xs;x++)
      {
        cIdx = y*xs+x;
        if (!im[cIdx])
          continue;
        psBI[nBlob].nParent = -1;
        psBI[nBlob].nChild1 = -1;
        psBI[nBlob].nChild2 = -1;
        psBI[nBlob].nNeibStart = nNeibList_u;
        psBI[nBlob].nNeibCount = 0;
//mia_error(0);//Check: if neibs should be enumerated in growing index order!
        if ((y>0)&&(im[cIdx-xs]))
        { // upper
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx-xs];
          psBI[nBlob].nNeibCount++;
        }
        if ((x>0)&&(im[cIdx-1]))
        { // left
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx-1];
          psBI[nBlob].nNeibCount++;
        }
        if ((x<xs-1)&&(im[cIdx+1]))
        { // right
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx+1];
          psBI[nBlob].nNeibCount++;
        }
        if ((y<ys-1)&&(im[cIdx+xs]))
        { // lower
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx+xs];
          psBI[nBlob].nNeibCount++;
        }
        if ((y>0)&&(x>0)&&(im[cIdx-xs-1]))
        { // upper-left
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx-xs-1];
          psBI[nBlob].nNeibCount++;
        }
        if ((y<ys-1)&&(x>0)&&(im[cIdx+xs-1]))
        { // lower-left
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx+xs-1];
          psBI[nBlob].nNeibCount++;
        }
        if ((y>0)&&(x<xs-1)&&(im[cIdx-xs+1]))
        { // upper-right
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx-xs+1];
          psBI[nBlob].nNeibCount++;
        }
        if ((y<ys-1)&&(x<xs-1)&&(im[cIdx+xs+1]))
        { // lower-right
          pnNeibList[nNeibList_u++] = pnHelperIdxs[cIdx+xs+1];
          psBI[nBlob].nNeibCount++;
        }
        psBD[nBlob].nS = 1;
        psBD[nBlob].nSx = x;
        psBD[nBlob].nSy = y;
        psBD[nBlob].nSxx = x*x;
        psBD[nBlob].nSxy = x*y;
        psBD[nBlob].nSyy = y*y;
        psBD[nBlob].nSgx = pgx[cIdx];
        psBD[nBlob].nSgy = pgy[cIdx];
        nBlob++;
      }
    // main processing
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
      dMinDistVal = 1e100;
      nMinDistIdx = -1;
      for (cNeib=0;cNeib<psBI[cBlob].nNeibCount;cNeib++)
      {
        nnei = pnNeibList[psBI[cBlob].nNeibStart+cNeib];
        if (cBlob>=nnei)
          continue;
        dDist = BlobDistance_StraightChunks_v1(psBD+cBlob,psBD+nnei);
//        if ((dMinDistVal>dDist)&&(dDist<2.))
        if ((dMinDistVal>dDist)&&(dDist<20.))
        {
          dMinDistVal = dDist;
          nMinDistIdx = cNeib;
        }
      }
      if (nMinDistIdx==-1)
        continue; // cannot merge now, next object
      // set 'cneib' to index in global object list
      cNeib = pnNeibList[psBI[cBlob].nNeibStart+nMinDistIdx];
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
      psBD[nBlob].nS   = psBD[cBlob].nS  +psBD[cNeib].nS;
      psBD[nBlob].nSx  = psBD[cBlob].nSx +psBD[cNeib].nSx;
      psBD[nBlob].nSy  = psBD[cBlob].nSy +psBD[cNeib].nSy;
      psBD[nBlob].nSxx = psBD[cBlob].nSxx+psBD[cNeib].nSxx;
      psBD[nBlob].nSxy = psBD[cBlob].nSxy+psBD[cNeib].nSxy;
      psBD[nBlob].nSyy = psBD[cBlob].nSyy+psBD[cNeib].nSyy;
      psBD[nBlob].nSgx = psBD[cBlob].nSgx+psBD[cNeib].nSgx;
      psBD[nBlob].nSgy = psBD[cBlob].nSgy+psBD[cNeib].nSgy;
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
      uint32 *dbgim,cObj,nObj;
      dbgim = (uint32*)malloc(xs*ys*4);
      memset(dbgim,0,xs*ys*4);
      for (cIdx=xs*ys-1;cIdx>=0;cIdx--)
        if (im[cIdx])
        {
          for (cIdx2=pnHelperIdxs[cIdx];psBI[cIdx2].nParent>=0;cIdx2=psBI[cIdx2].nParent);
          dbgim[cIdx] = ES_rand(cIdx2);
        }
      for (cObj=nObj=0;(int)cObj<nBlob;cObj++)
        if ((psBI[cObj].nParent<0)&&(psBD[cObj].nS>0))
          nObj++;
      sprintf(nama,"%d objects",nObj);
      DBGL_TYPE_TextInRGBA(dbgim,xs,ys,nama,0,0,0xff00,-1);
      strcpy(nama,nambeg);
      strcat(nama,"_2pnt.bmp");
      DBGL_FILE_SaveRGBA8Image(dbgim,xs,xs,ys,nama);
      free(dbgim);
    }*/
HIERCLUST_EnumerateLinear_Coords_exit:;
    if (pnHelperIdxs)
      free(pnHelperIdxs);
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
/*
typedef struct
{
  // moments
  int64 M;
  int64 Mx;
  int64 My;
  int64 Mxx;
  int64 Mxy;
  int64 Myy;
  // elliptic approximation
  double x;
  double y;
  double a;
  double b;
  double t;
  // point list
  int nPtStart;
//  int nPtCount;
} SBlobInfo2;
*/
MIA_RESULT_CODE IPL_HIERCLUST_ExtractPoints(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of clusters
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP), optional
  const SBlobInfo *psBI,    // IN:  blob infos
  const SBlobData3 *psBD,   // IN:  blob data
  int nBlob,                // IN:  number of blobs
  int xs,
  int ys)
{
  MIA_RESULT_CODE res = ERR_OK;
  int   *pnStack=NULL;
  int32 *pnIdxs=NULL;         // OUT: list of point indexes (len = total number of enumerated points)
  int32 *pnCnts=NULL;         // OUT: list of index counts (len = number of enumerated objects+1)
  uint32* pCoIm=NULL;         // OUT: 'colored' image of segments (32bpP)
  CHUNK_Linear* psCL=NULL;    // OUT: array of segments
  int nPix,cPix,cIdx,nStackDepth,nStackIdx,eIdx;
  int nPntIdx,nCntIdx,nCL=0;
  int32 S1,Sx,Sy,Sxx,Sxy,Syy,x,y;
  double co,si,rho,phi,err1,err2;

    // check arguments
    if ((ppsCL==NULL)||(ppnIdxs==NULL)||(ppnCnts==NULL)||(psBI==NULL)||(pnCL==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppsCL!=NULL)||(*ppnIdxs!=NULL)||(*ppnCnts!=NULL)||
        ((ppCoIm!=NULL)&&(*ppCoIm!=NULL)) )
      return ERR_GEN_INITIALISED_ALREADY;
    if ((nBlob<=0)||(xs<=0)||(ys<=0))
      return ERR_GEN_INVALID_PARAMETER;
    // estimate number of pixels as number of objects without children
    for (nPix=0;(psBI[nPix].nChild1<0)&&(nPix<nBlob);nPix++);
    // estimate number of pixels as mass of objects without parents
    nStackDepth = 0;
    cPix = 0;
    for (cIdx=nPix;cIdx<nBlob;cIdx++)
      if (psBI[cIdx].nParent<0)
      {
        cPix += (int)(psBD[cIdx].nS);
        if (nStackDepth<psBD[cIdx].nS)
          nStackDepth = (int)(psBD[cIdx].nS);
        nCL++;
      }
    /*if (nPix!=cPix)
    { // two ways of estimating pixel count gave different result
      res = ERR_GEN_INTERNAL;
      goto HIERCLUST_ExtractPoints_exit;
    }*/
    if (nPix==nBlob)
    { // no combined objects detected - nothing to do
      res = ERR_GEN_NO_DATA;
      goto HIERCLUST_ExtractPoints_exit;
    }
    // allocate stack
    if ((pnStack = (int*)malloc(nStackDepth*sizeof(*pnStack)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_ExtractPoints_exit;
    }
    // allocate list of point indexes
    if ((pnIdxs = (int*)malloc(nPix*sizeof(*pnIdxs)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_ExtractPoints_exit;
    }
    // allocate list of index counts
    if ((pnCnts = (int*)malloc((nBlob+1)*sizeof(*pnCnts)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_ExtractPoints_exit;
    }
    // allocate object list
    if ((psCL = (CHUNK_Linear*)malloc(nBlob*sizeof(*psCL)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto HIERCLUST_ExtractPoints_exit;
    }
    // main cycle - by objects
    nStackIdx = 0;
    nPntIdx = nCntIdx = 0;
    for (cIdx=nPix;cIdx<nBlob;cIdx++)
      if (psBI[cIdx].nParent<0)
      { // new cluster is found
        // fix start index of its points
        pnCnts[nCntIdx++] = nPntIdx;
        // enumerate points in the binary tree
        pnStack[nStackIdx++] = cIdx;
        while (nStackIdx)
        {
          eIdx = pnStack[--nStackIdx];
          if (psBI[eIdx].nChild1<0)
          { // leaf
//            pnIdxs[nPntIdx++] = eIdx;
            pnIdxs[nPntIdx++] = (int)(psBD[eIdx].nSy*xs+psBD[eIdx].nSx);
            continue;
          }
          // two branches
          pnStack[nStackIdx++] = psBI[eIdx].nChild2;
          pnStack[nStackIdx++] = psBI[eIdx].nChild1;
        } 
      }
    pnCnts[nCntIdx] = nPntIdx;
    // calculate object parameters
    for (cIdx=0;cIdx<nCntIdx;cIdx++)
    {
      S1 = Sx = Sy = Sxx = Sxy = Syy = 0;
      for (eIdx=pnCnts[cIdx];eIdx<pnCnts[cIdx+1];eIdx++)
      {
        x = pnIdxs[eIdx]%xs;
        y = pnIdxs[eIdx]/xs;
        S1++;
        Sx += x;
        Sy += y;
        Sxx += x*x;
        Sxy += x*y;
        Syy += y*y;
      }
      psCL[cIdx].S1 = S1;
      psCL[cIdx].id = cIdx+1;
      psCL[cIdx].xc = ((double)Sx)/S1;
      psCL[cIdx].yc = ((double)Sy)/S1;
      phi = .5*atan2(2*(S1*(double)Sxy-Sx*(double)Sy),
                     S1*((double)Sxx-Syy)-Sx*(double)Sx+Sy*(double)Sy);
      co = cos(phi);
      si = sin(phi);
      rho = (co*Sx+si*Sy)/S1;
      err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
      rho = (-si*Sx+co*Sy)/S1;
      err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho*rho*S1;
      psCL[cIdx].phi = phi; // angle
      psCL[cIdx].len = sqrt((double)S1);; // length (diameter)
      psCL[cIdx].err = ((err1<err2)?err1:err2)/S1; // error (MSD per point)
    }
    // picture
    if (ppCoIm)
    {
      if ((pCoIm = (uint32*)malloc(xs*ys*sizeof(*pCoIm)))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto HIERCLUST_ExtractPoints_exit;
      }
      memset(pCoIm,0,xs*ys*sizeof(*pCoIm));
      for (cIdx=0;cIdx<nCntIdx;cIdx++)
      {
        S1 = ES_rand(psCL[cIdx].id);
        for (eIdx=pnCnts[cIdx];eIdx<pnCnts[cIdx+1];eIdx++)
          pCoIm[pnIdxs[eIdx]] = S1;
      }
    }
HIERCLUST_ExtractPoints_exit:;
    // release resources
    if (pnStack)
      free(pnStack);
    if (res!=ERR_OK)
    {
      if (psCL)
        free(psCL);
      psCL = NULL;
      if (pnIdxs)
        free(pnIdxs);
      pnIdxs = NULL;
      if (pnCnts)
        free(pnCnts);
      pnCnts = NULL;
      if (pCoIm)
        free(pCoIm);
      pCoIm = NULL;
    }
    *ppsCL = psCL;
    *pnCL = nCL;
    *ppnIdxs = pnIdxs;
    *ppnCnts = pnCnts;
    *ppCoIm = pCoIm;
    return res;
}
/*
MIA_RESULT_CODE IPL_CONOB_ExtractObjects(
  SBlobInfo2** ppsBI2,      // OUT: extracted objects
  int** ppnPointList,       // OUT: list of extracted points
  int* pnBlob2,             // OUT: number of extracted blobs
  const SBlobInfo *psBI,    // IN:  source blob infos
  SBlobData1 *psBD,         // IN:  source blob data (not changed indeed)
  int nBlob1,               // IN:  number of source blobs
  int xs,                   // IN:  image size
  int ys,                   // 
  const char* nambeg)       // for debug, set it to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  SBlobInfo2* psBI2=NULL,*pbi2;
  int* pnPointList=NULL;
  int nPoint=0,cBlob,nBlob2=0,nPix,nBlob2_a=0,tmpint,cPix,cIdx,x,y;

    // check arguments
    if ((ppsBI2==NULL)||(ppnPointList==NULL)||(pnBlob2==NULL)||
        (psBI==NULL)||(psBD==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppsBI2!=NULL)||(*ppnPointList!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    if (nBlob1<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // estimate number of pixels (objects without children)
    for (nPix=0;(psBI[nPix].nChild1<0)&&(nPix<nBlob1);nPix++);
    if (nPix==nBlob1)
    { // no combined objects detected
      res = ERR_GEN_INTERNAL;
      goto IPL_CONOB_ExtractObjects_exit;
    }
    if (nPix!=xs*ys)
    { // numbers of pixels do not match
      res = ERR_GEN_INTERNAL;
      goto IPL_CONOB_ExtractObjects_exit;
    }
    // estimate numbers of objects and points in them
    for (cBlob=nPix;cBlob<nBlob1;cBlob++)
      if ((psBI[cBlob].nParent<0)&&(psBD[cBlob].nS>500)&&(psBD[cBlob].nS<5000))
      {
        if (nBlob2==nBlob2_a)
        {
          tmpint = (nBlob2_a)?(nBlob2_a*2):256;
          if ((psBI2 = (SBlobInfo2*)realloc(psBI2,tmpint*sizeof(psBI2[0])))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto IPL_CONOB_ExtractObjects_exit;
          }
          memset(psBI2+nBlob2_a,0,(tmpint-nBlob2_a)*sizeof(psBI2[0]));
          nBlob2_a = tmpint;
        }
        psBI2[nBlob2].nPtStart = nPoint;
        psBD[cBlob].nR = nBlob2;    // index
        psBD[cBlob].nG = -1;        // mark to quickly recognize
        // increase counters
        nPoint += psBD[cBlob].nS;
        nBlob2++;
      }
    if ((nBlob2==0)||(nPoint==0))
    { // no objects, but return OK
      if (psBI2)
        free(psBI2);
      psBI2 = NULL;
      nBlob2 = 0;
      goto IPL_CONOB_ExtractObjects_exit;
    }
    // reallocate to save space
    if ((psBI2 = (SBlobInfo2*)realloc(psBI2,nBlob2*sizeof(psBI2[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_ExtractObjects_exit;
    }
    // allocate point list
    if ((pnPointList = (int*)malloc(nPoint*sizeof(pnPointList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_ExtractObjects_exit;
    }
    // run through image, collect points
    for (cPix=0;cPix<nPix;cPix++)
    {
      // get index of top parent
      for (cIdx=cPix;psBI[cIdx].nParent>=0;cIdx=psBI[cIdx].nParent);
      if (psBD[cIdx].nG>=0)
        continue;
      // this is the point belonging to object of interest
      // get Cartesian coordinates from index
      x = cPix%xs;
      y = cPix/xs;
      // shortcut
      pbi2 = psBI2+psBD[cIdx].nR;
      // add point
      pnPointList[(int)(pbi2->M)] = cPix;
      // add moments
      pbi2->M++;
      pbi2->Mx += x;
      pbi2->My += y;
      pbi2->Mxx += x*x;
      pbi2->Mxy += x*y;
      pbi2->Myy += y*y;
    }
    // assign point counts (redundant, but helpful)
//    for (cBlob=0;cBlob<nBlob2-1;cBlob++)
  //    psBI2[cBlob].nPtCount = psBI2[cBlob+1].nPtStart-psBI2[cBlob].nPtStart;
    //psBI2[nBlob2-1].nPtCount = nPoint;
    // calculate approximating ellipse params
    for (cBlob=0;cBlob<nBlob2;cBlob++)
    {
      int64 Mxx,Mxy,Myy;
      double tilt,coco,sisi,cosi,a,b;
      pbi2 = psBI2+cBlob;
      pbi2->x = ((double)pbi2->Mx)/pbi2->M;
      pbi2->y = ((double)pbi2->My)/pbi2->M;
      // normalise second order moments
      Mxx = pbi2->M*pbi2->Mxx-pbi2->Mx*pbi2->Mx;
      Mxy = pbi2->M*pbi2->Mxy-pbi2->Mx*pbi2->My;
      Myy = pbi2->M*pbi2->Myy-pbi2->My*pbi2->My;
      if ((Mxx<=0)||(Myy<=0))
        return ERR_GEN_INVALID_PARAMETER;
      // calculate tilt, fairly simple with ATAN2
      tilt = atan2(2.*((double)Mxy),(double)(Mxx-Myy))/2.;
      if (tilt<0)
        tilt += 2*PI;
      if (tilt>PI)
        tilt -= PI;
      pbi2->t = tilt;
      //
      coco = cos(tilt);
      sisi = sin(tilt);
      cosi = 2*coco*sisi;
      coco *= coco;
      sisi *= sisi;
      a = coco*((double)Mxx)+cosi*((double)Mxy)+sisi*((double)Myy);
      b = sisi*((double)Mxx)-cosi*((double)Mxy)+coco*((double)Myy);
//      pbi2->a = sqrt(pbi2->M*a/b/PI);
//      pbi2->b = sqrt(pbi2->M*b/a/PI);
      pbi2->a = sqrt(PI*a)/pbi2->M;
      pbi2->b = sqrt(PI*b)/pbi2->M;

//      a = pbi2->a;
  //    b = pbi2->b;
    //  pbi2->a = sqrt(pbi2->M*a/PI/b);
      //pbi2->b = sqrt(pbi2->M*b/PI/a);
    }
IPL_CONOB_ExtractObjects_exit:;
    if (res!=ERR_OK)
    {
      if (psBI2)
        free(psBI2);
      psBI2 = NULL;
      if (pnPointList)
        free(pnPointList);
      pnPointList = NULL;
      nBlob2 = 0;
    }
    *ppsBI2 = psBI2;
    *ppnPointList = pnPointList;
    *pnBlob2 = nBlob2;
    return res;
}
*/
// get (a) linear (straight) chunks by (b) hierarchical clustering from (c) coordinates
MIA_RESULT_CODE IPL_CHUNK_Linear_HierClust_Coords(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr)             // IN:  object mass threshold
{
  MIA_RESULT_CODE res=ERR_OK;
  SBlobInfo *psBI=NULL;       // blob infos
  SBlobData3 *psBD=NULL;      // blob data
  int *pnNeibList=NULL;       // neighbor indices
  int nBlob;                  // number of blobs
  int nCL=0;                  // number of clusters
  CHUNK_Linear* psCL=NULL;    // array of segments
  int32* pnIdxs=NULL;         // list of point indexes (len = total number of enumerated points)
  int32* pnCnts=NULL;         // list of index counts (len = number of enumerated objects+1)
  uint32* pCoIm=NULL;         // 'colored' image of segments (32bpP)

    // check arguments
    if ((ppsCL==NULL)||(pnCL==NULL)||(ppnIdxs==NULL)||(ppnCnts==NULL)||(ppCoIm==NULL)||
        (im==NULL)||(gxim==NULL)||(gyim==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppsCL!=NULL)||(*ppnIdxs!=NULL)||(*ppnCnts!=NULL)||(*ppCoIm!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    // create cluster hierarchy
    res = IPL_HIERCLUST_EnumerateLinear_Coords(
      &psBI,        // OUT: blob infos
      &psBD,        // OUT: blob data
      &pnNeibList,  // OUT: neighbor indices
      &nBlob,       // OUT: number of blobs
      im,           // IN:  source image
      gxim,
      gyim,
      xs,           // IN:  image size
      ys);          //
    if (res!=ERR_OK)
      goto IPL_CHUNK_Linear_HierClust_Coords_exit;
    // gain points
    res = IPL_HIERCLUST_ExtractPoints(
      &psCL,    // OUT: array of segments
      &nCL,
      &pnIdxs,  // OUT: list of point indexes (len = total number of enumerated points)
      &pnCnts,  // OUT: list of index counts (len = number of enumerated objects+1)
      &pCoIm,   // OUT: 'colored' image of segments (32bpP), optional
      psBI,     // IN:  blob infos
      psBD,     // IN:  blob data
      nBlob,    // IN:  number of blobs (a.k.a. nCL)
      xs,
      ys);
    if (res!=ERR_OK)
      goto IPL_CHUNK_Linear_HierClust_Coords_exit;
IPL_CHUNK_Linear_HierClust_Coords_exit:;
    // free resources
    if (psBI)
      free(psBI);
    if (psBD)
      free(psBD);
    if (pnNeibList)
      free(pnNeibList);
    if (res!=ERR_OK)
    {
      if (psCL)
        free(psCL);
      psCL = NULL;
      if (pnIdxs)
        free(pnIdxs);
      pnIdxs = NULL;
      if (pnCnts)
        free(pnCnts);
      pnCnts = NULL;
      if (pCoIm)
        free(pCoIm);
      pCoIm = NULL;
    }
    *ppsCL = psCL;      // OUT: array of segments
    *pnCL = nCL;        // OUT: number of segments
    *ppnIdxs = pnIdxs;  // OUT: list of point indexes (len = total number of enumerated points)
    *ppnCnts = pnCnts;  // OUT: list of index counts (len = number of enumerated objects+1)
    *ppCoIm = pCoIm;    // OUT: 'colored' image of segments (32bpP)
    return res;
}
