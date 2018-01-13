/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  28 November 2012                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  construct chunks by iterative object merging           */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 116 // __FILENUM__TAG116

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "dbg.h"
#include "bpl.h"
#include "../ownipl.h"

double BlobDistance_StraightChunks_v1(
  const SBlobData3* bd1,
  const SBlobData3* bd2)
{
    if (bd1->nS+bd2->nS<3)
    { // small units, decide by gradient
      double G1,G2,P,gx1,gy1,gx2,gy2;
      gx1 = (double)(bd1->nSgx)/(double)(bd1->nS);
      gy1 = (double)(bd1->nSgy)/(double)(bd1->nS);
      gx2 = (double)(bd2->nSgx)/(double)(bd2->nS);
      gy2 = (double)(bd2->nSgy)/(double)(bd2->nS);
      G1 = gx1*gx1+gy1*gy1;
      G2 = gx2*gx2+gy2*gy2;
      if ((G1>0)&&(G2>0))
      {
        P  = gx1*gx2+gy1*gy2;
        P /= sqrt(G1*G2);
      }
      else
        P = -1.;
      return 1.-P;
    }
    { // big units, decide by MSD from straight
      double phi,S1,Sx,Sy,Sxx,Sxy,Syy,co,si,rho,err1,err2;
      S1  = (double)(bd1->nS  +bd2->nS);
      Sx  = (double)(bd1->nSx +bd2->nSx);
      Sy  = (double)(bd1->nSy +bd2->nSy);
      Sxx = (double)(bd1->nSxx+bd2->nSxx);
      Sxy = (double)(bd1->nSxy+bd2->nSxy);
      Syy = (double)(bd1->nSyy+bd2->nSyy);
      phi = .5*atan2(2*(S1*Sxy-Sx*Sy),S1*(Sxx-Syy)-Sx*Sx+Sy*Sy);
      co = cos(phi);
      si = sin(phi);
      rho = (co*Sx+si*Sy)/S1;
      err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
      rho = (-si*Sx+co*Sy)/S1;
      err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho*rho*S1;
      //if ((err1/S1>=.2)&&(err2/S1>=.2))
      return (err1<err2)?err1:err2;
    }
}

MIA_RESULT_CODE IPL_CONOB_Enumerate_StraightChunks_v1(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData3 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  const short* pgx,
  const short* pgy,
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg)       // IN:  for debug, set to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  SBlobInfo *psBI = NULL;   // blob infos
  SBlobData3 *psBD = NULL;  // blob data
  int *pnNeibList = NULL;   // neighbor indices
  int *pnHelperIdxs = NULL; // neighbor indices
  int x,y,cIdx,nNeibList_u,nNeibList_a,nPts;
//  char nama[FILENAME_MAX];
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
    // prepare helper array of indexes, count number of existing points
    if ((pnHelperIdxs = (int*)malloc(xs*ys*sizeof(pnHelperIdxs[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    }
    memset(pnHelperIdxs,-1,xs*ys*sizeof(pnHelperIdxs[0]));
    nPts = 0;
    for (cIdx=0;cIdx<xs*ys;cIdx++)
      if (im[cIdx])
        pnHelperIdxs[cIdx] = nPts++;
    if (!nPts)
      // no points, nothing to do
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    // allocate
    if ((psBI = (SBlobInfo*)malloc(sizeof(psBI[0])*nPts*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    }
    if ((psBD = (SBlobData3*)malloc(sizeof(psBD[0])*nPts*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    }
    nNeibList_u = 0;
    if ((pnNeibList = (int*)malloc((nNeibList_a = nPts*8)*sizeof(pnNeibList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
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
        if ((dMinDistVal>dDist)&&(dDist<2.))
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
if (nambeg)
{
  char nama[FILENAME_MAX];
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
}
IPL_CONOB_Enumerate_StraightChunks_v1_exit:;
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

staticFunc double BlobDistance_Archunks_v1(
  const SBlobData4* bd1,
  const SBlobData4* bd2)
{
  double dist;

    if (bd1->nS+bd2->nS<3)
/*    { // small units, decide by gradient
      double G1,G2,P,gx1,gy1,gx2,gy2;
      gx1 = (double)(bd1->nSgx)/(double)(bd1->nS);
      gy1 = (double)(bd1->nSgy)/(double)(bd1->nS);
      gx2 = (double)(bd2->nSgx)/(double)(bd2->nS);
      gy2 = (double)(bd2->nSgy)/(double)(bd2->nS);
      G1 = gx1*gx1+gy1*gy1;
      G2 = gx2*gx2+gy2*gy2;
      if ((G1>0)&&(G2>0))
      {
        P  = gx1*gx2+gy1*gy2;
        P /= sqrt(G1*G2);
      }
      else
        P = -1.;
      return 1.-P;
    }*/
    { // small units, decide by distance
      double x1,x2,y1,y2,m1,m2;
      m1 = (double)(bd1->nS);
      m2 = (double)(bd2->nS);
      x1 = bd1->nSx/m1;
      y1 = bd1->nSy/m1;
      x2 = bd2->nSx/m2;
      y2 = bd2->nSy/m2;
      x1 -= x2;
      y1 -= y2;
      dist = x1*x1+y1*y1;
      return dist/5.;
    }
    if (bd1->nS+bd2->nS<50)
    { // big units, decide by MSD from straight
      double phi,S1,Sx,Sy,Sxx,Sxy,Syy,co,si,rho,err1,err2,rho2;
      S1  = (double)(bd1->nS  +bd2->nS);
      Sx  = (double)(bd1->nSx +bd2->nSx);
      Sy  = (double)(bd1->nSy +bd2->nSy);
      Sxx = (double)(bd1->nSxx+bd2->nSxx);
      Sxy = (double)(bd1->nSxy+bd2->nSxy);
      Syy = (double)(bd1->nSyy+bd2->nSyy);
      phi = .5*atan2(2*(S1*Sxy-Sx*Sy),S1*(Sxx-Syy)-Sx*Sx+Sy*Sy);
      co = cos(phi);
      si = sin(phi);
      rho = (co*Sx+si*Sy)/S1;
      err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
      rho2 = (-si*Sx+co*Sy)/S1;
      err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho2*rho2*S1;
      dist = ((err1<err2)?err1:err2)/S1;
      return dist/5.;
    }

    { // enough point to verify circularity
      SLinearSystem pLinSys3;
      double u,x0,y0,err1,r;
      int64 S,Sx,Sy,Sxx,Sxy,Syy,Sxxx,Sxxy,Sxyy,Syyy,Sxxxx,Sxxyy,Syyyy;
      int64 dx,dy;
      BPL_LIN_AllocateSystem(&pLinSys3,3,3);
      S = bd1->nS+bd2->nS;
      Sx = bd1->nSx+bd2->nSx;
      Sy = bd1->nSy+bd2->nSy;
      Sxx = bd1->nSxx+bd2->nSxx;
      Sxy = bd1->nSxy+bd2->nSxy;
      Syy = bd1->nSyy+bd2->nSyy;
      Sxxx = bd1->nSxxx+bd2->nSxxx;
      Sxxy = bd1->nSxxy+bd2->nSxxy;
      Sxyy = bd1->nSxyy+bd2->nSxyy;
      Syyy = bd1->nSyyy+bd2->nSyyy;
      Sxxxx = bd1->nSxxxx+bd2->nSxxxx;
      Sxxyy = bd1->nSxxyy+bd2->nSxxyy;
      Syyyy = bd1->nSyyyy+bd2->nSyyyy;
      { // shift
        dx = -(Sx+S/2)/S;
        dy = -(Sy+S/2)/S;
        Sxxxx += 4*dx*Sxxx+6*dx*dx*Sxx+4*dx*dx*dx*Sx+dx*dx*dx*dx*S;
        Sxxyy += 2*dx*Sxyy+dx*dx*Syy+2*dy*Sxxy+4*dx*dy*Sxy+2*dx*dx*dy*Sy+dy*dy*Sxx+2*dx*dy*dy*Sx+dx*dx*dy*dy*S;
        Syyyy += 4*dy*Syyy+6*dy*dy*Syy+4*dy*dy*dy*Sy+dy*dy*dy*dy*S;
        Sxxx  += 3*dx*Sxx+3*dx*dx*Sx+dx*dx*dx*S;
        Sxxy  += 2*dx*Sxy+dx*dx*Sy+dy*Sxx+2*dx*dy*Sx+dx*dx*dy*S;
        Sxyy  += dx*Syy+2*dy*Sxy+2*dx*dy*Sy+dy*dy*Sx+dx*dy*dy*S;
        Syyy  += 3*dy*Syy+3*dy*dy*Sy+dy*dy*dy*S;
        Sxx   += 2*dx*Sx+dx*dx*S;
        Sxy   += dx*Sy+dy*Sx+dx*dy*S;
        Syy   += 2*dy*Sy+dy*dy*S;
        Sx    += dx*S;
        Sy    += dy*S;
      }
      pLinSys3.A[0][0] = (double)S;
      pLinSys3.A[0][1] = (double)Sx;
      pLinSys3.A[0][2] = (double)Sy;
      pLinSys3.A[1][0] = (double)Sx;
      pLinSys3.A[1][1] = (double)Sxx;
      pLinSys3.A[1][2] = (double)Sxy;
      pLinSys3.A[2][0] = (double)Sy;
      pLinSys3.A[2][1] = (double)Sxy;
      pLinSys3.A[2][2] = (double)Syy;
      pLinSys3.b[0] = (double)(Sxx+Syy);
      pLinSys3.b[1] = (double)(Sxxx+Sxyy);
      pLinSys3.b[2] = (double)(Sxxy+Syyy);
      BPL_LIN_SolveSystem(&pLinSys3);
      u  = .5*pLinSys3.x[0];
      x0 = .5*pLinSys3.x[1];
      y0 = .5*pLinSys3.x[2];
      BPL_LIN_FreeSystem(&pLinSys3);
      r = 2*u+x0*x0+y0*y0;
      if (r>0)
        r = sqrt(r);
/*      err1 = x0*x0+y0*y0;
      err1 = sqrt(2*u+x0*x0+y0*y0)-sqrt(err1);
      err1 *=10;*/
      err1 =   4*u*u*S
              +4*x0*x0*Sxx
              +4*y0*y0*Syy
              +8*u*x0*Sx
              +8*u*y0*Sy
              +8*x0*y0*Sxy
              -4*u*(Sxx+Syy)
              -4*x0*(Sxxx+Sxyy)
              -4*y0*(Sxxy+Syyy)
              +Sxxxx+Syyyy+2*Sxxyy;
      x0 -= dx;
      y0 -= dy;
      err1 = sqrt(err1/S);
      dist = err1/r;
      return 10*dist/sqrt((double)S);
    }
}

staticFunc double BlobDistance_Archunks_v2(
  const SBlobData4* bd1,
  const SBlobData4* bd2)
{
  double dist1,dist2;
  int64 S,Sx,Sy,Sxx,Sxy,Syy,Sxxx,Sxxy,Sxyy,Syyy,Sxxxx,Sxxyy,Syyyy;

    if (bd1->nS+bd2->nS<3)
    { // small units, decide by distance
      double x,y,m1,m2;
      m1 = (double)(bd1->nS);
      m2 = (double)(bd2->nS);
      x = bd1->nSx/m1-bd2->nSx/m2;
      y = bd1->nSy/m1-bd2->nSy/m2;
      dist1 = x*x+y*y;
      return dist1/5.;
    }
    { // calculate MSD from straight
      double phi,co,si,rho,err1,err2,rho2;
      S   = bd1->nS  +bd2->nS;
      Sx  = bd1->nSx +bd2->nSx;
      Sy  = bd1->nSy +bd2->nSy;
      Sxx = bd1->nSxx+bd2->nSxx;
      Sxy = bd1->nSxy+bd2->nSxy;
      Syy = bd1->nSyy+bd2->nSyy;
      phi = .5*atan2(2.*(S*Sxy-Sx*Sy),(double)(S*(Sxx-Syy)-Sx*Sx+Sy*Sy));
      co = cos(phi);
      si = sin(phi);
      rho = (co*Sx+si*Sy)/S;
      err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S;
      rho2 = (-si*Sx+co*Sy)/S;
      err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho2*rho2*S;
      dist1 = ((err1<err2)?err1:err2)/S;
      dist1 /= 5.;
    }
    if (bd1->nS+bd2->nS<10)
      return dist1;
    { // enough points to verify circularity
      double u,x0,y0,err1,r;
      int64 A,B,C,D,E,F;
      /*S = bd1->nS+bd2->nS;
      Sx = bd1->nSx+bd2->nSx;
      Sy = bd1->nSy+bd2->nSy;
      Sxx = bd1->nSxx+bd2->nSxx;
      Sxy = bd1->nSxy+bd2->nSxy;
      Syy = bd1->nSyy+bd2->nSyy;*/
      Sxxx = bd1->nSxxx+bd2->nSxxx;
      Sxxy = bd1->nSxxy+bd2->nSxxy;
      Sxyy = bd1->nSxyy+bd2->nSxyy;
      Syyy = bd1->nSyyy+bd2->nSyyy;
      Sxxxx = bd1->nSxxxx+bd2->nSxxxx;
      Sxxyy = bd1->nSxxyy+bd2->nSxxyy;
      Syyyy = bd1->nSyyyy+bd2->nSyyyy;
      /*{ // shift
        dx = -(Sx+S/2)/S;
        dy = -(Sy+S/2)/S;
        Sxxxx += 4*dx*Sxxx+6*dx*dx*Sxx+4*dx*dx*dx*Sx+dx*dx*dx*dx*S;
        Sxxyy += 2*dx*Sxyy+dx*dx*Syy+2*dy*Sxxy+4*dx*dy*Sxy+2*dx*dx*dy*Sy+dy*dy*Sxx+2*dx*dy*dy*Sx+dx*dx*dy*dy*S;
        Syyyy += 4*dy*Syyy+6*dy*dy*Syy+4*dy*dy*dy*Sy+dy*dy*dy*dy*S;
        Sxxx  += 3*dx*Sxx+3*dx*dx*Sx+dx*dx*dx*S;
        Sxxy  += 2*dx*Sxy+dx*dx*Sy+dy*Sxx+2*dx*dy*Sx+dx*dx*dy*S;
        Sxyy  += dx*Syy+2*dy*Sxy+2*dx*dy*Sy+dy*dy*Sx+dx*dy*dy*S;
        Syyy  += 3*dy*Syy+3*dy*dy*Sy+dy*dy*dy*S;
        Sxx   += 2*dx*Sx+dx*dx*S;
        Sxy   += dx*Sy+dy*Sx+dx*dy*S;
        Syy   += 2*dy*Sy+dy*dy*S;
        Sx    += dx*S;
        Sy    += dy*S;
      }*/
      A = S*Sxx-Sx*Sx;
      B = S*Sxy-Sx*Sy;
      C = S*Syy-Sy*Sy;
      D = A*C-B*B;
      if (D==0)
        return dist1;
      E = S*(Sxxx+Sxyy)-Sx*(Sxx+Syy);
      F = S*(Sxxy+Syyy)-Sy*(Sxx+Syy);
      x0 = ((double)E*(double)C-(double)F*(double)B)/D;
      y0 = ((double)A*(double)F-(double)E*(double)B)/D;
      u  = (Sxx+Syy-x0*Sx-y0*Sy)/S;
      x0 /= 2;
      y0 /= 2;
      u  /= 2;
      r = 2*u+x0*x0+y0*y0;
      if (r<=0)
        return dist1; // error, should not happen indeed
      r = sqrt(r);
      err1 =   4*u*u*S
              +4*x0*x0*Sxx
              +4*y0*y0*Syy
              +8*u*x0*Sx
              +8*u*y0*Sy
              +8*x0*y0*Sxy
              -4*u*(Sxx+Syy)
              -4*x0*(Sxxx+Sxyy)
              -4*y0*(Sxxy+Syyy)
              +Sxxxx+Syyyy+2*Sxxyy;
//      x0 -= dx;
  //    y0 -= dy;
      err1 = sqrt(err1/S);
      dist2 = err1/r;
      dist2 = 10*dist2/sqrt((double)S);
    }
    return (dist1<dist2)?dist1:dist2;
}

MIA_RESULT_CODE IPL_ELL_GetApproximatingCircle(
  double *px0,
  double *py0,
  double *pr,
  const SBlobData4* pD)
{
  int64 A,B,C,D,E,F;
  double x0,y0,u,r;

    A = pD->nS*pD->nSxx-pD->nSx*pD->nSx;
    B = pD->nS*pD->nSxy-pD->nSx*pD->nSy;
    C = pD->nS*pD->nSyy-pD->nSy*pD->nSy;
    D = A*C-B*B;
    if (D==0)
      return ERR_GEN_BAD_DATA;
    E = pD->nS*(pD->nSxxx+pD->nSxyy)-pD->nSx*(pD->nSxx+pD->nSyy);
    F = pD->nS*(pD->nSxxy+pD->nSyyy)-pD->nSy*(pD->nSxx+pD->nSyy);
    x0 = ((double)E*(double)C-(double)F*(double)B)/D;
    y0 = ((double)A*(double)F-(double)E*(double)B)/D;
    u  = (pD->nSxx+pD->nSyy-x0*pD->nSx-y0*pD->nSy)/pD->nS;
    x0 /= 2;
    y0 /= 2;
    u  /= 2;
    r = 2*u+x0*x0+y0*y0;
    if (r<=0)
      return ERR_GEN_BAD_DATA;
    r = sqrt(r);
    *px0 = x0;
    *py0 = y0;
    *pr = r;
    return ERR_OK;
}

MIA_RESULT_CODE IPL_CONOB_Enumerate_CircularChunks_v1(
  SBlobInfo  **ppsBI,       // OUT: blob infos
  SBlobData4 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  const short* pgx,
  const short* pgy,
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg)       // IN:  for debug, set to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  SBlobInfo *psBI = NULL;   // blob infos
  SBlobData4 *psBD = NULL;  // blob data
  int *pnNeibList = NULL;   // neighbor indices
  int *pnHelperIdxs = NULL; // neighbor indices
  int x,y,cIdx,nNeibList_u,nNeibList_a,nPts;
  //char nama[FILENAME_MAX];
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
    // prepare helper array of indexes, count number of existing points
    if ((pnHelperIdxs = (int*)malloc(xs*ys*sizeof(pnHelperIdxs[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    }
    memset(pnHelperIdxs,-1,xs*ys*sizeof(pnHelperIdxs[0]));
    nPts = 0;
    for (cIdx=0;cIdx<xs*ys;cIdx++)
      if (im[cIdx])
        pnHelperIdxs[cIdx] = nPts++;
    if (!nPts)
      // no points, nothing to do
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    // allocate
    if ((psBI = (SBlobInfo*)malloc(sizeof(psBI[0])*nPts*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    }
    if ((psBD = (SBlobData4*)malloc(sizeof(psBD[0])*nPts*2))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
    }
    nNeibList_u = 0;
    if ((pnNeibList = (int*)malloc((nNeibList_a = nPts*8)*sizeof(pnNeibList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_Enumerate_StraightChunks_v1_exit;
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
        psBD[nBlob].nSxxx = x*x*x;
        psBD[nBlob].nSxxy = x*x*y;
        psBD[nBlob].nSxyy = x*y*y;
        psBD[nBlob].nSyyy = y*y*y;
        psBD[nBlob].nSxxxx = (int64)x*x*x*x;
        psBD[nBlob].nSxxyy = (int64)x*x*y*y;
        psBD[nBlob].nSyyyy = (int64)y*y*y*y;
        //psBD[nBlob].nSgx = pgx[cIdx];
        //psBD[nBlob].nSgy = pgy[cIdx];
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
        dDist = BlobDistance_Archunks_v2(psBD+cBlob,psBD+nnei);
        if ((dMinDistVal>dDist)&&(dDist<1.))
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
//      psBD[nBlob].nSgx = psBD[cBlob].nSgx+psBD[cNeib].nSgx;
  //    psBD[nBlob].nSgy = psBD[cBlob].nSgy+psBD[cNeib].nSgy;
      psBD[nBlob].nSxxx = psBD[cBlob].nSxxx+psBD[cNeib].nSxxx;
      psBD[nBlob].nSxxy = psBD[cBlob].nSxxy+psBD[cNeib].nSxxy;
      psBD[nBlob].nSxyy = psBD[cBlob].nSxyy+psBD[cNeib].nSxyy;
      psBD[nBlob].nSyyy = psBD[cBlob].nSyyy+psBD[cNeib].nSyyy;
      psBD[nBlob].nSxxxx = psBD[cBlob].nSxxxx+psBD[cNeib].nSxxxx;
      psBD[nBlob].nSxxyy = psBD[cBlob].nSxxyy+psBD[cNeib].nSxxyy;
      psBD[nBlob].nSyyyy = psBD[cBlob].nSyyyy+psBD[cNeib].nSyyyy;
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
if (nambeg)
{
  char nama[FILENAME_MAX];
  uint32 *dbgim,cObj,nObj;
  dbgim = (uint32*)malloc(xs*ys*4);
  memset(dbgim,0,xs*ys*4);
  for (cIdx=xs*ys-1;cIdx>=0;cIdx--)
    if (im[cIdx])
    {
      for (cIdx2=pnHelperIdxs[cIdx];psBI[cIdx2].nParent>=0;cIdx2=psBI[cIdx2].nParent);
      dbgim[cIdx] = ES_rand(cIdx2);
    }
  for (cObj=0;(int)cObj<nBlob;cObj++)
    if ((psBI[cObj].nParent<0)&&(psBD[cObj].nS>0))
    {
      SBlobData4* pD = psBD+cObj;
      double x0,y0,r;//,u;
      if (IPL_ELL_GetApproximatingCircle(&x0,&y0,&r,pD)!=ERR_OK)
        continue;
      if (r<10)
        continue;
      if (r>3*pD->nS)
        continue;
      r = 3*pD->nS/r;
  //          DBGL_DRAW_PlusInRGBA(dbgim,xs,ys,
  //          (int)(x0+.5),(int)(y0+.5),(int)(r+.5),ES_rand(cObj));
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
}
IPL_CONOB_Enumerate_StraightChunks_v1_exit:;
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
MIA_RESULT_CODE IPL_CONOB_GetArchunkData(
  SArchunkData2** ppArch, // OUT: list of extracted arc chunx
  int* pnArch,            // OUT: number of them
  const SBlobInfo* psBI,  // IN:  source blob infos
  const SBlobData4* psBD, // IN:  blob data
  int nBlob)              // IN:  number of source blobs (in blob heap)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnPointList=NULL;
  int *pnPointCounts=NULL;
  int nTopBlobCount;
  int *pnPointCoords=NULL;
  int nTopBlobPts;
  SArchunkData2* pArch=NULL;
  int nArch=0,cArch,cBlob;
  int *pnBackTable=NULL;

    // check arguments
    if ((ppArch==NULL)||(pnArch==NULL)||(psBI==NULL)||(psBD==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*ppArch!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    // collect point indices
    if ((res = IPL_CONOB_GetBlobPointIndices(
        &pnPointList,       // OUT: list of extracted points
        &pnPointCounts,     // OUT: triples of <start_position,length,src_blob>
        &nTopBlobCount,        // OUT: number of resulting blobs
        psBI,    // IN:  source blob infos
        nBlob))!=ERR_OK)            // IN:  number of source blobs (in blob heap)
      goto IPL_CONOB_GetArchunkData_exit;
    nArch = nTopBlobCount;
    // collect points
    nTopBlobPts = pnPointCounts[3*(nTopBlobCount-1)+0]+pnPointCounts[3*(nTopBlobCount-1)+1];
    if ((res = IPL_CONOB_GetBlobPointCoords(
        &pnPointCoords,     // OUT: coordinates as pairs <x,y>
        psBD,   // IN:  blob data
        pnPointList,   // IN:  list of point indexes
        pnPointCounts, // IN:  triples of <start_position,length,src_blob>
        nTopBlobPts))!=ERR_OK)
      goto IPL_CONOB_GetArchunkData_exit;
    // allocate back-table (resolving sBI/sBD source array number to Archunk number)
    if ((pnBackTable = (int*)malloc(nBlob*sizeof(pnBackTable[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_GetArchunkData_exit;
    }
    // fill in back-table
    for (cBlob=0;cBlob<nBlob;cBlob++)
      pnBackTable[cBlob] = -1;
    for (cArch=0;cArch<nArch;cArch++)
      pnBackTable[pnPointCounts[3*cArch+2]] = cArch;
    // allocate target structure
    if ((pArch = (SArchunkData2*)malloc(nTopBlobCount*sizeof(pArch[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_GetArchunkData_exit;
    }
    for (cArch=0;cArch<nArch;cArch++)
    {
      pArch[cArch].moments = psBD[pnPointCounts[3*cArch+2]];
      if ((pArch[cArch].pnPointList = (int*)malloc(
            pnPointCounts[3*cArch+1]*sizeof(int)*2))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto IPL_CONOB_GetArchunkData_exit;
      }
      memcpy( pArch[cArch].pnPointList,
              pnPointCoords+pnPointCounts[3*cArch+0]*2,
              pnPointCoords+pnPointCounts[3*cArch+1]*2*sizeof(int));
      if ((pArch[cArch].pnNeibList = (int*)malloc(
            psBI[pnPointCounts[3*cArch+2]]->nNeibCount*sizeof(int)))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto IPL_CONOB_GetArchunkData_exit;
      }
      memcpy( pArch[cArch].pnNeibList,
              );
    }
IPL_CONOB_GetArchunkData_exit:;
    if (pnPointList)
      free(pnPointList);
    if (pnPointCounts)
      free(pnPointCounts);
    if (pnPointCoords)
      free(pnPointCoords);
    if (pnBackTable)
      free(pnBackTable);
    if (res!=ERR_OK)
    {
      if (pArch)
      {
        for (cArch=0;cArch<nArch;cArch++)
        {
          if (pArch[cArch].pnPointList)
            free(pArch[cArch].pnPointList);
          if (pArch[cArch].pnNeibList)
            free(pArch[cArch].pnNeibList);
        }
        free(pArch);
      }
      pArch = NULL;
      nArch = 0;
    }
    *ppArch = pArch;
    pnArch = nArch;
    return res;
}
*/
