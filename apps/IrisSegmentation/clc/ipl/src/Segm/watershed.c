/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  28 December 2014                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  watershed processing                                   */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 129 // __FILENUM__TAG129

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "bpl.h"
#include "dbg.h"
#include "cgl.h"
#include "../ownipl.h"

typedef struct
{
  int64 M,Mx,My,Mxx,Mxy,Myy;//,Mc; // moments
  int num;                      // number of area
  int color;                    // color of area
  int nIdxOrdered;              // enumeration of adjacent areas by order
  int nIdxDarker;               // enumeration of adjacent areas by brightness
  int nMinColor;
  int criteria;
} SArea;

typedef struct 
{
  int nIdxNext;   // index of the next element
  int nValue;     // value
} SIntList;



/* enumerate areas segment-wise - overwhelming complexity
#define NONAREA 0x7fffffff

MIA_RESULT_CODE watershed_enumareas(
  SArea** ppArea,
  int *pnArea,
  SIntList** ppIntlist,
  int* pnIntList,
  const uint8* im,
  int xs,
  int ys,
  int connect8, // 0/1 : using 4-connectivity/8-connectivity
  int *sHash)
{
  MIA_RESULT_CODE res;
  unsigned char *ptr,*ptr_b,*ptr_e,*ptr_ee;
  int i,qhead,qtail,qdivmask,x_c,y_c,idx,len,area_u=0,area_a=0,ilist_u=0,ilist_a=0;
  int *queue=NULL;
  SArea *pArea=NULL,*pAr;
  int idx_scan;
  SIntList* pIntlist=NULL;
  int x,y,idx_b,idx_e,idx_b0,idx_e0,curcol;

    // check arguments
    if ( (ppArea==NULL)||(pnArea==NULL)||(im==NULL)||
         (ppIntlist==NULL)||(pnIntList==NULL))
      return ERR_GEN_NULLPOINTER;
    if ( (xs<3)||(ys<3) )
      return ERR_GEN_NO_DATA;
    // allocate sufficient space for queue
    // not less than (xs+ys)*4, not more than (xs+ys)*8
    //qdivmask = (1<<(MIA_IntLog2Hi(xs+ys)+3))-1;
    MIA_IntLog2Hi_Def(qdivmask,xs+ys);
    qdivmask = (1<<(qdivmask+3))-1;
    queue = (int*)malloc(sizeof(*queue)*2*(qdivmask+1));
    if (queue==NULL)
      return ERR_GEN_NOMEMORY;
    // set s-hash
    memset(sHash,-1,xs*ys*sizeof(sHash[0]));
//    for (i=xs-1;i>=0;i--)
  //    sHash[i] = sHash[xs*(ys-1)+i] = NONAREA;
    //for (i=xs*(ys-2);i>0;i-=xs)
      //sHash[i] = sHash[i+xs-1] = NONAREA;
    // queue is empty now
    qhead = qtail = 0;
    // shift image pointer to non-zero region
//    finish_y = im+xs+1;
    // default result is overflow
    res = ERR_IPL_QUEUE_OVERFLOW;
    // scan image - cycle along y
//    for (finish_x=finish_y+xs*(ys-3);finish_x>=finish_y;finish_x-=xs)
      // scan image - cycle along x
  //    for (ptr_scan=finish_x+xs-3;ptr_scan>=finish_x;ptr_scan--)
    //for (idx_scan=0;idx_scan<xs*ys;idx_scan++)
    for (y=0;y<ys;y++)
      for (x=0;x<ys;x++)
      {
        idx_scan = y*xs+x;
        //if (*ptr_scan)
        if (sHash[idx_scan]<0)
        { // unmarked zone found
          // ask for space
          if (area_a==area_u)
            if (pArea = (SArea*)realloc(pArea,(
                          area_a=(area_a)?(area_a*2):1024)*sizeof(*pArea))==NULL)
            {
              res = ERR_GEN_NOMEMORY;
              goto endfunction;
            }
          // init vars concerning this area
          pAr = pArea+area_u;
          MIA_memset(pAr,0,sizeof(*pAr));
          pAr->color = im[idx_scan];
          pAr->num = area_u;
          pAr->nIdxDarker = -1;
          pAr->nIdxOrdered = -1;
          // locate beginning and end of segment
//          idx_b0 = y*xs;
  //        for (idx_b=idx_scan-1;(idx_b>=idx_b0)&&(im[idx_b]==pAr->color)&&(sHash[idx_b]<0);idx_b--);
          idx_b = idx_scan;
          idx_e0 = y*xs+xs-1;
          for (idx_e=idx_scan+1;(idx_e<=idx_e0)&&(im[idx_e]==pAr->color)&&(sHash[idx_e]<0);idx_e++);
          // remember the points
          queue[qtail  ] = idx_b;      // index of start
          queue[qtail+1] = 1+idx_e-idx_b;   // length
          // shift tail
          qtail = (qtail+2)&qdivmask;
          // paint the segment in s-hash
//          memset(ptr_b,0,1+ptr_e-ptr_b);
//          for (;ptr_b<=ptr_e;*ptr_b++=0);
          for (;idx_b<=idx_e;sHash[idx_b++]=area_u);
          // process queue
          while (qtail!=qhead)
          {
            // take segment from the queue
            idx = queue[qhead  ];
            len = queue[qhead+1];
            x_c = idx%xs;
            y_c = idx/ys;
            // shift head
            qhead = (qhead+2)&qdivmask;
            // update object parameters
            pAr->M   += len;
            pAr->My  += len*y_c;
            pAr->Mx  += len*(2*x_c+len-1)/2;
            pAr->Myy += len*y_c*y_c;
            pAr->Mxy += len*y_c*(2*x_c+len-1)/2;
            pAr->Mxx += len*x_c*(x_c+len-1)+len*(len-1)*(2*len-1)/6;
            // scan upper line
            ptr    = im+idx-str-((connect8)?1:0);
            ptr_ee = im+idx-str+len-((!connect8)?1:0);
            for (;ptr<=ptr_ee;ptr++)
              if (*ptr)
              {
                for (ptr_b=ptr;ptr_b[-1];ptr_b--);
                for (ptr_e=ptr;ptr_e[ 1];ptr_e++);
                // remember the points
                queue[qtail  ] = (int)(ptr_b-im);      // index of start
                queue[qtail+1] = (int)(1+ptr_e-ptr_b);   // length
                // shift tail
                qtail = (qtail+2)&qdivmask;
                if (qtail==qhead)
                  goto endfunction;
                // clear the segment
//                memset(ptr_b,0,1+ptr_e-ptr_b);
                for(;ptr_b<=ptr_e;*ptr_b++=0);
                // set pointer
                ptr = ptr_e+1;
              }
            // scan lower line
            ptr    = im+idx+str-((connect8)?1:0);
            ptr_ee = im+idx+str+len-((!connect8)?1:0);
            for (;ptr<=ptr_ee;ptr++)
              if (*ptr)
              {
                for (ptr_b=ptr;ptr_b[-1];ptr_b--);
                for (ptr_e=ptr;ptr_e[ 1];ptr_e++);
                // remember the points
                queue[qtail  ] = (int)(ptr_b-im);      // index of start
                queue[qtail+1] = (int)(1+ptr_e-ptr_b);   // length
                // shift tail
                qtail = (qtail+2)&qdivmask;
                if (qtail==qhead)
                  goto endfunction;
                // clear the segment
//                memset(ptr_b,0,1+ptr_e-ptr_b);
                for(;ptr_b<=ptr_e;*ptr_b++=0);
                // set pointer
                ptr = ptr_e+1;
              }
          }
          // object enumerated
          area_u++;
        }
      }
    // everything finished OK
    res = ERR_OK;
endfunction:
    // free resources
    if (queue)
      free(queue);
    if (res!=ERR_OK)
    {
      if (pArea!=NULL)
        free(pArea);
      pArea = NULL;
      area_u = 0;
      if (pIntlist!=NULL)
        free(pIntlist);
      pIntlist = 0;
      ilist_u = 0;
    }
    // output
    *ppArea = pArea;
    *pnArea = area_u;
    *ppIntlist = pIntlist;
    *pnIntList = ilist_u;
    return res;
}
*/

// enumerate areas point-wise
staticFunc MIA_RESULT_CODE watershed_enumareas(
  SArea** ppArea,
  int *pnArea,
  SIntList** ppIntlist,
  int* pnIntList,
  const uint8* im,
  int xs,
  int ys,
  int connect8, // 0/1 : using 4-connectivity/8-connectivity
  void *pvBuf)   // buffer of size: 3*xs*ys*sizeof(int)
{
  MIA_RESULT_CODE res;
  int qhead,qtail,x_c,y_c,idx,area_u=0,ilist_u=0;
  int *queue,*aHash,*iHash;  // area-hash and image-hash
  SArea *pArea=NULL,*pAr;
  SIntList* pIntlist=NULL;
  int x,y,neibnum,nbden,idx_n;
  const int nbdiffs[16] = { 0,-1,   -1, 0,    1,0,    0,1, 
                           -1,-1,    1,-1,   -1,1,    1,1};

    // allocate
    aHash = (int*)pvBuf;
    iHash = aHash+xs*ys;
    queue = iHash+xs*ys;
    if ((pIntlist = (SIntList*)malloc(xs*ys*4*sizeof(*pIntlist)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto endfunction;
    }
    if ((pArea = (SArea*)malloc(xs*ys*sizeof(*pArea)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto endfunction;
    }
    // set hashes
    MIA_memset(aHash,-1,xs*ys*sizeof(aHash[0]));
    MIA_memset(iHash,-1,xs*ys*sizeof(iHash[0]));
    for (idx=0;idx<xs;idx++)
      iHash[idx] = iHash[idx+xs*(ys-1)] = -2;
    for (idx=1;idx<ys-1;idx++)
      iHash[idx*xs] = iHash[idx*xs+xs-1] = -2;
    // default result is overflow
    res = ERR_IPL_QUEUE_OVERFLOW;
    // scan image - cycle along y
    for (y=1;y<ys-1;y++)
      // scan image - cycle along x
      for (x=1;x<xs-1;x++)
      {
        idx = y*xs+x;
        if (iHash[idx]<0)
        { // unmarked zone found
          // init vars concerning this area
          pAr = pArea+area_u;
          MIA_memset(pAr,0,sizeof(*pAr));
          pAr->color = im[idx];
          pAr->num = area_u;
          pAr->nIdxDarker = -1;
          pAr->nIdxOrdered = -1;
          // init the queue
          qtail = qhead = 0;
          queue[qtail++] = idx;
          // paint the point in image-hash
          iHash[idx] = area_u;
          // process queue
          while (qtail!=qhead)
          {
            // take point from the queue
            idx = queue[qhead++];
            x_c = idx%xs;
            y_c = idx/xs;
            // update object parameters
            pAr->M   += 1;
            pAr->My  += y_c;
            pAr->Mx  += x_c;
            pAr->Myy += y_c*y_c;
            pAr->Mxy += y_c*x_c;
            pAr->Mxx += x_c*x_c;
            // scan neighborhood
            for (nbden=0;nbden<4*(1+connect8);nbden++)
            {
              idx_n = idx+nbdiffs[nbden*2+0]+xs*nbdiffs[nbden*2+1];
              // upper
              neibnum = iHash[idx_n];
              if (neibnum>=0)
              { // this neighbor point was already enumerated
                if ((neibnum!=area_u)&&(aHash[neibnum]!=area_u))
                { // it is not current area, and was not enumerated as connection
                  pIntlist[ilist_u].nIdxNext = pAr->nIdxOrdered;
                  pIntlist[ilist_u].nValue = neibnum;
                  pAr->nIdxOrdered = ilist_u;
                  ilist_u++;
                  aHash[neibnum] = area_u;
                }
              }
              else
                // this neighbor point was not enumerated so far
                if ((neibnum!=-2)&&(im[idx_n]==pAr->color))
                { // neighbor point is same color and not border: join the area
                  queue[qtail++] = idx_n;
                  iHash[idx_n] = area_u;
                }
            }
          }
          // object enumerated
          area_u++;
        }
      }
    // compact the memory blocks, since might be allocated with gross overhead
    if ((pArea = (SArea*)realloc(pArea,area_u*sizeof(*pArea)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto endfunction;
    }
    if ((pIntlist = (SIntList*)realloc(pIntlist,ilist_u*sizeof(*pIntlist)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto endfunction;
    }
    // everything finished OK
    res = ERR_OK;
endfunction:
    if (res!=ERR_OK)
    {
      if (pArea!=NULL)
        free(pArea);
      pArea = NULL;
      area_u = 0;
      if (pIntlist!=NULL)
        free(pIntlist);
      pIntlist = 0;
      ilist_u = 0;
    }
    // output
    *ppArea = pArea;
    *pnArea = area_u;
    *ppIntlist = pIntlist;
    *pnIntList = ilist_u;
    return res;
}

// few checks of result of area enumeration function (debug only)
staticFunc MIA_RESULT_CODE watershed_verifyareas(
  const SArea* pArea,
  int nArea,
  const SIntList* pIntlist,
  int nIntList,
  int xs,
  int ys)
{
  int M,idx,idx_n,cIntList=0;

    // calculate total mass
    M = 0;
    for (idx=0;idx<nArea;idx++)
      M += (int)(pArea[idx].M);
    if (M!=(xs-2)*(ys-2))
      return ERR_GEN_INTERNAL;
    // check that color of adjacent areas are different
    for (idx=0;idx<nArea;idx++)
    {
      for (idx_n=pArea[idx].nIdxOrdered;idx_n>=0;idx_n=pIntlist[idx_n].nIdxNext)
      {
        if (pArea[idx].color==pArea[pIntlist[idx_n].nValue].color)
          break;
        cIntList++;
      }
      if (idx_n>=0)
        break;
    }
    if (idx<nArea)
      return ERR_GEN_INTERNAL;
    // check if number of references is correct
    if (cIntList!=nIntList)
      return ERR_GEN_INTERNAL;
    return ERR_OK;
}

// re-orient references as to be {bright}->{dark}
staticFunc MIA_RESULT_CODE watershed_orientrefs(
  SArea* pArea,
  int nArea,
  SIntList* pIntlist)
{
  int idx_s;  // index of self in area array
  int idx_n;  // index of neighbor in area array
  int ref_l;  // index of reference in int-list
  MIA_RESULT_CODE res = ERR_OK;

    for (idx_s=0;idx_s<nArea;idx_s++)
      while ((ref_l = pArea[idx_s].nIdxOrdered)>=0)
      {
        // drop item 'ref_l' from ordered list of self
        pArea[idx_s].nIdxOrdered = pIntlist[ref_l].nIdxNext;
        // index of the neighbor
        idx_n = pIntlist[ref_l].nValue;
        // where should this reference go?
        if (pArea[idx_s].color>pArea[idx_n].color)
        { // keep orientation: add reference to the list of self
          pIntlist[ref_l].nIdxNext = pArea[idx_s].nIdxDarker;
          pArea[idx_s].nIdxDarker = ref_l;
        }
        else
        { // change orientation: add reference to the list of neighbor
          pIntlist[ref_l].nIdxNext = pArea[idx_n].nIdxDarker;
          pArea[idx_n].nIdxDarker = ref_l;
          pIntlist[ref_l].nValue = idx_s; // reference from neighbor to self
        }
      }
    return res;
}

// check if references are {bright}->{dark}
staticFunc MIA_RESULT_CODE watershed_verifyrefs(
  SArea* pArea,
  int nArea,
  SIntList* pIntlist,
  int nIntList)
{
  int idx,idx_n,cIntList=0;

    for (idx=0;idx<nArea;idx++)
    {
      // check that nothing is left in ordered list
      if (pArea[idx].nIdxOrdered>=0)
        break;
      // check that colors of adjacent areas are lower
      for (idx_n=pArea[idx].nIdxDarker;idx_n>=0;idx_n=pIntlist[idx_n].nIdxNext)
      {
        if (pArea[idx].color<=pArea[pIntlist[idx_n].nValue].color)
          break;
        cIntList++;
      }
      if (idx_n>=0)
        break;
    }
    if (idx<nArea)
      return ERR_GEN_INTERNAL;
    // check if number of references is correct
    if (cIntList!=nIntList)
      return ERR_GEN_INTERNAL;
    return ERR_OK;
}

#define cmp_area_brt(o1,o2) (o1->color<o2->color)
IMPLEMENT_HEAPSORT(watershed_sortareabrt,SArea*,cmp_area_brt)

// sort area by brightness increasing
staticFunc MIA_RESULT_CODE watershed_sortbrightness(
  SArea** ppArea,
  int nArea,
  SIntList* pIntlist,
  int nIntList,
  void* pvBuf)
{
  SArea *pArea,*pArea_new=NULL;
  SArea** arrefs;
  int cArea,*reindex,cIntList;
  MIA_RESULT_CODE res = ERR_OK;

    // allocate
    arrefs = (SArea**)pvBuf;
    reindex = (int*)(arrefs+nArea);
    if ((pArea_new = (SArea*)malloc(nArea*sizeof(*pArea_new)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto endfunction;
    }
    // create array of references
    pArea = *ppArea;
    for (cArea=0;cArea<nArea;cArea++)
      arrefs[cArea] = pArea+cArea;
    // sort
    watershed_sortareabrt(arrefs,nArea);
    // copy bodies in new order, fill in reindex array
    for (cArea=0;cArea<nArea;cArea++)
    {
      memcpy(pArea_new+cArea,arrefs[cArea],sizeof(*pArea_new));
      reindex[arrefs[cArea]->num] = cArea;
    }
    // free unnecessary
    free(pArea);
    pArea = pArea_new;
    pArea_new = NULL;
    // re-indexation
    for (cArea=0;cArea<nArea;cArea++)
      pArea[cArea].num = reindex[pArea[cArea].num];
    for (cIntList=0;cIntList<nIntList;cIntList++)
      pIntlist[cIntList].nValue = reindex[pIntlist[cIntList].nValue];
    // output
    *ppArea = pArea;
endfunction:;
    return res;
}

// enumerate objects basing on areas structure
staticFunc MIA_RESULT_CODE watershed_enumobjects(
  SArea** ppChamps,
  int* pnChamps,
  int* pnRuns,
  SArea* pArea,
  int nArea,
  SIntList* pIntlist,
  void* pvBuf)
{
  int cArea;
  SArea *pAr;
  int *oHash; // object hash
  int *idxList,idxdarker,nNeibObj,idxobject,cNeibObj;
  int criteria,maxcriteria,maxidx,tmpint,minbrt,nChamps=0,nRuns=0;
  SArea* pChamps;

    // allocate
    oHash = (int*)pvBuf;
    memset(oHash,-1,sizeof(*oHash)*nArea*2);
    idxList = oHash+nArea*2;
    pChamps = (SArea*)(idxList+nArea);
    // run along all areas (which are sorted by brightness increasing)
    for (cArea=0;cArea<nArea;cArea++)
    {
      pAr = pArea+cArea;
//if (cArea==95)
  //tmpint = 0;
      if (pAr->nIdxDarker<0)
      { // no darker neighbors - this area is an object itself
        // transform area to object (in place)
        //pAr->Mc = pAr->M*pAr->color;
        pAr->nMinColor = pAr->color;
        continue;
      }
      // there are darker neighbors, enumerate them
      for (nNeibObj=0,idxdarker=pAr->nIdxDarker;idxdarker>=0;idxdarker=pIntlist[idxdarker].nIdxNext)
      {
        // 'ordered' index contains uniting object
        for (idxobject=pIntlist[idxdarker].nValue;pArea[idxobject].nIdxOrdered>=0;)
        {
//          idxobject = pArea[idxobject].nIdxOrdered;
          tmpint = pArea[idxobject].nIdxOrdered;
          pArea[idxobject].nIdxOrdered = cArea;
          idxobject = tmpint;
          nRuns++;
        }
//        if (oHash[idxobject]!=cArea)
        if ((oHash[idxobject]!=cArea)&&(idxobject!=cArea))
        { // not in hash - new object - enter into the list
          idxList[nNeibObj++] = idxobject;
          oHash[idxobject] = cArea;
        }
      }
      // now 'idxlist' contains 'nNeibObj' high-level neighbors
      // process each of neighbors for current brightness
      maxcriteria = 0;
      maxidx = -1;
      minbrt = 255;
      for (cNeibObj=0;cNeibObj<nNeibObj;cNeibObj++)
      {
        // set reference to current object
        tmpint = idxList[cNeibObj];
        pArea[tmpint].nIdxOrdered = cArea;
        pArea[tmpint].criteria = -1;
        // update current object counters adding neighbor
        pAr->M   += pArea[tmpint].M;
        pAr->Mx  += pArea[tmpint].Mx;
        pAr->My  += pArea[tmpint].My;
        pAr->Mxx += pArea[tmpint].Mxx;
        pAr->Mxy += pArea[tmpint].Mxy;
        pAr->Myy += pArea[tmpint].Myy;
        // determine minimum brightness
        if (minbrt>pArea[tmpint].nMinColor)
          minbrt = pArea[tmpint].nMinColor;
//        if (minbrt==0)
  //        break;
        // calculate criteria for neighbor
        criteria = pAr->color-pArea[tmpint].nMinColor;
        if (criteria<70)
          continue;
        tmpint = (int)(pAr->M);
        if (tmpint>300)
          continue;
        if (tmpint<9)
          tmpint = 9;
//        criteria = (criteria*criteria*criteria+tmpint/2)/tmpint;
        criteria = (criteria*criteria+tmpint/2)/tmpint;
//        if (criteria<100)
  //        continue;
        pArea[idxList[cNeibObj]].criteria = criteria;
        if (maxcriteria<criteria)
        {
          maxcriteria = criteria;
          maxidx = cNeibObj;
        }
      }
      // set minimum brightness
      pAr->nMinColor = minbrt;
      // is there any champion?
      if (maxidx>=0)
      { // yes, store it
        maxidx = idxList[maxidx]; // number of object
        memcpy(pChamps+nChamps,pArea+maxidx,sizeof(*pChamps));
        nChamps++;
      }
    }
    *ppChamps = pChamps;
    *pnChamps = nChamps;
    *pnRuns = nRuns;
    return ERR_OK;
}

//#define cmp_area_xy(o1,o2) ((o1.Mx*o2.M < o2.Mx*o1.M)?1:(o1.My*o2.M < o2.My*o1.M))
#define cmp_area_xy(o1,o2) ((o1.Mx==o2.Mx)?(o1.My<o2.My):(o1.Mx<o2.Mx))
IMPLEMENT_HEAPSORT(watershed_sortarea_xy,SArea,cmp_area_xy)
//#define cmp_area_yx(o1,o2) ((o1.My*o2.M < o2.My*o1.M)?1:(o1.Mx*o2.M < o2.Mx*o1.M))
#define cmp_area_yx(o1,o2) ((o1.My==o2.My)?(o1.Mx<o2.Mx):(o1.My<o2.My))
IMPLEMENT_HEAPSORT(watershed_sortarea_yx,SArea,cmp_area_yx)

#define DELTCH 5

// filter 'champions' by proximity
staticFunc MIA_RESULT_CODE watershed_filterchamps(
  SArea* pChamps,
  int* pnChamps,
  void* pvBuf)
{
  int nChamps,cChamps,nDst,cComp,nBest,bx,by,xx,yy;
  SArea* pChamps2;

    pChamps2 = (SArea*)pvBuf;
    nChamps = *pnChamps;
    // set Mx,My to xc,yc
    for (cChamps=0;cChamps<nChamps;cChamps++)
    {
      pChamps[cChamps].Mx = (pChamps[cChamps].Mx+pChamps[cChamps].M/2)/pChamps[cChamps].M;
      pChamps[cChamps].My = (pChamps[cChamps].My+pChamps[cChamps].M/2)/pChamps[cChamps].M;
      pChamps[cChamps].M = 1;
    }
    // sort in yx order
    watershed_sortarea_yx(pChamps,nChamps);
    for (nDst=cChamps=0;cChamps<nChamps;)
    {
      // take 'cchamps' item, compare it with all following
      nBest = cChamps;
//      bx = (int)((pChamps[nBest].Mx+pChamps[nBest].M/2)/pChamps[nBest].M);
  //    by = (int)((pChamps[nBest].My+pChamps[nBest].M/2)/pChamps[nBest].M);
      bx = (int)(pChamps[nBest].Mx);
      by = (int)(pChamps[nBest].My);
      for (cComp=nBest+1;cComp<nChamps;cComp++)
      {
//        xx = (int)((pChamps[cComp].Mx+pChamps[cComp].M/2)/pChamps[cComp].M);
  //      yy = (int)((pChamps[cComp].My+pChamps[cComp].M/2)/pChamps[cComp].M);
        xx = (int)(pChamps[cComp].Mx);
        yy = (int)(pChamps[cComp].My);
        if ((xx-bx>DELTCH)||(xx-bx<-DELTCH)||(yy-by>DELTCH)||(yy-by<-DELTCH))
          break;
        if (pChamps[nBest].criteria<pChamps[cComp].criteria)
          nBest = cComp;
      }
      cChamps = cComp;
      memcpy(pChamps2+nDst,pChamps+nBest,sizeof(*pChamps2));
      nDst++;
    }
    // copy back
    memcpy(pChamps,pChamps2,nDst*sizeof(*pChamps));
    nChamps = nDst;
    // sort in xy order
    watershed_sortarea_xy(pChamps,nChamps);
    for (nDst=cChamps=0;cChamps<nChamps;)
    {
      // take 'cchamps' item, compare it with all following
      nBest = cChamps;
//      bx = (int)((pChamps[nBest].Mx+pChamps[nBest].M/2)/pChamps[nBest].M);
  //    by = (int)((pChamps[nBest].My+pChamps[nBest].M/2)/pChamps[nBest].M);
      bx = (int)(pChamps[nBest].Mx);
      by = (int)(pChamps[nBest].My);
      for (cComp=nBest+1;cComp<nChamps;cComp++)
      {
//        xx = (int)((pChamps[cComp].Mx+pChamps[cComp].M/2)/pChamps[cComp].M);
  //      yy = (int)((pChamps[cComp].My+pChamps[cComp].M/2)/pChamps[cComp].M);
        xx = (int)(pChamps[cComp].Mx);
        yy = (int)(pChamps[cComp].My);
        if ((xx-bx>DELTCH)||(xx-bx<-DELTCH)||(yy-by>DELTCH)||(yy-by<-DELTCH))
          break;
        if (pChamps[nBest].criteria<pChamps[cComp].criteria)
          nBest = cComp;
      }
      cChamps = cComp;
      memcpy(pChamps2+nDst,pChamps+nBest,sizeof(*pChamps2));
      nDst++;
    }
    // copy back
    memcpy(pChamps,pChamps2,nDst*sizeof(*pChamps));
    nChamps = nDst;
    // output
    *pnChamps = nChamps;
    return ERR_OK;
}

// watershed algorithm
MIA_RESULT_CODE watershed(
  const uint8* im,
  int xs,
  int ys,
  const char* nambeg,
  uint64* stix,
  int*   ntix)
{
  MIA_RESULT_CODE res = ERR_OK;
//  int* sHash=NULL;
  SArea* pArea=NULL;
  int nArea;
  SIntList* pIntlist=NULL;
  int nIntList,nChamps,nRuns;
  void* pvBuf = NULL;
  SArea* pChamps;
  int64 tix;
  //int64 Senumareas,Sverifyareas,Sorientrefs,Sverifyrefs,Ssortbrightness,Senumobjects;
  //int Nenumareas,Nverifyareas,Norientrefs,Nverifyrefs,Nsortbrightness,Nenumobjects;

    // check arguments
    if (im==NULL)
      return ERR_GEN_NULLPOINTER;
    if ( (xs<0)||(ys<0) )
      return ERR_GEN_NO_DATA;
    // allocate
    if ((pvBuf = malloc(3*xs*ys*sizeof(int)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto watershed_exit;
    }
    if ((stix!=NULL)&&(ntix!=NULL))
    { // DEBUG version
      // 1: enumerate areas
tix = CGL_GetTix();
      if ((res = watershed_enumareas(&pArea,&nArea,&pIntlist,&nIntList,
                    im,xs,ys,1,pvBuf))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[0] += tix;
ntix[0] ++;
      // verify enumeration
tix = CGL_GetTix();
      if ((res = watershed_verifyareas(pArea,nArea,pIntlist,nIntList,xs,ys))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[1] += tix;
ntix[1] ++;
      // 2: update directions
tix = CGL_GetTix();
      if ((res = watershed_orientrefs(pArea,nArea,pIntlist))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[2] += tix;
ntix[2] ++;
      // verify directions
tix = CGL_GetTix();
      if ((res = watershed_verifyrefs(pArea,nArea,pIntlist,nIntList))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[3] += tix;
ntix[3] ++;
      // 3: sort increasing brightness
tix = CGL_GetTix();
      if ((res = watershed_sortbrightness(&pArea,nArea,pIntlist,nIntList,pvBuf))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[4] += tix;
ntix[4] ++;
      // verify sorting
tix = CGL_GetTix();
      if ((res = watershed_verifyrefs(pArea,nArea,pIntlist,nIntList))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[3] += tix;
ntix[3] ++;
      // 4: create objects
tix = CGL_GetTix();
      if ((res = watershed_enumobjects(&pChamps,&nChamps,&nRuns,pArea,nArea,pIntlist,pvBuf))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[5] += tix;
ntix[5] ++;
      // 5: filter champions
tix = CGL_GetTix();
      if ((res = watershed_filterchamps(pChamps,&nChamps,pvBuf))!=ERR_OK)
        goto watershed_exit;
tix = CGL_GetTix()-tix;
stix[6] += tix;
ntix[6] ++;
    }
    else
    { // RELEASE version
      // 1: enumerate areas
      if ((res = watershed_enumareas(&pArea,&nArea,&pIntlist,&nIntList,
                    im,xs,ys,1,pvBuf))!=ERR_OK)
        goto watershed_exit;
      // 2: update directions
      if ((res = watershed_orientrefs(pArea,nArea,pIntlist))!=ERR_OK)
        goto watershed_exit;
      // 3: sort increasing brightness
      if ((res = watershed_sortbrightness(&pArea,nArea,pIntlist,nIntList,pvBuf))!=ERR_OK)
        goto watershed_exit;
      // 4: create objects
      if ((res = watershed_enumobjects(&pChamps,&nChamps,&nRuns,pArea,nArea,pIntlist,pvBuf))!=ERR_OK)
        goto watershed_exit;
      // 5: filter champions
      if ((res = watershed_filterchamps(pChamps,&nChamps,pvBuf))!=ERR_OK)
        goto watershed_exit;
    }
    // draw
if (nambeg)
{
uint8 *dbgim=NULL,*invim;
int i;
char nama[FILENAME_MAX];

  dbgim = (uint8*)malloc(xs*ys*4);
  invim = dbgim+xs*ys*3;
  for (i=0;i<xs*ys;i++)
    invim[i] = (uint8)(255-im[i]);
  DBGL_PUMP_Make3Bpp(dbgim,xs*3,invim,xs,xs,ys);
//  DBGL_PUMP_Make3Bpp(dbgim,xs*3,im,xs,xs,ys);
  for (i=0;i<nChamps;i++)
    DBGL_DRAW_PlusInRGB24(dbgim,xs,ys,
      (int)(pChamps[i].Mx/pChamps[i].M),(int)(pChamps[i].My/pChamps[i].M),5,0xff0000);
  sprintf(nama,"%d",nChamps);
  DBGL_TYPE_TextInRGB(dbgim,xs,ys,nama,0,0,0xff00ff,-1);
  strcpy(nama,nambeg);
  strcat(nama,"_zfin.bmp");
  DBGL_FILE_SaveRGB8Image(dbgim,xs,xs,ys,nama);
  free(dbgim);
}
watershed_exit:;
    if (pvBuf)
      free(pvBuf);
    if (pArea)
      free(pArea);
    if (pIntlist)
      free(pIntlist);
    return res;
}
