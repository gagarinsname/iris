/*----------------------------------------------------------------------------*/
/*                                                                            */
/*      Created:  31 October 2007                                             */
/*      Revision: 1.1.00                                                      */
/*      Purpose:  clusterisation methods                                      */
/*      Authors:                                                              */
/*        Ivan Matveev                                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#define __FILENUM__ 197 // __FILENUM__TAG197

#include "ipl.h"
#include "bpl.h"

typedef struct
{
  int nClusterNumber;
  int nParent;
  int nChild1;
  int nChild2;
  unsigned int nMatchValue;
} SClusterTreeNode;

// sorting of ClusterTreeNode by MatchValue increasing
#define cmp_ctn_incr(s1,s2) (\
  ((s1).nMatchValue!=(s2).nMatchValue)? \
  ((s1).nMatchValue <(s2).nMatchValue): \
 (((s1).nChild1    !=(s2).nChild1)? \
  ((s1).nChild1     <(s2).nChild1): \
  ((s1).nChild2     <(s2).nChild2)  )   )
staticFunc IMPLEMENT_HEAPSORT(HS_ClusterTreeNode_MatchValue_inc,SClusterTreeNode,cmp_ctn_incr)

#define UINTINFI ((unsigned int)(-1))

staticFunc void ownIPL_PaintTree(
  SClusterTreeNode* psTree,
  int idx,
  int color)
{
    psTree[idx].nClusterNumber = color;
    if (psTree[idx].nChild1>=0)
      ownIPL_PaintTree(psTree,psTree[idx].nChild1,color);
    if (psTree[idx].nChild2>=0)
      ownIPL_PaintTree(psTree,psTree[idx].nChild2,color);
}

// clusterize elements with minimal-distance joining rule
MIA_RESULT_CODE IPL_SEGM_ClusteriseByMinDist(
  int* anClusterNumbers,    // OUT: cluster numbers
  int* pnClasses,           // IN/OUT: number of clusters required or
                            //         zero to judge itself/obtained
  int nElements,            // IN:  number of elements
  const void** pavFeatures, // IN:  features to be compared
  unsigned int (*pfnCompareFeatures)(const void*,const void*), 
                            // IN:  comparator. Must return -1 for infinity
  unsigned char* pucBuffer, // IN:  external memory buffer
  int* pnBufsize)           // IN/OUT: buffer size allocated/used
{
  int i,j,idx1,idx2,idx0;
  int sz,nMatchCount,nClusterCount,color,nClusterFinalCount;
  unsigned int max_dist,min_dist;
  SClusterTreeNode *psTree,*psMatches;

    // check arguments
    if (pnBufsize==NULL)
      return ERR_GEN_NULLPOINTER;
    if (nElements<=2)
      return ERR_GEN_INVALID_PARAMETER;
    // calculate size
    sz = ((nElements*(nElements-1))/2)*sizeof(psMatches[0])+
         2*nElements*sizeof(psTree[0])+16;
    if (pucBuffer==NULL)
    { // mem size query
      *pnBufsize = sz;
      return ERR_OK;
    }
    if (*pnBufsize<sz)
    {
      *pnBufsize = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    // check more arguments
    if ((anClusterNumbers==NULL)||(pnClasses==NULL)||
        (pavFeatures==NULL)||(pfnCompareFeatures==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*pnClasses>0)
      if ((*pnClasses>=nElements)||(*pnClasses<=1))
        return ERR_GEN_INVALID_PARAMETER;
    // allocate with 16-byte alignment
    psTree = (SClusterTreeNode*)((((ptr_t)pucBuffer)+15)&(~15));
    psMatches = psTree+2*nElements;
    // fill the table
    for (nMatchCount=i=0;i<nElements;i++)
      for (j=i+1;j<nElements;j++)
      {
        psMatches[nMatchCount].nChild1 = i;
        psMatches[nMatchCount].nChild2 = j;
        psMatches[nMatchCount].nMatchValue = pfnCompareFeatures(pavFeatures[i],pavFeatures[j]);
        nMatchCount++;
      }
    // sort the matching table
    HS_ClusterTreeNode_MatchValue_inc(psMatches,nMatchCount);
    // init three leaves
    for (i=0;i<nElements;i++)
    {
      psTree[i].nClusterNumber = i;
      psTree[i].nParent = -1;
      psTree[i].nChild1 = -1;
      psTree[i].nChild2 = -1;
      psTree[i].nMatchValue = (unsigned int)(-1);
    }
    // build the tree
    nClusterCount = nElements;
    nClusterFinalCount = 1;
    for (i=0;(i<nMatchCount)&&(nClusterCount>1);i++)
    {
      if (psMatches[i].nMatchValue==UINTINFI)
        // All valid is merged. Fix it and compare merging just for uniformity.
        nClusterFinalCount = nClusterCount;
      if (psTree[psMatches[i].nChild1].nClusterNumber==
          psTree[psMatches[i].nChild2].nClusterNumber)
        continue;   // compared items belong to same cluster
      // Items from different clusters found. Merge these clusters.
      // find roots of two clusters being merged
      for (idx1=psMatches[i].nChild1;psTree[idx1].nParent>=0;idx1=psTree[idx1].nParent);
      for (idx2=psMatches[i].nChild2;psTree[idx2].nParent>=0;idx2=psTree[idx2].nParent);
      // color of new node a a color of first branch
      color = psTree[idx1].nClusterNumber;
      // number of new node
      idx0 = nElements+nElements-nClusterCount;
      // create new node
      psTree[idx0].nClusterNumber = color;
      psTree[idx0].nParent = -1;
      psTree[idx0].nChild1 = idx1;
      psTree[idx0].nChild2 = idx2;
      psTree[idx0].nMatchValue = psMatches[i].nMatchValue;
      // link old nodes to new as children to parent
      psTree[idx1].nParent = idx0;
      psTree[idx2].nParent = idx0;
      // fill clusternumber in second cluster with value of first one
      ownIPL_PaintTree(psTree,idx2,color);
      // two clusters merged into one
      nClusterCount--;
    }
    nClusterCount = nClusterFinalCount;
    // determine number of clusters to be split to
    if ((nClusterCount==1)&&(*pnClasses>1))
      // all objects can be merged to one cluster, final number is given
      nClusterFinalCount = *pnClasses;
    else
      if ((nClusterCount>1)&&(*pnClasses<=1))
        // clusters cannot be merged and final number not given
        nClusterFinalCount = nClusterCount;
      else
        if ((nClusterCount>1)&&(*pnClasses>1)&&(*pnClasses<=nClusterCount))
          nClusterFinalCount = nClusterCount;
        else
        { // find number of clusters by a distance gap
          idx0 = max_dist = 0;
          min_dist = UINTINFI;
          for (i=0;i<nElements-nClusterCount-1;i++)
          {
            if ((int)max_dist<=(int)(psTree[nElements+i+1].nMatchValue-psTree[nElements+i].nMatchValue))
            {
              max_dist = psTree[nElements+i+1].nMatchValue-psTree[nElements+i].nMatchValue;
              idx0 = i;
            }
            if (min_dist>(psTree[nElements+i+1].nMatchValue-psTree[nElements+i].nMatchValue))
              min_dist = psTree[nElements+i+1].nMatchValue-psTree[nElements+i].nMatchValue;
          }
          if (min_dist==max_dist)
            nClusterFinalCount = 1;
          else
            nClusterFinalCount = nElements-idx0-1;
        }
    // paint tree to zero color
    ownIPL_PaintTree(psTree,nElements+nElements-1,0);
    // split by painting in different colors
    for (i=1;i<nClusterFinalCount;i++)
      ownIPL_PaintTree(psTree,psTree[nElements+nElements-i-1].nChild2,i);
    // output data
    for (i=0;i<nElements;i++)
      anClusterNumbers[i] = psTree[i].nClusterNumber;
    *pnClasses = nClusterFinalCount;
    return ERR_OK;
}
