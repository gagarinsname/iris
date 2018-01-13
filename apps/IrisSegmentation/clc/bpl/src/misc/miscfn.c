#include <string.h>
#include "bpl.h"

MIA_RESULT_CODE BPL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize)       // IN:  string length
{
  char *str = "Version 1.4.0.5 from 20.09.2014";

    if ((int)strlen(str)>=pnStringSize)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    strcpy(pcVersionString,str);
    return ERR_OK;
}

//============= interpolation functions (for DM/SM sequences) ==================

MIA_RESULT_CODE BPL_InterpolateSequence_SimpleIntDM(
  int* pOutVal,       // OUT: buffer for output sequence
  int nOutLen,        // IN:  output sequence length
  const int* pInVal,  // IN:  input sequence values
  int nInLen,         // IN:  input sequence length
  int nBegDup)        // IN:  number of beginning duplicates
{
  MIA_RESULT_CODE res=ERR_OK;
  double x,frac;
  int intg,nOutSample;

    // check arguments
    if ((pInVal==NULL)||(pOutVal==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nOutLen<=0)||(nInLen<=0))
      return ERR_GEN_NO_DATA;
    if ((nBegDup<0)||(nBegDup>=nOutLen))
      return ERR_GEN_INVALID_PARAMETER;
    // make linear interpolation
    for (nOutSample=0;nOutSample<nOutLen;nOutSample++)
    {
      // endings should meet
      x = ((double)(nOutSample*(nInLen-1)))/(nOutLen-1);
      intg = (int)x;
      frac = x-intg;
      // mix the data
      if (frac!=0)
        pOutVal[nOutSample] = (int)((1.-frac)*pInVal[intg]+frac*pInVal[intg+1]+.5);
      else
        pOutVal[nOutSample] = pInVal[intg];
    }
    if (nBegDup)
    {
      // shift data
      memmove(pOutVal+nBegDup,pOutVal,sizeof(pOutVal[0])*(nOutLen-nBegDup));
      // duplicate
      for (nOutSample=0;nOutSample<nBegDup;nOutSample++)
        pOutVal[nOutSample] = pInVal[0];
    }
    return res;
}

MIA_RESULT_CODE BPL_InterpolateSequenceSimple(
  double* pOutVal,      // OUT: buffer for output sequence
  int nOutLen,          // IN:  output sequence length
  const double* pInVal, // IN:  input sequence values
  int nInLen)           // IN:  input sequence length
{
  MIA_RESULT_CODE res=ERR_OK;
  double x,frac;
  int intg,nOutSample;

    // check arguments
    if ((pInVal==NULL)||(pOutVal==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nOutLen<2)||(nInLen<2))
      return ERR_GEN_NO_DATA;
    // linear
    for (nOutSample=0;nOutSample<nOutLen;nOutSample++)
    {
      // endings should meet
      x = ((double)(nOutSample*(nInLen-1)))/(nOutLen-1);
      intg = (int)x;
      frac = x-intg;
      // mix the data
      if (frac!=0)
        pOutVal[nOutSample] = (1.-frac)*pInVal[intg]+frac*pInVal[intg+1];
      else
        pOutVal[nOutSample] = pInVal[intg];
    }
    return res;
}

// piece-wise linear interpolation
MIA_RESULT_CODE BPL_InterpolateSequence(
  double* pOutVal,      // OUT: buffer for output sequence
  int nOutLen,          // IN:  output sequence length
  const double* pInVal, // IN:  input sequence values
  const double* pInTim, // IN:  input sequence times
  int nInLen)           // IN:  input sequence length
{
  MIA_RESULT_CODE res=ERR_OK;
  double dInputInterval,dCurrentTime;
  int nOutSample,nInSample;
  int nLastSmaller;

    // check arguments
    if ((pInVal==NULL)||(pOutVal==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nOutLen<2)||(nInLen<2))
      return ERR_GEN_NO_DATA;
    // check income
    for (nInSample=0;nInSample<nInLen-1;nInSample++)
      if (pInTim[nInSample+1]<=pInTim[nInSample])
        return ERR_GEN_BAD_DATA;  // later mark has smaller time
    // get input interval
    dInputInterval = pInTim[nInLen-1]-pInTim[0];
    nLastSmaller = 0;
    // simply copy first
    pOutVal[0] = pInVal[0];
    // cycle for body
    for (nOutSample=1;nOutSample<nOutLen-1;nOutSample++)
    {
      // calculate input time matching output tick
      dCurrentTime = nOutSample*dInputInterval/(nOutLen-1);
      // find the tick having biggest value, smaller than current input time
      for (;pInTim[nLastSmaller]<=dCurrentTime;nLastSmaller++);
      nLastSmaller--;
      // calculate the share
      pOutVal[nOutSample] = 
        ((pInTim[nLastSmaller+1]-dCurrentTime)*pInVal[nLastSmaller  ]+
         (dCurrentTime-pInTim[nLastSmaller]  )*pInVal[nLastSmaller+1])/
         (pInTim[nLastSmaller+1]-pInTim[nLastSmaller]);
    }
    // simply copy last
    pOutVal[nOutLen-1] = pInVal[nInLen-1];
    return res;
}

// piece-wise linear interpolation
MIA_RESULT_CODE BPL_InterpolateSequence_IntDM(
  int* pOutVal,       // OUT: buffer for output sequence
  int nOutLen,        // IN:  output sequence length
  const int* pInVal,  // IN:  input sequence values
  const int* pInTim,  // IN:  input sequence times
  int nInLen,         // IN:  input sequence length
  int nShift)         // IN:  positive - add frames from start, negative - to end
{
  MIA_RESULT_CODE res=ERR_OK;
  double dInputInterval,dCurrentTime;
  int nOutSample,nInSample;
  int nLastSmaller;

    // check arguments
    if ((pInVal==NULL)||(pOutVal==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nOutLen<2)||(nInLen<2))
      return ERR_GEN_NO_DATA;
    if ((nShift>=nOutLen)||(-nShift>=nOutLen))
      return ERR_GEN_INVALID_PARAMETER;
    // check income
    for (nInSample=0;nInSample<nInLen-1;nInSample++)
      if (pInTim[nInSample+1]<=pInTim[nInSample])
        return ERR_GEN_BAD_DATA;  // later mark has smaller time
    // get input interval
    dInputInterval = pInTim[nInLen-1]-pInTim[0];
    nLastSmaller = 0;
    // simply copy first
    pOutVal[0] = pInVal[0];
    // cycle for body
    for (nOutSample=1;nOutSample<nOutLen-1;nOutSample++)
    {
      // calculate input time matching output tick
      dCurrentTime = nOutSample*dInputInterval/(nOutLen-1);
      // find the tick having biggest value, smaller than current input time
      for (;pInTim[nLastSmaller]<=dCurrentTime;nLastSmaller++);
      nLastSmaller--;
      // calculate the share
      pOutVal[nOutSample] = (int)(.5+
        ((pInTim[nLastSmaller+1]-dCurrentTime)*pInVal[nLastSmaller  ]+
         (dCurrentTime-pInTim[nLastSmaller]  )*pInVal[nLastSmaller+1])/
         (pInTim[nLastSmaller+1]-pInTim[nLastSmaller]) );
    }
    // simply copy last
    pOutVal[nOutLen-1] = pInVal[nInLen-1];
    // make move
    if (nShift>0)
    { // add from start
      memmove(pOutVal+nShift,pOutVal,(nOutLen-nShift)*sizeof(pOutVal[0]));
      for (nOutSample=0;nOutSample<nShift;nOutSample++)
        pOutVal[nOutSample] = pOutVal[nShift];
    }
    else
      if (nShift<0)
      { // add to end
        nShift = -nShift;
        memmove(pOutVal,pOutVal+nShift,(nOutLen-nShift)*sizeof(pOutVal[0]));
        for (nOutSample=nOutLen-nShift;nOutSample<nOutLen;nOutSample++)
          pOutVal[nOutSample] = pOutVal[nOutLen-nShift-1];
      }
    return res;
}

int32 BPL_ReadInt32Bigendian(
  const void* ptr)
{
  uint32 v = ((uint32*)ptr)[0];
  return (int32)( ((v>>24)&0x000000ff)+((v>> 8)&0x0000ff00)+
                  ((v<< 8)&0x00ff0000)+((v<<24)&0xff000000) );
}

void BPL_InvertEndianness_uint32(
  void* ptr,
  int num)
{
  int i;
  uint32 v;

    for (i=0;i<num;i++)
    {
      v = ((uint32*)ptr)[i];
      v = ((v>>24)&0x000000ff)+((v>> 8)&0x0000ff00)+
          ((v<< 8)&0x00ff0000)+((v<<24)&0xff000000);
      ((uint32*)ptr)[i] = v;
    }
}
