/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  4 August 2010                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  flash processing                                       */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 8 // __FILENUM__TAG8

#include "ipl.h"
/*
void IPL_SEGM_KillBlik_v3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  void* buf,
  int* len)
{
  unsigned char *FromLef,*FromRig,*FromTop,*FromBot,v;
  int i,j;

    FromLef = (unsigned char*)malloc(xs*ys*4);
    FromRig = FromLef+xs*ys;
    FromTop = FromRig+xs*ys;
    FromBot = FromTop+xs*ys;
    for (i=0;i<ys;i++)
    {
      v = 255;
      for (j=0;j<xs;j++)
      {
        if (v>src[i*xs+j])
          v = src[i*xs+j];
        else
          if (v<src[i*xs+j])
            v++;
        FromLef[i*xs+j] = v;
      }
    }
    for (i=0;i<ys;i++)
    {
      v = 255;
      for (j=xs-1;j>=0;j--)
      {
        if (v>src[i*xs+j])
          v = src[i*xs+j];
        else
          if (v<src[i*xs+j])
            v++;
        FromRig[i*xs+j] = v;
      }
    }
    for (j=0;j<xs;j++)
    {
      v = 255;
      for (i=0;i<ys;i++)
      {
        if (v>src[i*xs+j])
          v = src[i*xs+j];
        else
          if (v<src[i*xs+j])
            v++;
        FromTop[i*xs+j] = v;
      }
    }
    for (j=0;j<xs;j++)
    {
      v = 255;
      for (i=ys-1;i>=0;i--)
      {
        if (v>src[i*xs+j])
          v = src[i*xs+j];
        else
          if (v<src[i*xs+j])
            v++;
        FromBot[i*xs+j] = v;
      }
    }
    for (i=0;i<ys;i++)
      for (j=0;j<xs;j++)
      {
        v   = FromLef[i*xs+j];
        if (v<FromRig[i*xs+j])
          v = FromRig[i*xs+j];
        if (v<FromTop[i*xs+j])
          v = FromTop[i*xs+j];
        if (v<FromBot[i*xs+j])
          v = FromBot[i*xs+j];
        dst[i*xs+j] = v;
      }
    free(FromLef);
}
*/
MIA_RESULT_CODE IPL_SEGM_KillBlik(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  void* buf,
  int* len)
{
  unsigned char *LefRig,*str,v;
  const unsigned char *ptr;
  int i,j,sz;

    if ((xs&3)||(ys&3))
      return ERR_GEN_INVALID_PARAMETER;
    if (len==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    // calc size
    sz = xs;
    if (sz<ys)
      sz = ys;
    sz += xs*ys;
    if (buf==NULL)
    {
      *len = sz;
      return ERR_OK;
    }
    if (*len<sz)
    {
      *len = sz;
      return ERR_GEN_NOMEMORY;
    }
    *len = sz;
    // check more args
    if ((src==NULL)||(dst==NULL))
      return ERR_GEN_NULLPOINTER;
    // allocate
    LefRig = (unsigned char*)buf;
    // form top-bot image
    for (i=0;i<xs;i++)
    {
      ptr = src+i;
      str = LefRig+i;
      // slide from top
      v = 255;
      for (j=0;j<ys;j++)
      {
        if (v>*ptr)
          v = *ptr;
        else
          if (v<*ptr)
            v++;
        *str = v;
        ptr += xs;
        str += xs;
      }
      // slide from bot
      v = 255;
      for (j=0;j<ys;j++)
      {
        ptr -= xs;
        str -= xs;
        if (v>*ptr)
          v = *ptr;
        else
          if (v<*ptr)
            v++;
        if (*str<v)
          *str = v;
      }
    }
    // run lef-rig and form dst
    str = LefRig+xs*ys;
    for (i=0;i<ys;i++)
    {
      ptr = src+i*xs;
      // slide from left
      v = 255;
      for (j=0;j<xs;j++)
      {
        if (v>ptr[j])
          v = ptr[j];
        else
          if (v<ptr[j])
            v++;
        str[j] = v;
      }
      // slide from right
      v = 255;
      for (j=xs-1;j>=0;j--)
      {
        if (v>ptr[j])
          v = ptr[j];
        else
          if (v<ptr[j])
            v++;
        if (str[j]<v)
          str[j] = v;
      }
      // str is filled with lef-rig, compare with top-bot and save to dst
      ptr = LefRig+xs*i;
      for (j=0;j<xs;j++)
        dst[j] = (unsigned char)((str[j]>ptr[j])?str[j]:ptr[j]);
      dst += xs;
    }
    dst -= xs*ys;
    return ERR_OK;
}
