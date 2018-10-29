/*----------------------------------------------------------------------------*/
/*                                                                            */
/*      Created:  30 November 2000                                            */
/*      Modified: 13 July 2001 - adapted                                      */
/*      Modified: 13 April 2012 - removed ancient code                        */
/*      Modified: 19 September 2014 - added conobj optimized by segments      */
/*      Revision: 2.2.00                                                      */
/*      Purpose:  connected objects                                           */
/*      Authors:                                                              */
/*        Ivan Matveev                                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#define __FILENUM__ 199 // __FILENUM__TAG199

#include <math.h>
#include <memory.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include "ownipl.h"
#include "bpl.h"

#define INSPOINT(dx,dy,dY)                      \
  if (qptr[(dx)+(dY)])                          \
  {                                             \
    qptr[(dx)+(dY)] = 0;                        \
    queue[qtail] = ((y_c+(dy))<<16)+(x_c+(dx)); \
    qtail = (qtail+1)&qdivmask;                 \
    if (qtail==qhead)                           \
      goto endfunction;                         \
  }

// Non-optimized enumeration function (operating with pixels)
MIA_RESULT_CODE IPL_CONOB_EnumerateConnectedObjects_pixeled(
  SConObjPrimaryParams* pCOPP,  // OUT: buffer for object parameters
  int* pObjnum,        // IN/OUT: number of objects allocated/found
  unsigned char* im,            // input image (completely zeroed inside the function)
  int str,             // image stride in elements
  int xs,              // image size
  int ys,              //
  int MassThr_low,     // lower object mass threshold (-1 - disable)
  int MassThr_high,    // higher object mass threshold (-1 - disable)
  int connect8)        // if(TRUE) connect 8 neibours else - 4
{
  MIA_RESULT_CODE res;
  unsigned char *ptr,*finish_y,*finish_x,*qptr;
  int i,qhead,qtail,currentobject,currentcolor,qdivmask,x_c,y_c;
  int *queue=NULL;

    // check arguments
    if ( (pCOPP==NULL)||(pObjnum==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*pObjnum==0)
      return ERR_GEN_INVALID_PARAMETER;
    if ( (xs<3)||(ys<3) )
      return ERR_GEN_NO_DATA;
    // allocate sufficient space for queue
    // not less than (xs+ys)*4, not more than (xs+ys)*8
    //qdivmask = (1<<(MIA_IntLog2Hi(xs+ys)+3))-1;
    MIA_IntLog2Hi_Def(qdivmask,xs+ys);
    qdivmask = (1<<(qdivmask+3))-1;
    queue = (int*)malloc(sizeof(*queue)*(qdivmask+1));
    if (queue==NULL)
      return ERR_GEN_NOMEMORY;
    // bound image by zeroes
    memset(im,0,str*sizeof(im[0]));
    memset(im+str*(ys-1),0,str*sizeof(im[0]));
    for (i=str*(ys-2);i>0;i-=str)
      im[i] = im[i+xs-1] = 0;
    // queue is empty now
    qhead = qtail = 0;
    // no objects found yet
    currentobject = 0;
    // first color
    currentcolor = 0;
    // shift image pointer to non-zero region
    finish_y = im+str+1;
    // default result is overflow
    res = ERR_IPL_QUEUE_OVERFLOW;
    // scan image - cycle along y
    for (finish_x=finish_y+str*(ys-3);finish_x>=finish_y;finish_x-=str)
      // scan image - cycle along x
      for (ptr=finish_x+xs-3;ptr>=finish_x;ptr--)
        if (*ptr)
        { // new object found!
          // new color
          currentcolor++;
          // zero vars concerning this object
          memset(pCOPP,0,sizeof(pCOPP[0]));
          // remember the initial point
          pCOPP->c = currentcolor;
          ((unsigned short*)(queue+qtail))[0] = (unsigned short)(pCOPP->xb = (unsigned short)((ptr-im)%str));
          ((unsigned short*)(queue+qtail))[1] = (unsigned short)(pCOPP->yb = (unsigned short)((ptr-im)/str));
          // clear the point from image
          *ptr = 0;
          // shift tail
          qtail = (qtail+1)&qdivmask;
          // process queue
          while (qtail!=qhead)
          {
            // get a point from queue
            i = queue[qhead];
            x_c = i&0xffff;
            y_c = i>>16;
            // shift head
            qhead = (qhead+1)&qdivmask;
            // update object metrics
            pCOPP->M   += 1;
            pCOPP->MX  += x_c;
            pCOPP->MY  += y_c;
            pCOPP->MXX += x_c*x_c;
            pCOPP->MXY += x_c*y_c;
            pCOPP->MYY += y_c*y_c;
            // calculate image pointer
            qptr = im+y_c*str+x_c;
            // add neighbors to queue
            if(connect8)
            {
              INSPOINT(-1,-1,-str);    // 8-connectivity only
              INSPOINT(-1, 1, str);
              INSPOINT( 1,-1,-str);
              INSPOINT( 1, 1, str);
            }
            INSPOINT(-1, 0, 0);
            INSPOINT( 0,-1,-str);
            INSPOINT( 0, 1, str);
            INSPOINT( 1, 0, 0);
          }
          // object enumerated
          if (((pCOPP->M>=MassThr_low )||(MassThr_low ==-1))&&
              ((pCOPP->M<=MassThr_high)||(MassThr_high==-1)) )
          {
            pCOPP++;
            currentobject++;
          }
          if (currentobject>=(*pObjnum))
          {
            res = ERR_IPL_TOO_MANY_OBJECTS;
            goto endfunction;
          }
        }
    // everything finished OK
    res = ERR_OK;
endfunction:
    if (queue)
      free(queue);
    *pObjnum = currentobject;
    return res;
}

// enumeration function optimized by segments
MIA_RESULT_CODE IPL_CONOB_EnumerateConnectedObjects(
  SConObjPrimaryParams* pCOPP,  // OUT: buffer for object parameters
  int* pObjnum,         // IN/OUT: number of objects allocated/found
  unsigned char* im,    // input image (completely zeroed inside the function)
  int str,              // image stride in elements
  int xs,               // image size
  int ys,               //
  int MassThr_low,      // lower object mass threshold (-1 - disable)
  int MassThr_high,     // higher object mass threshold (-1 - disable)
  int connect8)         // if(TRUE) connect 8 neibors else - 4
{
  MIA_RESULT_CODE res;
  unsigned char *ptr,*finish_y,*finish_x,*ptr_b,*ptr_e,*ptr_ee,*ptr_scan;
  int i,qhead,qtail,currentobject,currentcolor,qdivmask,x_c,y_c,idx,len;
  int *queue=NULL;

    // check arguments
    if ( (pCOPP==NULL)||(pObjnum==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*pObjnum==0)
      return ERR_GEN_INVALID_PARAMETER;
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
    // bound image by zeroes
    memset(im,0,str*sizeof(im[0]));
    memset(im+str*(ys-1),0,str*sizeof(im[0]));
    for (i=str*(ys-2);i>0;i-=str)
      im[i] = im[i+xs-1] = 0;
    // queue is empty now
    qhead = qtail = 0;
    // no objects found yet
    currentobject = 0;
    // first color
    currentcolor = 0;
    // shift image pointer to non-zero region
    finish_y = im+str+1;
//    finish_y = im+str*(ys-1)+xs-1;
    // default result is overflow
    res = ERR_IPL_QUEUE_OVERFLOW;
    // scan image - cycle along y
    for (finish_x=finish_y+str*(ys-3);finish_x>=finish_y;finish_x-=str)
//    for (finish_x=im+str+xs-1;finish_x<finish_y;finish_x+=str)
      // scan image - cycle along x
      for (ptr_scan=finish_x+xs-3;ptr_scan>=finish_x;ptr_scan--)
//      for (ptr_scan=finish_x-xs+2;ptr_scan<finish_x;ptr_scan++)
        if (*ptr_scan)
        { // new object found!
          // new color
          currentcolor++;
          // zero vars concerning this object
          //memset(pCOPP,0,sizeof(pCOPP[0]));
          pCOPP->c = currentcolor;
          pCOPP->xb = (unsigned short)((ptr_scan-im)%str);
          pCOPP->yb = (unsigned short)((ptr_scan-im)/str);
          pCOPP->M = pCOPP->MX = pCOPP->MY = 0;
          pCOPP->MXX = pCOPP->MXY = pCOPP->MYY = 0;
          // locate beginning and end of segment
          for (ptr_b=ptr_scan;ptr_b[-1];ptr_b--);
          for (ptr_e=ptr_scan;ptr_e[ 1];ptr_e++);
          // remember the points
          queue[qtail  ] = (int)(ptr_b-im);      // index of start
          queue[qtail+1] = (int)(1+ptr_e-ptr_b);   // length
          // shift tail
          qtail = (qtail+2)&qdivmask;
          // clear the segment
//          memset(ptr_b,0,1+ptr_e-ptr_b);
          for (;ptr_b<=ptr_e;*ptr_b++=0);
          // process queue
          while (qtail!=qhead)
          {
            // take segment from the queue
            idx = queue[qhead  ];
            len = queue[qhead+1];
            x_c = idx%str;
            y_c = idx/str;
            // shift head
            qhead = (qhead+2)&qdivmask;
            // update object parameters
            pCOPP->M   += len;
            pCOPP->MY  += len*y_c;
            pCOPP->MX  += len*(2*x_c+len-1)/2;
            pCOPP->MYY += len*y_c*y_c;
            pCOPP->MXY += len*y_c*(2*x_c+len-1)/2;
            pCOPP->MXX += len*x_c*(x_c+len-1)+len*(len-1)*(2*len-1)/6;
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
          if (((pCOPP->M>=MassThr_low )||(MassThr_low ==-1))&&
              ((pCOPP->M<=MassThr_high)||(MassThr_high==-1)) )
          {
            pCOPP++;
            currentobject++;
            if (currentobject>=(*pObjnum))
            {
              res = ERR_IPL_TOO_MANY_OBJECTS;
              goto endfunction;
            }
          }
        }
    // everything finished OK
    res = ERR_OK;
endfunction:
    if (queue)
      free(queue);
    *pObjnum = currentobject;
    return res;
}

#define INSPOINT_FLOODFILL(dx,dy,dY,c)          \
  if (!qptr[(dx)+(dY)])                         \
  {                                             \
    qptr[(dx)+(dY)] = c;                        \
    queue[qtail] = ((y_c+(dy))<<16)+(x_c+(dx)); \
    qtail = (qtail+1)&qdivmask;                 \
    if (qtail==qhead)                           \
      goto endfunction;                         \
  }

// Non-optimized function
MIA_RESULT_CODE IPL_CONOB_FloodFill(
  unsigned char* im,              // input image (completely zeroed inside the function)
  unsigned int str,               // image stride in elements
  unsigned int xs,                // image size
  unsigned int ys,                //
  unsigned int xc,                // center to begin fill
  unsigned int yc,                //
  unsigned int connect8,          // if(TRUE) connect 8 neibours else - 4
  unsigned int color)             // filling color
{
  MIA_RESULT_CODE res=ERR_OK;
  unsigned char *ptr,*qptr;
  unsigned int i,qhead,qtail,qdivmask,x_c,y_c;
  unsigned int *queue=NULL;

    // check arguments
    if (im==NULL)
      return ERR_GEN_NULLPOINTER;
    if ( (xs<3)||(ys<3) )
      return ERR_GEN_NO_DATA;
    if (!(color&0xff))
      return ERR_GEN_INVALID_PARAMETER;
    // allocate sufficient space for queue
    // not less than (xs+ys)*8, not more than (xs+ys)*4
//    qdivmask = (1<<(MIA_IntLog2Hi(xs+ys)+3))-1;
    MIA_IntLog2Hi_Def(qdivmask,xs+ys);
    qdivmask = (1<<(qdivmask+3))-1;
    queue = (unsigned int*)malloc(sizeof(int)*(qdivmask+1));
    if (queue==NULL)
      return ERR_GEN_NOMEMORY;
    // bound image by max
    memset(im,255,str*sizeof(im[0]));
    memset(im+str*(ys-1),255,str*sizeof(im[0]));
    for (i=str*(ys-2);i>0;i-=str)
      im[i] = im[i+xs-1] = 255;
    // queue is empty now
    qhead = qtail = 0;
    // default result is overflow
    ptr = im+yc*str+xc;
    if (ptr[0])
      goto endfunction;
    // new object found!
    ((unsigned short*)(queue+qtail))[0] = (unsigned short)((ptr-im)%str);
    ((unsigned short*)(queue+qtail))[1] = (unsigned short)((ptr-im)/str);
    // fill the point
    *ptr = (unsigned char)color;
    // shift tail
    qtail = (qtail+1)&qdivmask;
    // process queue
    while (qtail!=qhead)
    {
      // get a point from queue
      i = queue[qhead];
      x_c = i&0xffff;
      y_c = i>>16;
      // shift head
      qhead = (qhead+1)&qdivmask;
      // calculate image pointer
      qptr = im+y_c*str+x_c;
      // add neighbors to queue
      if (connect8)
      {
        INSPOINT_FLOODFILL(-1,-1,-(int)str,(unsigned char)color);    // 8-connectivity only
        INSPOINT_FLOODFILL(-1, 1,      str,(unsigned char)color);
        INSPOINT_FLOODFILL( 1,-1,-(int)str,(unsigned char)color);
        INSPOINT_FLOODFILL( 1, 1,      str,(unsigned char)color);
      }
      INSPOINT_FLOODFILL(-1, 0,        0,(unsigned char)color);
      INSPOINT_FLOODFILL( 0,-1,-(int)str,(unsigned char)color);
      INSPOINT_FLOODFILL( 0, 1,      str,(unsigned char)color);
      INSPOINT_FLOODFILL( 1, 0,        0,(unsigned char)color);
    }
endfunction:;
    if (queue)
      free(queue);
    // bound image by min
    memset(im,0,str*sizeof(im[0]));
    memset(im+str*(ys-1),0,str*sizeof(im[0]));
    for (i=str*(ys-2);i>0;i-=str)
      im[i] = im[i+xs-1] = 0;
    return res;
}

#define INSPOINT_PAINTCLEAR(dx,dy,dY)           \
  if (paintptr[(dx)+(dY)])                      \
  {                                             \
    paintptr[(dx)+(dY)] = 0;                    \
    queue[qtail] = ((y_c+(dy))<<16)+(x_c+(dx)); \
    qtail = (qtail+1)&qdivmask;                 \
  }

// Non-optimized painting function
MIA_RESULT_CODE IPL_CONOB_PaintConnectedObjects(
  unsigned int* painting,        // OUT: painted buffer
  unsigned int dststr,           // stride in elements
  SConObjPrimaryParams* pCOPP, // OUT: buffer for object parameters
  unsigned int* pObjnum,         // IN:  size of this buffer
                                  // OUT: number of objects found
  unsigned char* im,              // input image (completely zeroed inside the function)
  unsigned int str,              // image stride in elements
  unsigned int xs,               // image size
  unsigned int ys,               //
  unsigned int MassThr,          // object mass threshold
  unsigned int connect8)         // if(TRUE) connect 8 neibours else - 4
{
  MIA_RESULT_CODE res;
  unsigned char *ptr,*finish_y,*finish_x,*qptr;
  unsigned int i,qhead,qtail,currentobject,currentcolor,qdivmask,x_c,y_c,*queue;
  unsigned int *paintptr;

    // check arguments
    if ( (pCOPP==NULL)||(pObjnum==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*pObjnum==0)
      return ERR_GEN_INVALID_PARAMETER;
    if ( (xs<3)||(ys<3) )
      return ERR_GEN_NO_DATA;
    // allocate sufficient space for queue
    // not less than (xs+ys)*4, not more than (xs+ys)*8
//    qdivmask = (1<<(MIA_IntLog2Hi(xs+ys)+3))-1;
    MIA_IntLog2Hi_Def(qdivmask,xs+ys);
    qdivmask = (1<<(qdivmask+3))-1;
    queue = (unsigned int*)malloc(sizeof(int)*(qdivmask+1));
    if (queue==NULL)
      return ERR_GEN_NOMEMORY;
    // clear painting
    memset(painting,0,ys*dststr*sizeof(painting[0]));
    // bound image by zeroes
    memset(im,0,str*sizeof(im[0]));
    memset(im+str*(ys-1),0,str*sizeof(im[0]));
    for (i=str*(ys-2);i>0;i-=str)
      im[i] = im[i+xs-1] = 0;
    // queue is empty now
    qhead = qtail = 0;
    // no objects found yet
    currentobject = 0;
    // first color
    currentcolor = 0;
    // shift image pointer to non-zero region
    finish_y = im+str+1;
    // default result is overflow
    res = ERR_IPL_QUEUE_OVERFLOW;
    // scan image - cycle along y
    for (finish_x=finish_y+str*(ys-3);finish_x>=finish_y;finish_x-=str)
      // scan image - cycle along x
      for (ptr=finish_x+xs-3;ptr>=finish_x;ptr--)
        if (*ptr)
        { // new object found!
          // new color
          currentcolor++;
          // zero vars concerning this object
          memset(pCOPP,0,sizeof(pCOPP[0]));
          // remember the initial point
          pCOPP->c = currentcolor;
          ((unsigned short*)(queue+qtail))[0] = (unsigned short)(pCOPP->xb = (unsigned short)((ptr-im)%str));
          ((unsigned short*)(queue+qtail))[1] = (unsigned short)(pCOPP->yb = (unsigned short)((ptr-im)/str));
          // clear the point from image
          *ptr = 0;
          // shift tail
          qtail = (qtail+1)&qdivmask;
          // process queue
          while (qtail!=qhead)
          {
            // get a point from queue
            i = queue[qhead];
            x_c = i&0xffff;
            y_c = i>>16;
            // shift head
            qhead = (qhead+1)&qdivmask;
            // update object metrics
            pCOPP->M   += 1;
            pCOPP->MX  += x_c;
            pCOPP->MY  += y_c;
            pCOPP->MXX += x_c*x_c;
            pCOPP->MXY += x_c*y_c;
            pCOPP->MYY += y_c*y_c;
            // calculate image pointer
            qptr = im+y_c*str+x_c;
            // paint
            painting[y_c*dststr+x_c] = currentcolor;
            // add neighbors to queue
            if(connect8)
            {
              INSPOINT(-1,-1,-(int)str);    // 8-connectivity only
              INSPOINT(-1, 1,       str);
              INSPOINT( 1,-1,-(int)str);
              INSPOINT( 1, 1,       str);
            }
            INSPOINT(-1, 0,         0);
            INSPOINT( 0,-1,-(int)str);
            INSPOINT( 0, 1,       str);
            INSPOINT( 1, 0,         0);
          }
          // object enumerated
          if (pCOPP->M>=(int)MassThr)
          {
            pCOPP++;
            currentobject++;
          }
          else
          { // object is too small. erase from painting
            // start the queue
            ((unsigned short*)(queue+qtail))[0] = (unsigned short)pCOPP->xb;
            ((unsigned short*)(queue+qtail))[1] = (unsigned short)pCOPP->yb;
            // clear the point from image
            painting[pCOPP->yb*dststr+pCOPP->xb] = 0;
            // shift tail
            qtail = (qtail+1)&qdivmask;
            // process queue
            while (qtail!=qhead)
            {
              // get a point from queue
              i = queue[qhead];
              x_c = i&0xffff;
              y_c = i>>16;
              // shift head
              qhead = (qhead+1)&qdivmask;
              // calculate image pointer
              paintptr = painting+y_c*str+x_c;
              // add neighbors to queue
              if(connect8)
              {
                INSPOINT_PAINTCLEAR(-1,-1,-(int)str);    // 8-connectivity only
                INSPOINT_PAINTCLEAR(-1, 1,       str);
                INSPOINT_PAINTCLEAR( 1,-1,-(int)str);
                INSPOINT_PAINTCLEAR( 1, 1,       str);
              }
              INSPOINT_PAINTCLEAR(-1, 0,         0);
              INSPOINT_PAINTCLEAR( 0,-1,-(int)str);
              INSPOINT_PAINTCLEAR( 0, 1,       str);
              INSPOINT_PAINTCLEAR( 1, 0,         0);
            }
          }
          if (currentobject>=(*pObjnum))
          {
            res = ERR_IPL_TOO_MANY_OBJECTS;
            goto endfunction;
          }
        }
    // everything finished OK
    res = ERR_OK;
endfunction:
    free(queue);
    *pObjnum = currentobject;
    return res;
}

#define INSPOINT_THR(dx,dy,dY,Thr)              \
  if (qptr[(dx)+(dY)]>=Thr)                     \
  {                                             \
    qptr[(dx)+(dY)] = 0;                        \
    queue[qtail] = ((y_c+(dy))<<16)+(x_c+(dx)); \
    qtail = (qtail+1)&qdivmask;                 \
    if (qtail==qhead)                           \
      goto endfunction;                         \
  }

// Non-optimized painting function
MIA_RESULT_CODE IPL_CONOB_PaintConnectedObjectsThresholded(
  unsigned int* painting,       // OUT: painted buffer
  unsigned int dststr,          // stride in elements
  SConObjPrimaryParams* pCOPP,  // OUT: buffer for object parameters
  unsigned int* pObjnum,        // IN/OUT:  size of buffer/number of objects found
  unsigned char* im,            // input image (completely zeroed inside the function)
  unsigned int str,             // image stride in elements
  unsigned int xs,              // image size
  unsigned int ys,              //
  unsigned int MassThr,         // object mass threshold
  unsigned int connect8,        // if(TRUE) connect 8 neibors else - 4
  uint8 ThrActivate,
  uint8 ThrContinue)
{
  MIA_RESULT_CODE res;
  unsigned char *ptr,*finish_y,*finish_x,*qptr;
  unsigned int i,qhead,qtail,currentobject,currentcolor,qdivmask,x_c,y_c,*queue;
  unsigned int *paintptr;

    // check arguments
    if ( (pCOPP==NULL)||(pObjnum==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*pObjnum==0)
      return ERR_GEN_INVALID_PARAMETER;
    if ( (xs<3)||(ys<3) )
      return ERR_GEN_NO_DATA;
    // allocate sufficient space for queue
    // not less than (xs+ys)*4, not more than (xs+ys)*8
//    qdivmask = (1<<(MIA_IntLog2Hi(xs+ys)+3))-1;
    MIA_IntLog2Hi_Def(qdivmask,xs+ys);
    qdivmask = (1<<(qdivmask+3))-1;
    queue = (unsigned int*)malloc(sizeof(int)*(qdivmask+1));
    if (queue==NULL)
      return ERR_GEN_NOMEMORY;
    // clear painting
    memset(painting,0,ys*dststr*sizeof(painting[0]));
    // bound image by zeroes
    memset(im,0,str*sizeof(im[0]));
    memset(im+str*(ys-1),0,str*sizeof(im[0]));
    for (i=str*(ys-2);i>0;i-=str)
      im[i] = im[i+xs-1] = 0;
    // queue is empty now
    qhead = qtail = 0;
    // no objects found yet
    currentobject = 0;
    // first color
    currentcolor = 0;
    // shift image pointer to non-zero region
    finish_y = im+str+1;
    // default result is overflow
    res = ERR_IPL_QUEUE_OVERFLOW;
    // scan image - cycle along y
    for (finish_x=finish_y+str*(ys-3);finish_x>=finish_y;finish_x-=str)
      // scan image - cycle along x
      for (ptr=finish_x+xs-3;ptr>=finish_x;ptr--)
        if (*ptr>=ThrActivate)
        { // new object found!
          // new color
          currentcolor++;
          // zero vars concerning this object
          memset(pCOPP,0,sizeof(pCOPP[0]));
          // remember the initial point
          pCOPP->c = currentcolor;
          ((unsigned short*)(queue+qtail))[0] = (unsigned short)(pCOPP->xb = (unsigned short)((ptr-im)%str));
          ((unsigned short*)(queue+qtail))[1] = (unsigned short)(pCOPP->yb = (unsigned short)((ptr-im)/str));
          // clear the point from image
          *ptr = 0;
          // shift tail
          qtail = (qtail+1)&qdivmask;
          // process queue
          while (qtail!=qhead)
          {
            // get a point from queue
            i = queue[qhead];
            x_c = i&0xffff;
            y_c = i>>16;
            // shift head
            qhead = (qhead+1)&qdivmask;
            // update object metrics
            pCOPP->M   += 1;
            pCOPP->MX  += x_c;
            pCOPP->MY  += y_c;
            pCOPP->MXX += x_c*x_c;
            pCOPP->MXY += x_c*y_c;
            pCOPP->MYY += y_c*y_c;
            // calculate image pointer
            qptr = im+y_c*str+x_c;
            // paint
            painting[y_c*dststr+x_c] = currentcolor;
            // add neighbors to queue
            if(connect8)
            { // 8-connectivity only
              INSPOINT_THR(-1,-1,-(int)str,ThrContinue);
              INSPOINT_THR(-1, 1,      str,ThrContinue);
              INSPOINT_THR( 1,-1,-(int)str,ThrContinue);
              INSPOINT_THR( 1, 1,      str,ThrContinue);
            }
            INSPOINT_THR(-1, 0,        0,ThrContinue);
            INSPOINT_THR( 0,-1,-(int)str,ThrContinue);
            INSPOINT_THR( 0, 1,      str,ThrContinue);
            INSPOINT_THR( 1, 0,        0,ThrContinue);
          }
          // object enumerated
          if (pCOPP->M>=(int)MassThr)
          {
            pCOPP++;
            currentobject++;
          }
          else
          { // object is too small, erase from painting
            // start the queue
            ((unsigned short*)(queue+qtail))[0] = (unsigned short)pCOPP->xb;
            ((unsigned short*)(queue+qtail))[1] = (unsigned short)pCOPP->yb;
            // clear the point from image
            painting[pCOPP->yb*dststr+pCOPP->xb] = 0;
            // shift tail
            qtail = (qtail+1)&qdivmask;
            // process queue
            while (qtail!=qhead)
            {
              // get a point from queue
              i = queue[qhead];
              x_c = i&0xffff;
              y_c = i>>16;
              // shift head
              qhead = (qhead+1)&qdivmask;
              // calculate image pointer
              paintptr = painting+y_c*str+x_c;
              // add neighbors to queue
              if (connect8)
              {
                INSPOINT_PAINTCLEAR(-1,-1,-(int)str);    // 8-connectivity only
                INSPOINT_PAINTCLEAR(-1, 1,       str);
                INSPOINT_PAINTCLEAR( 1,-1,-(int)str);
                INSPOINT_PAINTCLEAR( 1, 1,       str);
              }
              INSPOINT_PAINTCLEAR(-1, 0,         0);
              INSPOINT_PAINTCLEAR( 0,-1,-(int)str);
              INSPOINT_PAINTCLEAR( 0, 1,       str);
              INSPOINT_PAINTCLEAR( 1, 0,         0);
            }
          }
          if (currentobject>=(*pObjnum))
          {
            res = ERR_IPL_TOO_MANY_OBJECTS;
            goto endfunction;
          }
        }
    // everything finished OK
    res = ERR_OK;
endfunction:
    free(queue);
    *pObjnum = currentobject;
    return res;
}

// up-to-low recursion
MIA_RESULT_CODE IPL_CONOB_GetBlobPointIndices(
  int** ppnPointList,       // OUT: list of extracted points
  int** ppnPointCounts,     // OUT: triples of <start_position,length,src_blob>
  int*  pnBlobCount,        // OUT: number of resulting blobs
  const SBlobInfo *psBI,    // IN:  source blob infos
  int nBlobSrc)             // IN:  number of source blobs (in blob heap)
{
  MIA_RESULT_CODE res = ERR_OK;
  int *pnPointList=NULL,*wayback=NULL,*pnPointCounts=NULL;
  int nWayIdx,nCurIdx,nState,nPix,cObj,cTopObj=0,nTopObj,nListIdx;

    // check arguments
    if ((ppnPointList==NULL)||(ppnPointCounts==NULL)||(pnBlobCount==NULL)||(psBI==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppnPointList!=NULL)||(*ppnPointCounts!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    if (nBlobSrc<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // estimate number of pixels (objects without children, at the start of list)
    for (nPix=0;(psBI[nPix].nChild1<0)&&(nPix<nBlobSrc);nPix++);
    // estimate number of top objects (without parent)
    for (cObj=nTopObj=0;cObj<nBlobSrc;cObj++)
      nTopObj += ((psBI[cObj].nParent<0)?1:0);
    // allocate memory for point list
    if ((pnPointList = (int*)malloc(nPix*sizeof(pnPointList[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_ExtractBlobPoints_exit;
    }
    // allocate memory for point list counters
    if ((pnPointCounts = (int*)malloc(nTopObj*3*sizeof(pnPointCounts[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_ExtractBlobPoints_exit;
    }
    // allocate memory to track recursion
    if ((wayback = (int*)malloc(nPix*sizeof(wayback[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_CONOB_ExtractBlobPoints_exit;
    }
    // main cycle - by top-level objects
    nListIdx = 0;
    for (cTopObj=cObj=0;cObj<nBlobSrc;cObj++)
      if (psBI[cObj].nParent<0)
      { // top-level object, enumerate
        pnPointCounts[3*cTopObj+0] = nListIdx;
        // recursion
        nWayIdx = 0;
        nCurIdx = cObj;
        nState = 1;
        while ((nWayIdx>0)||(nState==1))
        {
          if (nState==1)
          { // came to curidx for the first time
            if (psBI[nCurIdx].nChild1>=0)
            { // we have a child, go to it
              wayback[nWayIdx++] = nCurIdx;
              nCurIdx = psBI[nCurIdx].nChild1;
              continue;
            }
            // no children, this is a final point
            pnPointList[nListIdx++] = nCurIdx;
            nState = 2;
            continue;
          }
          // we are here only is state=2 && wayback not empty
          nCurIdx = wayback[--nWayIdx];
          // go to child 2
          nCurIdx = psBI[nCurIdx].nChild2;
          nState = 1;
        }
        pnPointCounts[3*cTopObj+1] = nListIdx-pnPointCounts[3*cTopObj+0];
        pnPointCounts[3*cTopObj+2] = cObj;
        if (pnPointCounts[3*cTopObj+1]>=10)
          cTopObj++;
        else
          nListIdx -= pnPointCounts[3*cTopObj+1];
      }
    //
IPL_CONOB_ExtractBlobPoints_exit:;
    if (wayback)
      free(wayback);
    if (res!=ERR_OK)
    {
      if (pnPointList)
        free(pnPointList);
      pnPointList = NULL;
      if (pnPointCounts)
        free(pnPointCounts);
      pnPointCounts = NULL;
      cTopObj = 0;
    }
    *ppnPointList = pnPointList;
    *ppnPointCounts = pnPointCounts;
    *pnBlobCount = cTopObj;
    return res;
}

// collect points to handy array
MIA_RESULT_CODE IPL_CONOB_GetBlobPointCoords(
  int** ppnPointCoords,     // OUT: coordinates as pairs <x,y>
  const SBlobData4* psBD,   // IN:  blob data
  const int* pnPointList,   // IN:  list of point indexes
  const int* pnPointCounts, // IN:  triples of <start_position,length,src_blob>
  int nPoints)
{
  int *pnPointCoords=NULL;
  int cCoordIdx;

    pnPointCoords = (int*)malloc(2*sizeof(pnPointCoords[0])*nPoints);
    for (cCoordIdx=0;cCoordIdx<nPoints;cCoordIdx++)
    {
      pnPointCoords[cCoordIdx*2+0] = (int)(psBD[pnPointList[cCoordIdx]].nSx);
      pnPointCoords[cCoordIdx*2+1] = (int)(psBD[pnPointList[cCoordIdx]].nSy);
    }
    *ppnPointCoords = pnPointCoords;
    return ERR_OK;
}
