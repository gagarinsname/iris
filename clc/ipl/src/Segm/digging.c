#define __FILENUM__ 6 // __FILENUM__TAG6

// segmentation for dsp
#include "ipl.h"

//== find max digging stuff
#define MASK_CLIMB_STEP_MAX(delta)  \
  if (maxptr[delta]>*maxptr)              \
  {                                       \
    maxptr += (delta);                    \
    continue;                             \
  }

#define ADD_TO_QUEUE_EDGE_8_MAX(delta,TYPE)                     \
  if (maxptr[delta]==max_value)                                  \
    {                                                            \
      *(queue8_tail++)= (ptr_t)(maxptr+(delta));          \
      maxptr[delta] = 0;                                         \
      if (queue8_tail>=queue8_finish)                            \
        return ENHANCED_ERROR_CODE(0x202);                       \
    }                                                            \
  else if (maxptr[delta]>max_value)                              \
    {                                                            \
      while (queue8_tail>queue8_start)                           \
      {                                                          \
        addr=(TYPE*) *(--queue8_tail);                  \
        *addr=max_value;                                         \
      };                                                         \
      maxptr+=(delta);                                           \
      queue8_head  =queue8_start;                                \
      queue8_tail  =queue8_start;                                \
      M = MX = MY = 0;                                           \
      max_value = *maxptr;                                       \
      *queue8_tail++ = (ptr_t) (maxptr);                  \
      *maxptr = 0;                                               \
      continue;                                                  \
    }

#define ADD_TO_QUEUE_EDGE_INTERNAL_MAX(number)             \
  array_val=*addr;                                               \
  if (array_val)                                                 \
    {                                                            \
      if (array_val<maxthreshold)                                \
        *(addr) = 0;                                             \
      else if (array_val<=value)                                 \
        {                                                        \
          *queue_tail++=number;                                  \
          *queue_tail++=array_val;                               \
          *queue_tail++=(ptr_t)addr;                    \
          if (queue_tail>=queue_finish)                          \
            queue_tail=queue_start;                              \
          *(addr) = 0;                                           \
        };                                                       \
    }

#define ADD_TO_QUEUE_EDGE_0_MAX                            \
  addr=maxptr-1;                                                 \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(0)

#define ADD_TO_QUEUE_EDGE_1_MAX                            \
  addr=maxptr-w-1;                                               \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(1)

#define ADD_TO_QUEUE_EDGE_2_MAX                            \
  addr=maxptr-w;                                                 \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(2)

#define ADD_TO_QUEUE_EDGE_3_MAX                            \
  addr=maxptr-w+1;                                               \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(3)

#define ADD_TO_QUEUE_EDGE_4_MAX                            \
  addr=maxptr+1;                                                 \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(4)

#define ADD_TO_QUEUE_EDGE_5_MAX                            \
  addr=maxptr+w+1;                                               \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(5)

#define ADD_TO_QUEUE_EDGE_6_MAX                            \
  addr=maxptr+w;                                                 \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(6)

#define ADD_TO_QUEUE_EDGE_7_MAX                            \
  addr=maxptr+w-1;                                               \
  ADD_TO_QUEUE_EDGE_INTERNAL_MAX(7)

MIA_RESULT_CODE IPL_SEGM_FindMaximaDigging_uint8(
  SValuePosition* maxima,       // OUT:    destination array
  int* size,           // IN/OUT: size of array in elements/number of extremums found
  uint8* array,         // IN:     image to be searched for maxima (zeroed on exit)
  int w,                        // IN:     image size
  int h,                        //
  uint8 maxthreshold,   // IN:     threshold of max
  void* extbuftmp)              // IN:     temporary buffer
{
  int queue_size,queue8_size;
  ptr_t *queue,*queue8;
  int  y,maxima_counter,edge,M,MX,MY;
  ptr_t *queue8_start,*queue8_finish,*queue8_tail,*queue8_head;
  ptr_t *queue_start,*queue_finish,*queue_tail,*queue_head;
  uint8 array_val,value,max_value;
  uint8 *ptr,*maxptr,*finish,*addr;

    queue_size = queue8_size = w*h/(2*sizeof(ptr_t));
    queue = (ptr_t*)extbuftmp;
    queue8 = queue+queue_size;
    queue8_start  = queue8;
    queue8_finish = queue8_start+queue8_size;
    queue8_head   = queue8_start;
    queue8_tail   = queue8_start;
    queue_start  = queue;
    queue_finish = queue_start+queue_size;
    queue_head   = queue_start;
    queue_tail   = queue_start;

    // clear borders
    memset(array,0,w*sizeof(array[0]));
    memset(array+w*(h-1),0,w*sizeof(array[0]));
    for (y=(h-1)*w;y>0;y-=w)
      array[y] = 0;
    for (y=h*w-1;y>w;y-=w)
      array[y] = 0;
    // init
    maxima_counter = 0;
    ptr = array+w+1;
    finish = array+w*(h-1)-1;
    // run along array
    for (;;)
    {
      while ((ptr<finish)&&(*ptr==0))
        ptr++;
      if (ptr==finish)
      { // end up here
        // output element count
        *size = maxima_counter;
        // sort by value descending
        return ERR_OK;
      }
      // we've found a non-zero point!
      // is it a big enough?
      if (*ptr<maxthreshold)
      { // no
        *ptr = 0;
        continue;
      }
      // we've found good point!!!
      maxptr = ptr;
      // gradient climbing for max (8-connectivity)
      for (;;)
      {
        MASK_CLIMB_STEP_MAX(1);
        MASK_CLIMB_STEP_MAX(w);
        MASK_CLIMB_STEP_MAX(w+1);
        MASK_CLIMB_STEP_MAX(w-1);
        MASK_CLIMB_STEP_MAX(-1);
        MASK_CLIMB_STEP_MAX(-w);
        MASK_CLIMB_STEP_MAX(-w+1);
        MASK_CLIMB_STEP_MAX(-w-1);
        break;
      }
      // we have found a local maxima point
      // however it can contain several points - enumerate them.
      queue8_head = queue8_start;
      queue8_tail = queue8_start;
      M = MX = MY = 0;
      max_value = *maxptr;
      *queue8_tail++ = (ptr_t)(maxptr);
      *maxptr = 0;
      do
      {
        maxptr = (uint8*) *queue8_head++;
        M++;
        MX += (int)((maxptr-array)%w);
        MY += (int)((maxptr-array)/w);
        if (queue8_head==queue8_finish)
          return ENHANCED_ERROR_CODE(0x202);
        ADD_TO_QUEUE_EDGE_8_MAX(1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(w,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(w+1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(w-1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-w,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-w+1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-w-1,uint8);
      }
      while (queue8_tail!=queue8_head);
      // all the points of maxima are now in queue elements [0..qtail=qhead)
      // output to maxima array
      maxima[maxima_counter].value = max_value;
      maxima[maxima_counter].x = (MX+M/2)/M;
      maxima[maxima_counter].y = (MY+M/2)/M;
      maxima_counter++;
      if (maxima_counter==*size)
        return ERR_GEN_INSUFFICIENT_BUFFER;
      // set qhead=0 and start digging
      queue8_head=queue8_start;
      queue_head  =queue_start;
      queue_tail  =queue_start;
      do
      {
        value = max_value;
        maxptr = (uint8*) *queue8_head++;
        ADD_TO_QUEUE_EDGE_0_MAX;
        ADD_TO_QUEUE_EDGE_1_MAX;
        ADD_TO_QUEUE_EDGE_2_MAX;
        ADD_TO_QUEUE_EDGE_3_MAX;
        ADD_TO_QUEUE_EDGE_4_MAX;
        ADD_TO_QUEUE_EDGE_5_MAX;
        ADD_TO_QUEUE_EDGE_6_MAX;
        ADD_TO_QUEUE_EDGE_7_MAX;
        while (queue_head!=queue_tail)
        {
          edge  =  (uint8)(*queue_head++);
          value =  (uint8)(*queue_head++);
          maxptr = (uint8*) *queue_head++;
          if (queue_head>=queue_finish)
            queue_head=queue_start;
          switch (edge)
          {
            case 0:
Case_Entry_0: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr=maxptr-1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_0;
                  }
              }
              break;
            case 1:
Case_Entry_1: ADD_TO_QUEUE_EDGE_0_MAX;
              ADD_TO_QUEUE_EDGE_2_MAX;
              ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr-w-1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_1;
                  }
              }
              break;
            case 2:
Case_Entry_2: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_3_MAX;
              addr=maxptr-w;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_2;
                  }
              }
              break;
            case 3:
Case_Entry_3: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_2_MAX;
              ADD_TO_QUEUE_EDGE_4_MAX;
              ADD_TO_QUEUE_EDGE_5_MAX;
              addr=maxptr-w+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_3;
                  }
              }
              break;
            case 4:
Case_Entry_4: ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_5_MAX;
              addr = maxptr+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_4;
                  }
              }
              break;
            case 5:
Case_Entry_5: ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_4_MAX;
              ADD_TO_QUEUE_EDGE_6_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr+w+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_5;
                  }
              }
              break;
            case 6:
Case_Entry_6: ADD_TO_QUEUE_EDGE_5_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr+w;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_6;
                  }
              }
              break;
            case 7:
Case_Entry_7: ADD_TO_QUEUE_EDGE_5_MAX;
              ADD_TO_QUEUE_EDGE_6_MAX;
              ADD_TO_QUEUE_EDGE_0_MAX;
              ADD_TO_QUEUE_EDGE_1_MAX;
              addr = maxptr+w-1;
              array_val = *addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_7;
                  }
              }
              break;
          }
        }
        queue_head = queue_start;
        queue_tail = queue_start;
      }
      while (queue8_head!=queue8_tail);
    }
}

MIA_RESULT_CODE IPL_SEGM_FindMaximaDigging_uint8_ext(
  SValuePosition* maxima,       // OUT:    destination array
  int* size,           // IN/OUT: size of array in elements/number of extremums found
  uint8* array,         // IN:     image to be searched for maxima (zeroed on exit)
  int w,                        // IN:     image size
  int h,                        //
  uint8 maxthreshold, // IN: minimum possible brightness of local max to be detected
  uint8 mindrop_thr,  // IN: minimum possible brightness drop to be detected
  int maxobjsize,     // IN: maximum allowed object size
  void* extbuftmp)              // IN:     temporary buffer
{
  int queue_size,queue8_size;
  ptr_t *queue,*queue8;
  int  y,maxima_counter,edge,M,MX,MY;
  ptr_t *queue8_start,*queue8_finish,*queue8_tail,*queue8_head;
  ptr_t *queue_start,*queue_finish,*queue_tail,*queue_head;
  uint8 array_val,value,max_value;
  uint8 *ptr,*maxptr,*finish,*addr;
  uint8 minbrt,objsize;

    queue_size = queue8_size = w*h/(2*sizeof(ptr_t));
    queue = (ptr_t*)extbuftmp;
    queue8 = queue+queue_size;
    queue8_start  = queue8;
    queue8_finish = queue8_start+queue8_size;
    queue8_head   = queue8_start;
    queue8_tail   = queue8_start;
    queue_start  = queue;
    queue_finish = queue_start+queue_size;
    queue_head   = queue_start;
    queue_tail   = queue_start;

    // clear borders
    memset(array,0,w*sizeof(array[0]));
    memset(array+w*(h-1),0,w*sizeof(array[0]));
    for (y=(h-1)*w;y>0;y-=w)
      array[y] = 0;
    for (y=h*w-1;y>w;y-=w)
      array[y] = 0;
    // init
    maxima_counter = 0;
    ptr = array+w+1;
    finish = array+w*(h-1)-1;
    // run along array
    for (;;)
    {
      while ((ptr<finish)&&(*ptr==0))
        ptr++;
      if (ptr==finish)
      { // end up here
        // output element count
        *size = maxima_counter;
        // sort by value descending
        return ERR_OK;
      }
      // we've found a non-zero point!
      // is it a big enough?
      if (*ptr<maxthreshold)
      { // no
        *ptr = 0;
        continue;
      }
      // we've found good point!!!
      maxptr = ptr;
      // gradient climbing for max (8-connectivity)
      for (;;)
      {
        MASK_CLIMB_STEP_MAX(1);
        MASK_CLIMB_STEP_MAX(w);
        MASK_CLIMB_STEP_MAX(w+1);
        MASK_CLIMB_STEP_MAX(w-1);
        MASK_CLIMB_STEP_MAX(-1);
        MASK_CLIMB_STEP_MAX(-w);
        MASK_CLIMB_STEP_MAX(-w+1);
        MASK_CLIMB_STEP_MAX(-w-1);
        break;
      }
      // we have found a local maxima point
      // however it can contain several points - enumerate them.
      queue8_head = queue8_start;
      queue8_tail = queue8_start;
      M = MX = MY = 0;
      max_value = *maxptr;
      *queue8_tail++ = (ptr_t)(maxptr);
      *maxptr = 0;
      do
      {
        maxptr = (uint8*) *queue8_head++;
        M++;
        MX += (int)((maxptr-array)%w);
        MY += (int)((maxptr-array)/w);
        if (queue8_head==queue8_finish)
          return ENHANCED_ERROR_CODE(0x202);
        ADD_TO_QUEUE_EDGE_8_MAX(1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(w,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(w+1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(w-1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-w,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-w+1,uint8);
        ADD_TO_QUEUE_EDGE_8_MAX(-w-1,uint8);
      }
      while (queue8_tail!=queue8_head);
      // all the points of maxima are now in queue elements [0..qtail=qhead)
      // output to maxima array
      maxima[maxima_counter].value = max_value;
      maxima[maxima_counter].x = (MX+M/2)/M;
      maxima[maxima_counter].y = (MY+M/2)/M;
      maxima_counter++;
      if (maxima_counter==*size)
        return ERR_GEN_INSUFFICIENT_BUFFER;
      // set qhead=0 and start digging
      queue8_head = queue8_start;
      queue_head  = queue_start;
      queue_tail  = queue_start;
      minbrt = max_value;
      objsize = 0;
      do
      {
        value = max_value;
        maxptr = (uint8*) *queue8_head++;
        ADD_TO_QUEUE_EDGE_0_MAX;
        ADD_TO_QUEUE_EDGE_1_MAX;
        ADD_TO_QUEUE_EDGE_2_MAX;
        ADD_TO_QUEUE_EDGE_3_MAX;
        ADD_TO_QUEUE_EDGE_4_MAX;
        ADD_TO_QUEUE_EDGE_5_MAX;
        ADD_TO_QUEUE_EDGE_6_MAX;
        ADD_TO_QUEUE_EDGE_7_MAX;
        while (queue_head!=queue_tail)
        {
          edge  =  (uint8)(*queue_head++);
          value =  (uint8)(*queue_head++);
          maxptr = (uint8*) *queue_head++;
          if (queue_head>=queue_finish)
            queue_head=queue_start;
          if (minbrt>value)
            minbrt = value;
          objsize++;
          switch (edge)
          {
            case 0:
Case_Entry_0: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr=maxptr-1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_0;
                  }
              }
              break;
            case 1:
Case_Entry_1: ADD_TO_QUEUE_EDGE_0_MAX;
              ADD_TO_QUEUE_EDGE_2_MAX;
              ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr-w-1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_1;
                  }
              }
              break;
            case 2:
Case_Entry_2: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_3_MAX;
              addr=maxptr-w;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_2;
                  }
              }
              break;
            case 3:
Case_Entry_3: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_2_MAX;
              ADD_TO_QUEUE_EDGE_4_MAX;
              ADD_TO_QUEUE_EDGE_5_MAX;
              addr=maxptr-w+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_3;
                  }
              }
              break;
            case 4:
Case_Entry_4: ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_5_MAX;
              addr = maxptr+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_4;
                  }
              }
              break;
            case 5:
Case_Entry_5: ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_4_MAX;
              ADD_TO_QUEUE_EDGE_6_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr+w+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_5;
                  }
              }
              break;
            case 6:
Case_Entry_6: ADD_TO_QUEUE_EDGE_5_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr+w;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_6;
                  }
              }
              break;
            case 7:
Case_Entry_7: ADD_TO_QUEUE_EDGE_5_MAX;
              ADD_TO_QUEUE_EDGE_6_MAX;
              ADD_TO_QUEUE_EDGE_0_MAX;
              ADD_TO_QUEUE_EDGE_1_MAX;
              addr = maxptr+w-1;
              array_val = *addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_7;
                  }
              }
              break;
          }
        }
        queue_head = queue_start;
        queue_tail = queue_start;
      }
      while (queue8_head!=queue8_tail);
      if ((maxima[maxima_counter-1].value-minbrt<mindrop_thr)||
          (objsize>maxobjsize))
        maxima_counter--;
    }
}

MIA_RESULT_CODE IPL_SEGM_FindMaximaDigging_uint32(
  SValuePosition* maxima,       // OUT:    destination array
  int* size,           // IN/OUT: size of array in elements/number of extremums found
  uint32* array,         // IN:     image to be searched for maxima (zeroed on exit)
  int w,                        // IN:     image size
  int h,                        //
  uint32 maxthreshold,   // IN:     threshold of max
  void* extbuftmp)              // IN:     temporary buffer
{
  int queue_size,queue8_size;
  ptr_t *queue,*queue8;
  int  y,maxima_counter,edge,M,MX,MY;
  ptr_t *queue8_start,*queue8_finish,*queue8_tail,*queue8_head;
  ptr_t *queue_start,*queue_finish,*queue_tail,*queue_head;
  uint32 array_val,value,max_value;
  uint32 *ptr,*maxptr,*finish,*addr;

    queue_size = queue8_size = w*h/(2*sizeof(ptr_t));
    queue = (ptr_t*)extbuftmp;
    queue8 = queue+queue_size;
    queue8_start  = queue8;
    queue8_finish = queue8_start+queue8_size;
    queue8_head   = queue8_start;
    queue8_tail   = queue8_start;
    queue_start  = queue;
    queue_finish = queue_start+queue_size;
    queue_head   = queue_start;
    queue_tail   = queue_start;

    // clear borders
    memset(array,0,w*sizeof(array[0]));
    memset(array+w*(h-1),0,w*sizeof(array[0]));
    for (y=(h-1)*w;y>0;y-=w)
      array[y] = 0;
    for (y=h*w-1;y>w;y-=w)
      array[y] = 0;
    // init
    maxima_counter = 0;
    ptr = array+w+1;
    finish = array+w*(h-1)-1;
    // run along array
    for (;;)
    {
      while ((ptr<finish)&&(*ptr==0))
        ptr++;
      if (ptr==finish)
      { // end up here
        // output element count
        *size = maxima_counter;
        // sort by value descending
        return ERR_OK;
      }
      // we've found a non-zero point!
      // is it a big enough?
      if (*ptr<maxthreshold)
      { // no
        *ptr = 0;
        continue;
      }
      // we've found good point!!!
      maxptr = ptr;
      // gradient climbing for max (8-connectivity)
      for (;;)
      {
        MASK_CLIMB_STEP_MAX(1);
        MASK_CLIMB_STEP_MAX(w);
        MASK_CLIMB_STEP_MAX(w+1);
        MASK_CLIMB_STEP_MAX(w-1);
        MASK_CLIMB_STEP_MAX(-1);
        MASK_CLIMB_STEP_MAX(-w);
        MASK_CLIMB_STEP_MAX(-w+1);
        MASK_CLIMB_STEP_MAX(-w-1);
        break;
      }
      // we have found a local maxima point
      // however it can contain several points - enumerate them.
      queue8_head = queue8_start;
      queue8_tail = queue8_start;
      M = MX = MY = 0;
      max_value = *maxptr;
      *queue8_tail++ = (ptr_t)(maxptr);
      *maxptr = 0;
      do
      {
        maxptr = (uint32*) *queue8_head++;
        M++;
        MX += (int)((maxptr-array)%w);
        MY += (int)((maxptr-array)/w);
        if (queue8_head==queue8_finish)
          return ENHANCED_ERROR_CODE(0x202);
        ADD_TO_QUEUE_EDGE_8_MAX(1,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(w,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(w+1,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(w-1,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(-1,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(-w,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(-w+1,uint32);
        ADD_TO_QUEUE_EDGE_8_MAX(-w-1,uint32);
      }
      while (queue8_tail!=queue8_head);
      // all the points of maxima are now in queue elements [0..qtail=qhead)
      // output to maxima array
      maxima[maxima_counter].value = max_value;
      maxima[maxima_counter].x = (MX+M/2)/M;
      maxima[maxima_counter].y = (MY+M/2)/M;
      maxima_counter++;
      if (maxima_counter==*size)
        return ERR_GEN_INSUFFICIENT_BUFFER;
      // set qhead=0 and start digging
      queue8_head=queue8_start;
      queue_head  =queue_start;
      queue_tail  =queue_start;
      do
      {
        value = max_value;
        maxptr = (uint32*) *queue8_head++;
        ADD_TO_QUEUE_EDGE_0_MAX;
        ADD_TO_QUEUE_EDGE_1_MAX;
        ADD_TO_QUEUE_EDGE_2_MAX;
        ADD_TO_QUEUE_EDGE_3_MAX;
        ADD_TO_QUEUE_EDGE_4_MAX;
        ADD_TO_QUEUE_EDGE_5_MAX;
        ADD_TO_QUEUE_EDGE_6_MAX;
        ADD_TO_QUEUE_EDGE_7_MAX;
        while (queue_head!=queue_tail)
        {
          edge  =  (uint32)(*queue_head++);
          value =  (uint32)(*queue_head++);
          maxptr = (uint32*) *queue_head++;
          if (queue_head>=queue_finish)
            queue_head=queue_start;
          switch (edge)
          {
            case 0:
Case_Entry_0: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr=maxptr-1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_0;
                  }
              }
              break;
            case 1:
Case_Entry_1: ADD_TO_QUEUE_EDGE_0_MAX;
              ADD_TO_QUEUE_EDGE_2_MAX;
              ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr-w-1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_1;
                  }
              }
              break;
            case 2:
Case_Entry_2: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_3_MAX;
              addr=maxptr-w;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_2;
                  }
              }
              break;
            case 3:
Case_Entry_3: ADD_TO_QUEUE_EDGE_1_MAX;
              ADD_TO_QUEUE_EDGE_2_MAX;
              ADD_TO_QUEUE_EDGE_4_MAX;
              ADD_TO_QUEUE_EDGE_5_MAX;
              addr=maxptr-w+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr=addr;
                    *(addr) = 0;
                    goto Case_Entry_3;
                  }
              }
              break;
            case 4:
Case_Entry_4: ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_5_MAX;
              addr = maxptr+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value =array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_4;
                  }
              }
              break;
            case 5:
Case_Entry_5: ADD_TO_QUEUE_EDGE_3_MAX;
              ADD_TO_QUEUE_EDGE_4_MAX;
              ADD_TO_QUEUE_EDGE_6_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr+w+1;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_5;
                  }
              }
              break;
            case 6:
Case_Entry_6: ADD_TO_QUEUE_EDGE_5_MAX;
              ADD_TO_QUEUE_EDGE_7_MAX;
              addr = maxptr+w;
              array_val=*addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_6;
                  }
              }
              break;
            case 7:
Case_Entry_7: ADD_TO_QUEUE_EDGE_5_MAX;
              ADD_TO_QUEUE_EDGE_6_MAX;
              ADD_TO_QUEUE_EDGE_0_MAX;
              ADD_TO_QUEUE_EDGE_1_MAX;
              addr = maxptr+w-1;
              array_val = *addr;
              if (array_val)
              {
                if (array_val<maxthreshold)
                  *(addr) = 0;
                else
                  if (array_val<=value)
                  {
                    value = array_val;
                    maxptr = addr;
                    *(addr) = 0;
                    goto Case_Entry_7;
                  }
              }
              break;
          }
        }
        queue_head = queue_start;
        queue_tail = queue_start;
      }
      while (queue8_head!=queue8_tail);
    }
}
