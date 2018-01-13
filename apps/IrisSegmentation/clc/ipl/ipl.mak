#============================= library file lists ===============================
DIR_IPL = src

CFILES_IPL = \
  $(DIR_IPL)/Ajust/affine.c \
  $(DIR_IPL)/Ajust/lines.c \
  $(DIR_IPL)/Bit/BitMorph.c \
  $(DIR_IPL)/CVL/BlStrtch.c \
  $(DIR_IPL)/CVL/DrawCirc.c \
  $(DIR_IPL)/CVL/GradProc.c \
  $(DIR_IPL)/CVL/pol.c \
  $(DIR_IPL)/CVL/ROIProc.c \
  $(DIR_IPL)/CVL/Scale.c \
  $(DIR_IPL)/CVL/Sobel.c \
  $(DIR_IPL)/Geom/elldet.c \
  $(DIR_IPL)/Geom/geom.c \
  $(DIR_IPL)/Hist/hist.c \
  $(DIR_IPL)/Pumps/canny.c \
  $(DIR_IPL)/Pumps/diffus.c \
  $(DIR_IPL)/Pumps/edgedraw.c \
  $(DIR_IPL)/Pumps/filt.c \
  $(DIR_IPL)/Pumps/ordstat.c \
  $(DIR_IPL)/Pumps/pmp.c \
  $(DIR_IPL)/Segm/blob.c \
  $(DIR_IPL)/Segm/cascade.c \
  $(DIR_IPL)/Segm/chunx1.c \
  $(DIR_IPL)/Segm/chunx2.c \
  $(DIR_IPL)/Segm/chunx_lin_hc.c \
  $(DIR_IPL)/Segm/chunx_lin_ff.c \
  $(DIR_IPL)/Segm/cluster.c \
  $(DIR_IPL)/Segm/conobj.c \
  $(DIR_IPL)/Segm/digging.c \
  $(DIR_IPL)/Segm/flash.c \
  $(DIR_IPL)/Segm/proj.c \
  $(DIR_IPL)/Spectral/fft.c \
  $(DIR_IPL)/Spectral/fht.c \
  $(DIR_IPL)/Spectral/fwt.c \
  $(DIR_IPL)/Spectral/pyr.c \
  $(DIR_IPL)/iplmisc.c \

CFILES_LIB = $(CFILES_IPL)

DEPFILES_IPL = \
  $(DIR_IPL)/Ajust/affine.d \
  $(DIR_IPL)/Ajust/lines.d \
  $(DIR_IPL)/Bit/BitMorph.d \
  $(DIR_IPL)/CVL/BlStrtch.d \
  $(DIR_IPL)/CVL/DrawCirc.d \
  $(DIR_IPL)/CVL/GradProc.d \
  $(DIR_IPL)/CVL/pol.d \
  $(DIR_IPL)/CVL/ROIProc.d \
  $(DIR_IPL)/CVL/Scale.d \
  $(DIR_IPL)/CVL/Sobel.d \
  $(DIR_IPL)/Geom/elldet.d \
  $(DIR_IPL)/Geom/geom.d \
  $(DIR_IPL)/Hist/hist.d \
  $(DIR_IPL)/Pumps/canny.d \
  $(DIR_IPL)/Pumps/diffus.d \
  $(DIR_IPL)/Pumps/edgedraw.d \
  $(DIR_IPL)/Pumps/filt.d \
  $(DIR_IPL)/Pumps/ordstat.d \
  $(DIR_IPL)/Pumps/pmp.d \
  $(DIR_IPL)/Segm/blob.d \
  $(DIR_IPL)/Segm/cascade.d \
  $(DIR_IPL)/Segm/chunx1.d \
  $(DIR_IPL)/Segm/chunx2.d \
  $(DIR_IPL)/Segm/chunx_lin_hc.d \
  $(DIR_IPL)/Segm/chunx_lin_ff.d \
  $(DIR_IPL)/Segm/cluster.d \
  $(DIR_IPL)/Segm/conobj.d \
  $(DIR_IPL)/Segm/digging.d \
  $(DIR_IPL)/Segm/flash.d \
  $(DIR_IPL)/Segm/flash.d \
  $(DIR_IPL)/Spectral/fft.d \
  $(DIR_IPL)/Spectral/fht.d \
  $(DIR_IPL)/Spectral/fwt.d \
  $(DIR_IPL)/Spectral/pyr.d \
  $(DIR_IPL)/iplmisc.d \

DEPFILES_LIB = $(DEPFILES_IPL)

OBJFILES_IPL = \
  $(DIR_IPL)/Ajust/affine.$(OBJ_EXT) \
  $(DIR_IPL)/Ajust/lines.$(OBJ_EXT) \
  $(DIR_IPL)/Bit/BitMorph.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/BlStrtch.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/DrawCirc.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/GradProc.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/pol.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/ROIProc.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/Scale.$(OBJ_EXT) \
  $(DIR_IPL)/CVL/Sobel.$(OBJ_EXT) \
  $(DIR_IPL)/Geom/elldet.$(OBJ_EXT) \
  $(DIR_IPL)/Geom/geom.$(OBJ_EXT) \
  $(DIR_IPL)/Hist/hist.$(OBJ_EXT) \
  $(DIR_IPL)/Pumps/canny.$(OBJ_EXT) \
  $(DIR_IPL)/Pumps/diffus.$(OBJ_EXT) \
  $(DIR_IPL)/Pumps/edgedraw.$(OBJ_EXT) \
  $(DIR_IPL)/Pumps/filt.$(OBJ_EXT) \
  $(DIR_IPL)/Pumps/ordstat.$(OBJ_EXT) \
  $(DIR_IPL)/Pumps/pmp.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/blob.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/cascade.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/chunx1.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/chunx2.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/chunx_lin_hc.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/chunx_lin_ff.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/cluster.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/conobj.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/digging.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/flash.$(OBJ_EXT) \
  $(DIR_IPL)/Segm/proj.$(OBJ_EXT) \
  $(DIR_IPL)/Spectral/fft.$(OBJ_EXT) \
  $(DIR_IPL)/Spectral/fht.$(OBJ_EXT) \
  $(DIR_IPL)/Spectral/fwt.$(OBJ_EXT) \
  $(DIR_IPL)/Spectral/pyr.$(OBJ_EXT) \
  $(DIR_IPL)/iplmisc.$(OBJ_EXT) \

OBJFILES_LIB = $(OBJFILES_IPL)

#============================= library files ===============================
DIR_OUT = ../lib

LIBFILE_NAME = libIPL_$(SUFFIX).a

LIBFILE_PATH = $(DIR_OUT)/$(LIBFILE_NAME)

GENERATEDFILES_ALL = $(DEPFILES_LIB) $(OBJFILES_LIB) $(LIBFILE_PATH)

#====================== options ==============================
ifeq ($(bits),64)
 ifeq ($(compiler),intel)
  OBJ_EXT = o64i
  SUFFIX = icc_x64
  CC = $(INTELCC)/bin/intel64/icc 
  CCOPTS = -m64
 else
  OBJ_EXT = o64g
  SUFFIX = gcc_x64
  CC = gcc
  CCOPTS = -m64
 endif
else
 # default bit is 32
 ifeq ($(compiler),intel)
  OBJ_EXT = o32i
  SUFFIX = icc_x32
  CC = $(INTELCC)/bin/ia32/icc 
  CCOPTS = -m32
 else
  # default compiler is GNU
  OBJ_EXT = o32g
  SUFFIX = gcc_x32
  CC = gcc
  CCOPTS = -m32
 endif
endif

RM := rm -rf

#======================= dependencies ==========================
%.$(OBJ_EXT): %.c
	$(CC) $(CCOPTS) -c -O3 -MMD -MP -fPIC -I ../../alc/inc -I ../inc -I ../../api/inc -o $@ $<
	@echo ' '

all: printsign $(DIR_OUT) $(LIBFILE_PATH)
	@echo ' '

$(LIBFILE_PATH): $(OBJFILES_LIB)
	ar -r $@ $(OBJFILES_LIB)

printsign:
	@echo ----------------------------------------------------------------
	@echo ------------------ Building IPL library ------------------------
	@echo ----------------------------------------------------------------
	
$(DIR_OUT):
	mkdir $(DIR_OUT)
	@echo ' '

clean:
	$(RM) $(GENERATEDFILES_ALL)
	@echo ' '
