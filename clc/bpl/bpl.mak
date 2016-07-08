#============================= library file lists ===============================
DIR_BPL = ./src

CFILES_BPL = \
  $(DIR_BPL)/errate/errate.c \
  $(DIR_BPL)/errate/statbld.c \
  $(DIR_BPL)/linalg/eigen.c \
  $(DIR_BPL)/linalg/ellproc.c \
  $(DIR_BPL)/linalg/invert.c \
  $(DIR_BPL)/linalg/linsys.c \
  $(DIR_BPL)/linalg/matrbase.c \
  $(DIR_BPL)/linalg/matrchk.c \
  $(DIR_BPL)/linalg/matrsolv.c \
  $(DIR_BPL)/misc/crc32.c \
  $(DIR_BPL)/misc/LevMar.c \
  $(DIR_BPL)/misc/mathmix.c \
  $(DIR_BPL)/misc/miscfn.c \
  $(DIR_BPL)/pca/SLDBuilder.c \
  $(DIR_BPL)/pca/SNDBuilder.c \
  $(DIR_BPL)/pca/SNDClassifier.c \
  $(DIR_BPL)/pca/SPCAnalyser.c \
  $(DIR_BPL)/pca/SPCBuilder.c \
  $(DIR_BPL)/pca/svd.c

CFILES_LIB = $(CFILES_BPL)

DEPFILES_BPL = \
  $(DIR_BPL)/errate/errate.d \
  $(DIR_BPL)/errate/statbld.d \
  $(DIR_BPL)/linalg/eigen.d \
  $(DIR_BPL)/linalg/ellproc.d \
  $(DIR_BPL)/linalg/invert.d \
  $(DIR_BPL)/linalg/linsys.d \
  $(DIR_BPL)/linalg/matrbase.d \
  $(DIR_BPL)/linalg/matrchk.d \
  $(DIR_BPL)/linalg/matrsolv.d \
  $(DIR_BPL)/misc/crc32.d \
  $(DIR_BPL)/misc/LevMar.d \
  $(DIR_BPL)/misc/mathmix.d \
  $(DIR_BPL)/misc/miscfn.d \
  $(DIR_BPL)/pca/SLDBuilder.d \
  $(DIR_BPL)/pca/SNDBuilder.d \
  $(DIR_BPL)/pca/SNDClassifier.d \
  $(DIR_BPL)/pca/SPCAnalyser.d \
  $(DIR_BPL)/pca/SPCBuilder.d \
  $(DIR_BPL)/pca/svd.d

DEPFILES_LIB = $(DEPFILES_BPL)

OBJFILES_BPL = \
  $(DIR_BPL)/errate/errate.$(OBJ_EXT) \
  $(DIR_BPL)/errate/statbld.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/eigen.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/ellproc.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/invert.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/linsys.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/matrbase.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/matrchk.$(OBJ_EXT) \
  $(DIR_BPL)/linalg/matrsolv.$(OBJ_EXT) \
  $(DIR_BPL)/misc/crc32.$(OBJ_EXT) \
  $(DIR_BPL)/misc/LevMar.$(OBJ_EXT) \
  $(DIR_BPL)/misc/mathmix.$(OBJ_EXT) \
  $(DIR_BPL)/misc/miscfn.$(OBJ_EXT) \
  $(DIR_BPL)/pca/SLDBuilder.$(OBJ_EXT) \
  $(DIR_BPL)/pca/SNDBuilder.$(OBJ_EXT) \
  $(DIR_BPL)/pca/SNDClassifier.$(OBJ_EXT) \
  $(DIR_BPL)/pca/SPCAnalyser.$(OBJ_EXT) \
  $(DIR_BPL)/pca/SPCBuilder.$(OBJ_EXT) \
  $(DIR_BPL)/pca/svd.$(OBJ_EXT)

OBJFILES_LIB = $(OBJFILES_BPL)

#============================= library files ===============================
DIR_OUT = ../lib

LIBFILE_NAME = libBPL_$(SUFFIX).a

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
	$(CC) $(CCOPTS) -c -O3 -MMD -MP -fPIC -I ../inc -o $@ $<
	@echo ' '

all: printsign $(DIR_OUT) $(LIBFILE_PATH)
	@echo ' '

$(LIBFILE_PATH): $(OBJFILES_LIB)
	ar -r $@ $(OBJFILES_LIB)

printsign:
	@echo ----------------------------------------------------------------
	@echo ------------------ Building BPL library ------------------------
	@echo ----------------------------------------------------------------
	
$(DIR_OUT):
	mkdir $(DIR_OUT)
	@echo ' '

clean:
	$(RM) $(GENERATEDFILES_ALL)
	@echo ' '
