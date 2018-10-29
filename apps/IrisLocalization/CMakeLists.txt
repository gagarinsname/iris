#
# Iris Recognition Library
#

# apps/IrisLocalization

set(LIB_NAME IrisLocalization)

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC
	${SRC_DIR}/CoreMatveev/canny.c
    ${SRC_DIR}/CoreMatveev/conobj.c
	${SRC_DIR}/BaseRadii.cpp
    ${SRC_DIR}/imageio.c
    ${SRC_DIR}/ImageProcessing.c
    ${SRC_DIR}/ImageProcessing.h
    ${SRC_DIR}/imPlot.cpp
    ${SRC_DIR}/PairedGradients.c
    ${SRC_DIR}/ProjectionProcessing.cpp
)
set(HEADERS
	${SRC_DIR}/CoreMatveev/bdl.h
	${SRC_DIR}/CoreMatveev/bpl.h
    ${SRC_DIR}/CoreMatveev/dbg.h
    ${SRC_DIR}/CoreMatveev/errorcodes.h
	${SRC_DIR}/CoreMatveev/ipl.h
    ${SRC_DIR}/CoreMatveev/ownipl.h
	${SRC_DIR}/CoreMatveev/ffl.h
    ${SRC_DIR}/BaseRadii.h
    ${SRC_DIR}/BorderDetection.h
    ${SRC_DIR}/Geometry.h
    ${SRC_DIR}/GlobalParams.h
    ${SRC_DIR}/ImageProcessing.h
    ${SRC_DIR}/imio.h
    ${SRC_DIR}/imPlot.h
    ${SRC_DIR}/PairedGradients.h
    ${SRC_DIR}/stddefs.h
    ${SRC_DIR}/tables.h
    ${SRC_DIR}/unitypes.h
)

include_directories (
    ${SRC_DIR}
)

add_library(${LIB_NAME} STATIC ${SRC} ${HEADERS})

target_include_directories(${LIB_NAME} INTERFACE ${SRC_DIR})

SET_PROPERTY(TARGET ${LIB_NAME} PROPERTY FOLDER apps)