#
# Iris Recognition Library
#
# Authors: Iurii Efimov <yuri.efimov@phystech.edu>
# 

# root

cmake_minimum_required (VERSION 2.8.11)

set(PROJECT_NAME IrisRecognition)

project(${PROJECT_NAME})

set(bitness 32)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(bitness 64)
endif()

message(STATUS "Bitness: ${bitness}")


FIND_PACKAGE(PythonInterp)

if (PYTHONINTERP_FOUND)
  if (UNIX AND NOT APPLE)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  else()	
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  endif()
else()
    message("Python not found")
endif()

option(IRL_PYTHON_WRAPPERS "Build Cython wrapper for IrisRecognitionLibrary" ON)

# --- Python Support ---
# include(cmake/DetectPython.cmake)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMake")

if (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(SEND_ERROR "In-source builds are not allowed.")
endif()

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

if (WIN32)
    add_definitions(/MP) # multiprocessor compilation

    string(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

    string(REPLACE "/MDd" "/MTd" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
    string(REPLACE "/MD" "/MT" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
 
    # store debug information on .obj instead of .pdb
    string(REGEX REPLACE "/Z." "/Z7" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    string(REGEX REPLACE "/Z." "/Z7" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

elseif(ANDROID)
    set(common_cflags "-std=c99 -fsigned-char")
    set(common_cxxflags "-std=c++11 -fsigned-char")

    if (bitness EQUAL 64)
        enable_language(CXX C ASM)

        set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${common_cflags} -fPIE -flax-vector-conversions -Ofast -mfpu=neon")
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${common_cxxflags} -fPIE -flax-vector-conversions -Ofast -mfpu=neon ")
        set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -fPIE -pie -Ofast -mfpu=neon")
    
        set (CMAKE_ASM_FLAGS_RELEASE "${CMAKE_ASM_FLAGS_RELEASE}  -x assembler-with-cpp -Ofast -mcpu=krait  -mfpu=neon")

        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${common_cflags} -fPIE -flax-vector-conversions -Ofast")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${common_cxxflags} -fPIE -flax-vector-conversions -Ofast")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIE -pie -Ofast")
    
        set (CMAKE_ASM_FLAGS "${CMAKE_ASM_FLAGS}  -x assembler-with-cpp -Ofast -mcpu=krait  -mfpu=neon")

        add_definitions(-DNEON=1)

    else()
        enable_language(CXX C ASM)
        set(cflags "-fPIE -flax-vector-conversions -Ofast  -mfpu=neon")

        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${common_cflags} ${cflags}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${common_cxxflags} ${cflags} -std=c++11")

        set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -fPIE -pie -Ofast")

        set (CMAKE_ASM_FLAGS_RELEASE "${CMAKE_ASM_FLAGS_RELEASE}  -x assembler-with-cpp -Ofast -mfpu=neon -mcpu=krait")

        add_definitions(-D_USE_ASM_=1)
        add_definitions(-DNEON=1)
    endif()

elseif(LINUX_ARM)
    enable_language(CXX C ASM)

    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fPIE -flax-vector-conversions -Ofast -mfpu=neon -mcpu=cortex-a15")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIE -flax-vector-conversions -Ofast -mfpu=neon -mfpu=neon -mcpu=cortex-a15")
    set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -fPIE -pie -Ofast -mfpu=neon -mcpu=cortex-a15")

    set (CMAKE_ASM_FLAGS_RELEASE "${CMAKE_ASM_FLAGS_RELEASE}  -x assembler-with-cpp -Ofast -mfpu=neon -mcpu=cortex-a15")

    set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mfpu=neon -flax-vector-conversions -mcpu=cortex-a15")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mfpu=neon -flax-vector-conversions -mcpu=cortex-a15")

    add_definitions(-D_USE_ASM_=1)
    add_definitions(-DNEON=1)
elseif(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -msse3 -std=c++11")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fPIC")
endif()

add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set(EXECUTABLE_OUTPUT_PATH   ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH      ${PROJECT_BINARY_DIR}/lib)
set(RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# set search paths

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

link_directories(${LIBRARY_OUTPUT_PATH})

include(cmake/add_optional_subdir.cmake)

add_subdirectory(apps)
# add_optional_subdirectory(test)

# write build command to file, useful for automatic building
build_command(bldcmd)
file(WRITE "${CMAKE_BINARY_DIR}/build.cmd" 
    "cd ${CMAKE_BINARY_DIR}\n" 
    ${bldcmd})
