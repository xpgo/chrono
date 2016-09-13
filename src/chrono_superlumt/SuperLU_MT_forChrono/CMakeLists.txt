######################################################################
#
# CMakeLists.txt for SUPERLU_MT
#
######################################################################

# Required version
cmake_minimum_required(VERSION 3.4)

# Project Version	
project(SuperLU_MT NONE)
set(VERSION_MAJOR "3")
set(VERSION_MINOR "0")
set(VERSION_BugFix "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# Usual initialization stuff
#
######################################################################
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????
set(CMAKE_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_PREFIX}/include")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

enable_language(C)
enable_language(CXX)
enable_language (Fortran)
set(NOFORTRAN FALSE)
set(SUPERLU_VERSION "${PROJECT_VERSION}")
set(SUPERLU_REV "${PROJECT_REV}")

# setup options
option(build_internal_blas   "Build the CBLAS library" ON)
option(BUILD_SHARED_LIBS   "Build shared libraries" OFF)
option(enable_matlabmex "Build the Matlab mex library" OFF)
option(enable_tests     "Build tests" OFF)
option(enable_doc       "Build doxygen documentation" OFF)
option(enable_single    "Enable single precision library" OFF)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" OFF)
option(enable_complex16 "Enable complex16 precision library" OFF)
option(enable_examples  "Build examples" OFF)

# multithreading settings
SET(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads QUIET)
find_package(OpenMP QUIET)

SET(SUPERLU_DEFS "Add_;PRNTlevel=0")
if (CMAKE_USE_PTHREADS_INIT)
	MESSAGE(STATUS "Using pthread multithreading")
	SET(MULTITHREADING_LIBRARIES "Threads::Threads")
	LIST(APPEND SUPERLU_DEFS "__PTHREAD")
elseif(OPENMP_FOUND)
	MESSAGE(STATUS "Using OpenMP multithreading")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	LIST(APPEND SUPERLU_DEFS "__OPENMP")
else()
	MESSAGE(WARNING "No multithreading library has been found")
endif()

# shared libraries tweak
if (BUILD_SHARED_LIBS)
	SET(WINDOWS_EXPORT_ALL_SYMBOLS)
	if (MSVC)
		MESSAGE(WARNING "Shared libraries dont't work yet on MSVC. Switching to static.")
		SET(BUILD_SHARED_LIBS OFF)
	endif(MSVC)
endif (BUILD_SHARED_LIBS)


######################################################################
#
# Find packages
#
######################################################################
#
#--------------------- BLAS ---------------------
if(build_internal_blas)
  add_subdirectory(CBLAS)
  set(BLAS_LIBRARIES blas)
  if (BUILD_SHARED_LIBS)  # export to be referenced by downstream makefile
      set(BLAS_LIB_EXPORT ${CMAKE_SOURCE_DIR}/build/CBLAS/libblas.so)
  else(BUILD_SHARED_LIBS)
      set(BLAS_LIB_EXPORT ${CMAKE_SOURCE_DIR}/build/CBLAS/libblas.a)
  endif(BUILD_SHARED_LIBS)
else(build_internal_blas)
	MESSAGE(STATUS "Looking for BLAS with Find")
    find_package(BLAS)
    if (BLAS_FOUND)
	  set(CMAKE_C_FLAGS "-DUSE_VENDOR_BLAS ${CMAKE_C_FLAGS}")
	  set(BLAS_INCLUDE_DIRS "${BLAS_INCLUDE_DIR}")
	else(BLAS_FOUND)
		message(FATAL_ERROR "-- Did not find any specify BLAS and not able to build internal BLAS ...")
    endif(BLAS_FOUND)
endif()

MESSAGE(STATUS "BLAS_INCLUDE_DIRS is ${BLAS_INCLUDE_DIRS}")
MESSAGE(STATUS "BLAS_LIBRARIES is ${BLAS_LIBRARIES}")


######################################################################
#
# Include directories
#
######################################################################
SET(SUPERLU_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/SRC")


######################################################################
#
# Add subdirectories
#
######################################################################

include_directories(${BLAS_INCLUDE_DIRS})

add_subdirectory(SRC)

include_directories(${SUPERLU_INCLUDE_DIRS})

if (MSVC)
	LIST(APPEND SUPERLU_DEFS "_CRT_SECURE_NO_WARNINGS")
endif()

if(enable_tests)
  message(FATAL_ERROR "Testing build requested but not implemented.")
  # enable_testing()
  # add_subdirectory(TESTING)
endif()

if(enable_doc)
  message(FATAL_ERROR "Documentation build requested but not implemented.")
  #implement doxygen
endif()

if(enable_examples)
  add_subdirectory(EXAMPLE)
endif()
 
