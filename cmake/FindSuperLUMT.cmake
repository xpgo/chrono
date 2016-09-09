# .. cmake_module::
#
#    Module that checks whether SuperLU_MT is available and usable.
#    SuperLU_MT must be a version released after the year 2005.
#
#    Variables used by this module which you may want to set:
#
#    :ref:`SUPERLUMT_ROOT`
#       Path list to search for SuperLU_MT
#
#    Sets the follwing variables:
#
#    :code:`SUPERLUMT_FOUND`
#       True if SuperLU_MT available and usable.
#
#
#    :code:`SUPERLUMT_INCLUDE_DIRS`
#       Path to the SuperLU_MT include dirs.
#
#    :code:`SUPERLUMT_LIBRARIES`
#       Name to the SuperLU_MT library.
#
# .. cmake_variable:: SUPERLUMT_ROOT
#
#    You may set this variable to have :ref:`FindSuperLU_MT` look
#    for the SuperLU_MT package in the given path before inspecting
#    system paths.
#

# look for header files, only at positions given by the user
find_path(SUPERLUMT_INCLUDE_DIR
  NAMES supermatrix.h
  PATHS ${SUPERLUMT_PREFIX} ${SUPERLUMT_ROOT}
  PATH_SUFFIXES "superlu_mt" "SuperLU_MT" "include/superlu_mt" "include" "SRC"
  NO_DEFAULT_PATH
)

# look for header files, including default paths
find_path(SUPERLUMT_INCLUDE_DIR
  NAMES supermatrix.h
  PATH_SUFFIXES "superlu_mt" "SuperLU_MT" "include/superlu_mt" "include" "SRC"
)

# look for library, only at positions given by the user
find_library(SUPERLUMT_LIBRARY
  NAMES "superlu_mt_3.0" "superlu_mt"
  PATHS ${SUPERLUMT_PREFIX} ${SUPERLUMT_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  NO_DEFAULT_PATH
)

# look for library files, including default paths
find_library(SUPERLUMT_LIBRARY
  NAMES "superlu_mt_3.0" "superlu_mt"
  PATH_SUFFIXES "lib" "lib32" "lib64"
)


# check version specific macros
include(CheckCSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND
# if the searches above were not successful
# Without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they are set to NOTFOUND.
# Please set them or make sure they are set and tested correctly in the CMake files:"
#
if(SUPERLUMT_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${SUPERLUMT_INCLUDE_DIR})
endif(SUPERLUMT_INCLUDE_DIR)
if(SUPERLUMT_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${SUPERLUMT_LIBRARY})
endif(SUPERLUMT_LIBRARY)
if(BLAS_LIBRARIES)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${BLAS_LIBRARIES})
endif(BLAS_LIBRARIES)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "SuperLUMT"
  DEFAULT_MSG
  SUPERLUMT_INCLUDE_DIR
  SUPERLUMT_LIBRARY
)

mark_as_advanced(SUPERLUMT_INCLUDE_DIR SUPERLUMT_LIBRARY)

# if both headers and library are found, store results
if(SUPERLUMT_FOUND)
  if (NOT SUPERLUMT_ROOT)
	SET (SUPERLUMT_ROOT "${SUPERLUMT_LIBRARIES}/../")
  endif()
  set(SUPERLUMT_INCLUDE_DIRS ${SUPERLUMT_INCLUDE_DIR})
  set(SUPERLUMT_LIBRARIES    ${SUPERLUMT_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of ${SUPERLUMT_WITH_VERSION} succeeded:\n"
    "Include directory: ${SUPERLUMT_INCLUDE_DIRS}\n"
    "Library directory: ${SUPERLUMT_LIBRARIES}\n\n")
else(SUPERLUMT_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of SuperLU_MT failed:\n"
    "Include directory: ${SUPERLUMT_INCLUDE_DIRS}\n"
    "Library directory: ${SUPERLUMT_LIBRARIES}\n\n")
endif(SUPERLUMT_FOUND)

# set HAVE_SUPERLUMT for config.h
set(HAVE_SUPERLUMT ${SUPERLUMT_FOUND})