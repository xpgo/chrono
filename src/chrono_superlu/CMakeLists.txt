#=============================================================================
# CMake configuration file for the Chrono SuperLU module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_SUPERLU "Enable the Chrono SuperLU module" OFF)

IF(NOT ENABLE_MODULE_SUPERLU)
    MARK_AS_ADVANCED(FORCE SUPERLU_ROOT)
	MARK_AS_ADVANCED(FORCE BLAS_ROOT)
    RETURN()
ENDIF()

MESSAGE(STATUS "==== Chrono SuperLU module ====")


MARK_AS_ADVANCED(CLEAR SUPERLU_ROOT)
MARK_AS_ADVANCED(CLEAR BLAS_ROOT)

# ------------------------------------------------------------------------------
# Threading
# ------------------------------------------------------------------------------

SET(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads QUIET)
find_package(OpenMP QUIET)

SET(CH_SUPERLU_DEFS)

if (CMAKE_USE_PTHREADS_INIT)
	MESSAGE(STATUS "Using pthread multithreading")
	SET(MULTITHREADING_LIBRARIES "Threads::Threads")
	LIST(APPEND CH_SUPERLU_DEFS "__PTHREAD")
elseif(OPENMP_FOUND)
	MESSAGE(STATUS "Using OpenMP multithreading")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	LIST(APPEND CH_SUPERLU_DEFS "__OPENMP")
else()
	MESSAGE(WARNING "No multithreading library has been found")
endif()


# ------------------------------------------------------------------------------
# Dependencies for SuperLU module
# ------------------------------------------------------------------------------
message(STATUS "Finding SuperLU libraries...")

# Find SuperLU
set(SUPERLU_ROOT "C:/Program Files/SuperLU_MT" CACHE PATH "Location of SuperLU installation")
find_package(SUPERLU)
message(STATUS "SUPERLU_ROOT is ${SUPERLU_ROOT}")
if(SUPERLU_FOUND)
	MARK_AS_ADVANCED(FORCE SUPERLU_ROOT)
endif()
message(STATUS "   SUPERLU include dirs:  ${SUPERLU_INCLUDE_DIRS}")
message(STATUS "   SUPERLU libraries:     ${SUPERLU_LIBRARIES}")



# Find BLAS (eventually looking into SuperLU folder)
set(BLAS_ROOT "${SUPERLU_ROOT}" CACHE PATH "Location of BLAS installation")
find_package(BLAS QUIET)
if(BLAS_FOUND)
	MARK_AS_ADVANCED(FORCE BLAS_ROOT)
else()
# custom, easy, straightforward search of blas library
	find_library(BLAS_LIBRARIES "blas"
		PATHS "${BLAS_ROOT}"
		PATH_SUFFIXES "lib"
		)
		
	find_path(BLAS_INCLUDE_DIRS "f2c.h"
		PATHS "${BLAS_ROOT}"
		PATH_SUFFIXES "include"
		)
		
  if (BLAS_LIBRARIES AND BLAS_INCLUDE_DIR)
	SET(BLAS_FOUND TRUE)
  endif (BLAS_LIBRARIES AND BLAS_INCLUDE_DIR)
endif()
message(STATUS "   BLAS include dirs:     ${BLAS_INCLUDE_DIRS}")
message(STATUS "   BLAS library:          ${BLAS_LIBRARIES}")



# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the SUPERLU module
# ------------------------------------------------------------------------------

SET(CH_SUPERLU_INCLUDES ${SUPERLU_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CH_SUPERLU_INCLUDES})
SET(CH_SUPERLU_INCLUDES "${CH_SUPERLU_INCLUDES}" PARENT_SCOPE)
SET(CH_SUPERLU_DEFS "${CH_SUPERLU_DEFS}" PARENT_SCOPE)

# ------------------------------------------------------------------------------
# List all files in the Chrono SUPERLU module
# ------------------------------------------------------------------------------

SET(ChronoEngine_SuperLU_HEADERS
  ChApiSuperLU.h
  ChSolverSuperLU.h
  ChSuperLUEngine.h
)

SET(ChronoEngine_SuperLU_SOURCES
	ChSuperLUEngine.cpp
)

SOURCE_GROUP("Headers" FILES ${ChronoEngine_SuperLU_HEADERS})
SOURCE_GROUP("Sources" FILES ${ChronoEngine_SuperLU_SOURCES})


# ------------------------------------------------------------------------------
# Add the ChronoEngine_SuperLU library
# ------------------------------------------------------------------------------

ADD_LIBRARY(ChronoEngine_SuperLU SHARED
            ${ChronoEngine_SuperLU_SOURCES}
            ${ChronoEngine_SuperLU_HEADERS})

MESSAGE(STATUS "CH_SUPERLU_DEFS is ${CH_SUPERLU_DEFS}")
SET_TARGET_PROPERTIES(ChronoEngine_SuperLU PROPERTIES
                      COMPILE_FLAGS "${COMPILER_FLAGS}"
                      LINK_FLAGS "${LINKER_FLAGS}"
                      COMPILE_DEFINITIONS "CH_API_COMPILE_SUPERLU;${CH_SUPERLU_DEFS}")

TARGET_LINK_LIBRARIES(ChronoEngine_SuperLU
                      ChronoEngine
                      ${SUPERLU_LIBRARIES}
					  ${BLAS_LIBRARIES}
					  ${MULTITHREADING_LIBRARIES}
					  )

INSTALL(TARGETS ChronoEngine_SuperLU
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

INSTALL(FILES ${ChronoEngine_SuperLU_HEADERS} 
        DESTINATION include/chrono_superlu)
