#=============================================================================
# CMake configuration file for the Chrono SuperLUMT module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_SUPERLUMT "Enable the Chrono SuperLUMT module" OFF)

IF(NOT ENABLE_MODULE_SUPERLUMT)
    MARK_AS_ADVANCED(FORCE SUPERLUMT_ROOT)
	MARK_AS_ADVANCED(FORCE BLAS_ROOT)
    RETURN()
ENDIF()

MESSAGE(STATUS "==== Chrono SuperLUMT module ====")


MARK_AS_ADVANCED(CLEAR SUPERLUMT_ROOT)
MARK_AS_ADVANCED(CLEAR BLAS_ROOT)

# ------------------------------------------------------------------------------
# Threading
# ------------------------------------------------------------------------------

SET(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads QUIET)
find_package(OpenMP QUIET)

SET(CH_SUPERLUMT_DEFS)
if (CMAKE_USE_PTHREADS_INIT)
	MESSAGE(STATUS "Using pthread multithreading")
	SET(MULTITHREADING_LIBRARIES "Threads::Threads")
	LIST(APPEND CH_SUPERLUMT_DEFS "__PTHREAD")
elseif(OPENMP_FOUND)
	MESSAGE(STATUS "Using OpenMP multithreading")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	LIST(APPEND CH_SUPERLUMT_DEFS "__OPENMP")
else()
	MESSAGE(WARNING "No multithreading library has been found")
endif()


# ------------------------------------------------------------------------------
# Dependencies for SuperLUMT module
# ------------------------------------------------------------------------------
message(STATUS "Finding SuperLUMT libraries...")

# Find SuperLUMT
set(SUPERLUMT_ROOT "C:/Program Files/SuperLU_MT" CACHE PATH "Location of SuperLUMT installation")
find_package(SUPERLUMT)
message(STATUS "SUPERLUMT_ROOT is ${SUPERLUMT_ROOT}")
if(SUPERLUMT_FOUND)
	MARK_AS_ADVANCED(FORCE SUPERLUMT_ROOT)
endif()
message(STATUS "   SUPERLUMT include dirs:  ${SUPERLUMT_INCLUDE_DIRS}")
message(STATUS "   SUPERLUMT libraries:     ${SUPERLUMT_LIBRARIES}")



# Find BLAS (eventually looking into SuperLUMT folder)
set(BLAS_ROOT "${SUPERLUMT_ROOT}" CACHE PATH "Location of BLAS installation")
find_package(BLAS QUIET)
if(BLAS_FOUND)
	MARK_AS_ADVANCED(FORCE BLAS_ROOT)
else() # custom, easy, straightforward search of blas library
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
	MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES)
	MARK_AS_ADVANCED(FORCE BLAS_INCLUDE_DIRS)
  endif (BLAS_LIBRARIES AND BLAS_INCLUDE_DIR)
endif()
message(STATUS "   BLAS include dirs:     ${BLAS_INCLUDE_DIRS}")
message(STATUS "   BLAS library:          ${BLAS_LIBRARIES}")



# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the SuperLUMT module
# ------------------------------------------------------------------------------

SET(CH_SUPERLUMT_INCLUDES ${SUPERLUMT_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CH_SUPERLUMT_INCLUDES})
SET(CH_SUPERLUMT_INCLUDES "${CH_SUPERLUMT_INCLUDES}" PARENT_SCOPE)
SET(CH_SUPERLUMT_DEFS "${CH_SUPERLUMT_DEFS}" PARENT_SCOPE)

# ------------------------------------------------------------------------------
# List all files in the Chrono SuperLUMT module
# ------------------------------------------------------------------------------

SET(ChronoEngine_SuperLUMT_HEADERS
  ChApiSuperLUMT.h
  ChSolverSuperLUMT.h
  ChSuperLUMTEngine.h
)

SET(ChronoEngine_SuperLUMT_SOURCES
	ChSuperLUMTEngine.cpp
)

SOURCE_GROUP("Headers" FILES ${ChronoEngine_SuperLUMT_HEADERS})
SOURCE_GROUP("Sources" FILES ${ChronoEngine_SuperLUMT_SOURCES})


# ------------------------------------------------------------------------------
# Add the ChronoEngine_SuperLUMT library
# ------------------------------------------------------------------------------

ADD_LIBRARY(ChronoEngine_SuperLUMT SHARED
            ${ChronoEngine_SuperLUMT_SOURCES}
            ${ChronoEngine_SuperLUMT_HEADERS})

MESSAGE(STATUS "CH_SUPERLUMT_DEFS is ${CH_SUPERLUMT_DEFS}")
SET_TARGET_PROPERTIES(ChronoEngine_SuperLUMT PROPERTIES
                      COMPILE_FLAGS "${COMPILER_FLAGS}"
                      LINK_FLAGS "${LINKER_FLAGS}"
                      COMPILE_DEFINITIONS "CH_API_COMPILE_SUPERLUMT;${CH_SUPERLUMT_DEFS}")

TARGET_LINK_LIBRARIES(ChronoEngine_SuperLUMT
                      ChronoEngine
                      ${SUPERLUMT_LIBRARIES}
					  ${BLAS_LIBRARIES}
					  ${MULTITHREADING_LIBRARIES}
					  )

INSTALL(TARGETS ChronoEngine_SuperLUMT
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

INSTALL(FILES ${ChronoEngine_SuperLUMT_HEADERS} 
        DESTINATION include/chrono_superlumt)
