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


SET(CH_SUPERLUMT_DEFS "PRNTlevel=0 Add_")
if (CMAKE_USE_PTHREADS_INIT)
	LIST(APPEND CH_SUPERLUMT_DEFS "__PTHREAD")
elseif(OPENMP_FOUND)
	LIST(APPEND CH_SUPERLUMT_DEFS "__OPENMP")
else()
	MESSAGE(WARNING "No multithreading library has been found")
endif()


# ------------------------------------------------------------------------------
# Dependencies for SuperLUMT module
# ------------------------------------------------------------------------------
message(STATUS "Finding SuperLUMT libraries...")

# Find SuperLUMT
UNSET(SUPERLUMT_FOUND)
set(SUPERLUMT_ROOT "C:/Program Files/SuperLU_MT" CACHE PATH "Location of SuperLUMT installation")
find_package(SUPERLUMT)
message(STATUS "   SUPERLUMT_ROOT is ${SUPERLUMT_ROOT}")
if(SUPERLUMT_FOUND)
	MARK_AS_ADVANCED(FORCE SUPERLUMT_ROOT)
	MARK_AS_ADVANCED(FORCE SUPERLUMT_INCLUDE_DIRS)
	MARK_AS_ADVANCED(FORCE SUPERLUMT_LIBRARIES)
else()
	message(STATUS "Looking for SuperLU_MT manually")
	set(SUPERLUMT_LIBRARIES "C:/Program Files/SuperLU_MT/lib/superlu_mt.lib" CACHE FILE "Location of SuperLU_MT library")
	set(SUPERLUMT_INCLUDE_DIRS "C:/Program Files/SuperLU_MT/include" CACHE PATH "SuperLU_MT include folder")
	find_path(SUPERLUMT_INCLUDE_DIRS supermatrix.h PATHS ${SUPERLUMT_INCLUDE_DIRS} NO_DEFAULT_PATH)
	if(EXISTS "${SUPERLUMT_LIBRARIES}")
		SET(SUPERLUMT_FOUND)
		MARK_AS_ADVANCED(FORCE SUPERLUMT_ROOT)
		MARK_AS_ADVANCED(FORCE SUPERLUMT_INCLUDE_DIRS)
		MARK_AS_ADVANCED(FORCE SUPERLUMT_LIBRARIES)
	endif()
endif()
message(STATUS "   SUPERLUMT include dirs:  ${SUPERLUMT_INCLUDE_DIRS}")
message(STATUS "   SUPERLUMT libraries:     ${SUPERLUMT_LIBRARIES}")


# Find BLAS (eventually looking into SuperLUMT folder)
UNSET(BLAS_FOUND)
UNSET(BLAS_LIBRARIES)
set(BLAS_ROOT "${SUPERLUMT_ROOT}" CACHE PATH "Location of BLAS installation")
find_package(BLAS QUIET)
if(BLAS_FOUND)
	MARK_AS_ADVANCED(FORCE BLAS_ROOT)
	MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES)
	MARK_AS_ADVANCED(FORCE BLAS_INCLUDE_DIRS)
else() # custom, easy, straightforward search of blas library
	find_library(BLAS_LIBRARIES "blas"
		PATHS "${BLAS_ROOT}"
		PATH_SUFFIXES "lib"
		)
		
	find_path(BLAS_INCLUDE_DIRS "f2c.h"
		PATHS "${BLAS_ROOT}"
		PATH_SUFFIXES "include"
		)
	
  if (BLAS_LIBRARIES AND BLAS_INCLUDE_DIRS)
	SET(BLAS_FOUND TRUE)
	MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES)
	MARK_AS_ADVANCED(FORCE BLAS_INCLUDE_DIRS)
  endif (BLAS_LIBRARIES AND BLAS_INCLUDE_DIRS)
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
