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
# Dependencies for SuperLU module
# ------------------------------------------------------------------------------

message(STATUS "Finding SuperLU libraries...")
find_package(SUPERLU)
message(STATUS "SUPERLU_ROOT is ${SUPERLU_ROOT}")

message(STATUS "   SUPERLU include dirs:  ${SUPERLU_INCLUDE_DIRS}")
message(STATUS "   SUPERLU libraries:     ${SUPERLU_LIBRARIES}")
message(STATUS "   BLAS include dirs:     ${BLAS_INCLUDE_DIRS}")
message(STATUS "   BLAS library:          ${BLAS_LIBRARIES}")

if(NOT SUPERLU_LIBRARIES)
	set(SUPERLU_ROOT "C:/Program Files/SuperLU" CACHE PATH "Location of SuperLU installation")
else()
	MARK_AS_ADVANCED(FORCE SUPERLU_ROOT)
endif()

if(NOT BLAS_LIBRARIES)
	set(BLAS_ROOT "${SUPERLU_ROOT}" CACHE PATH "Location of SuperLU installation")
else()
	MARK_AS_ADVANCED(FORCE BLAS_ROOT)
endif()

# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the SUPERLU module
# ------------------------------------------------------------------------------

SET(CH_SUPERLU_INCLUDES ${SUPERLU_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CH_SUPERLU_INCLUDES})
SET(CH_SUPERLU_INCLUDES "${CH_SUPERLU_INCLUDES}" PARENT_SCOPE)

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
    # ChSolverSuperLU.cpp
)

SOURCE_GROUP("Headers" FILES ${ChronoEngine_SuperLU_HEADERS})
SOURCE_GROUP("Sources" FILES ${ChronoEngine_SuperLU_SOURCES})


# ------------------------------------------------------------------------------
# Add the ChronoEngine_SuperLU library
# ------------------------------------------------------------------------------

ADD_LIBRARY(ChronoEngine_SuperLU SHARED
            ${ChronoEngine_SuperLU_SOURCES}
            ${ChronoEngine_SuperLU_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_SuperLU PROPERTIES
                      COMPILE_FLAGS "${COMPILER_FLAGS}"
                      LINK_FLAGS "${LINKER_FLAGS}"
                      COMPILE_DEFINITIONS "CH_API_COMPILE_SUPERLU")

TARGET_LINK_LIBRARIES(ChronoEngine_SuperLU
                      ChronoEngine
                      ${SUPERLU_LIBRARIES}
					  ${BLAS_LIBRARIES}
					  )

INSTALL(TARGETS ChronoEngine_SuperLU
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

INSTALL(FILES ${ChronoEngine_SuperLU_HEADERS} 
        DESTINATION include/chrono_superlu)
