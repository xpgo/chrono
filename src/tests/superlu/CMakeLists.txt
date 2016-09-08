#--------------------------------------------------------------
# Add executables

SET(TESTS
  test_SUPERLU_Engine
)

MESSAGE(STATUS "Additional test programs for Chrono_SuperLU module...")

INCLUDE_DIRECTORIES("${CH_SUPERLU_INCLUDES}")

FOREACH(PROGRAM ${TESTS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER tests
		COMPILE_DEFINITIONS "${CH_SUPERLU_DEFS}"
        COMPILE_FLAGS "${CH_CXX_FLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ChronoEngine ChronoEngine_SuperLU)

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)
