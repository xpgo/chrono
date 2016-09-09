#--------------------------------------------------------------
# Add executables

SET(TESTS
  test_SUPERLUMT_Engine
)

MESSAGE(STATUS "Additional test programs for Chrono_SuperLUMT module...")

INCLUDE_DIRECTORIES("${CH_SUPERLUMT_INCLUDES}")

FOREACH(PROGRAM ${TESTS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER tests
		COMPILE_DEFINITIONS "${CH_SUPERLUMT_DEFS}"
        COMPILE_FLAGS "${CH_CXX_FLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ChronoEngine ChronoEngine_SuperLUMT)

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)
