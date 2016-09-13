SET(SEXAMPLE
pslinsol
pslinsolx
pslinsolx1
pslinsolx2
psrepeat
psspmd
)


SET(DEXAMPLE
pdlinsol
pdlinsolx
pdlinsolx1
pdlinsolx2
pdrepeat
pdspmd
)

SET(CEXAMPLE
pclinsol
pclinsolx
pclinsolx1
pclinsolx2
pcrepeat
pcspmd
)

SET(ZEXAMPLE
pzlinsol
pzlinsolx
pzlinsolx1
pzlinsolx2
pzrepeat
pzspmd
)

SET (ACTIVE_TARGET)
if(enable_single)
	LIST(APPEND ACTIVE_TARGET "${SEXAMPLE}")
endif(enable_single)

if(enable_double)
	LIST(APPEND ACTIVE_TARGET "${DEXAMPLE}")
endif(enable_double)

if(enable_complex)
	LIST(APPEND ACTIVE_TARGET "${CEXAMPLE}")
endif(enable_complex)

if(enable_complex16)
	LIST(APPEND ACTIVE_TARGET "${ZEXAMPLE}")
endif(enable_complex16)

if (WIN32)
	SET(GETOPT_HEAD getopt.h)
endif()

MESSAGE(STATUS "Build examples...")
foreach(PROGRAM ${ACTIVE_TARGET})
    message(STATUS "  ...add ${PROGRAM}")

    add_executable(${PROGRAM}  "${PROGRAM}.c" ${GETOPT_HEAD})
	
	set_target_properties(${PROGRAM} PROPERTIES
		  COMPILE_DEFINITIONS "${SUPERLU_DEFS}"
		  )
	
    target_link_libraries(${PROGRAM} superlu_mt ${BLAS_LIBRARIES})

endforeach()