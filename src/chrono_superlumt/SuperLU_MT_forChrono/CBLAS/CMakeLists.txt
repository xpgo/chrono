set(BLAS_HEAD
	superlu_f2c.h
    f2c.h
    slu_mt_Cnames.h
)

set(BLAS_SRC_COMMON
	../SRC/lsame.c
	../SRC/xerbla.c
)

SET(SBLAS1
	isamax.c
	sasum.c
	saxpy.c
	scopy.c
	sdot.c
	snrm2.c
	srot.c
	sscal.c
)

SET(SBLAS2
	sgemv.c
	sger.c
	ssymv.c
	ssyr2.c
	strsv.c
)

SET(DBLAS1
	idamax.c
	dasum.c
	daxpy.c
	dcopy.c
	ddot.c
	dnrm2.c
	drot.c
	dscal.c
)

SET(DBLAS2
	dgemv.c
	dger.c
	dsymv.c
	dsyr2.c
	dtrsv.c
)

SET(CBLAS1
	caxpy.c
	ccopy.c
	cscal.c
	icamax.c
	scasum.c
	scnrm2.c
)

SET(CBLAS2
	cgemv.c
	cgerc.c
	chemv.c
	cher2.c
	ctrsv.c
)

SET(ZBLAS1
	dcabs1.c
	dzasum.c
	dznrm2.c
	izamax.c
	zaxpy.c
	zcopy.c
	zscal.c
)

SET(ZBLAS2
	zgemv.c
	zgerc.c
	zhemv.c
	zher2.c
	ztrsv.c
)

# TODO: Are there xBLAS3 somewhere???


######### Sources and headers groups #########
SET(BLAS_SRC "${BLAS_SRC_COMMON}") # all sources that will be compiled must go here

if(enable_single)
  LIST(APPEND BLAS_SRC "${SBLAS1}" "${SBLAS2}")
endif(enable_single)

if(enable_double)
  LIST(APPEND BLAS_SRC "${DBLAS1}" "${DBLAS2}")
endif(enable_double)

if(enable_complex)
  LIST(APPEND BLAS_SRC "${CBLAS1}" "${CBLAS2}")
endif(enable_complex)

if(enable_complex16)
  LIST(APPEND BLAS_SRC "${ZBLAS1}" "${ZBLAS2}")
endif(enable_complex16)




############### Build library ###############
add_library(blas ${BLAS_SRC} ${BLAS_HEAD})
target_link_libraries(blas ${MULTITHREADING_LIBRARIES})
# include (GenerateExportHeader)

# GENERATE_EXPORT_HEADER(blas
             # BASE_NAME blas
             # EXPORT_MACRO_NAME blas_EXPORT
             # EXPORT_FILE_NAME blas_Export.h
             # STATIC_DEFINE blas_BUILT_AS_STATIC
# )

############### Install library ###############
install(TARGETS blas
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES ${BLAS_HEAD}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

############### Export variables as FindBLAS would do ###############
SET(BLAS_LIBRARIES "blas" PARENT_SCOPE)
SET(BLAS_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
