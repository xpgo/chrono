############### Headers list ###############
set(ALLHEAD 
    colamd.h
	pxgstrf_synch.h
	slu_dcomplex.h
	slu_mt_Cnames.h
	slu_mt_machines.h
	slu_mt_util.h
	slu_scomplex.h
	supermatrix.h
	)

SET(SLUHEAD
	slu_mt_sdefs.h
)

SET(DLUHEAD
	slu_mt_ddefs.h
)

SET(CLUHEAD
	slu_mt_cdefs.h
)

SET(ZLUHEAD
	slu_mt_zdefs.h
)


############### Sources list ###############
SET(ALLAUX
	await.c
	cholnzcnt.c
	colamd.c
	dclock.c
	get_perm_c.c
	heap_relax_snode.c
	# lsame.c
	mmd.c
	pmemory.c
	pxgstrf_finalize.c
	pxgstrf_mark_busy_descends.c
	pxgstrf_pruneL.c
	pxgstrf_relax_snode.c
	pxgstrf_scheduler.c
	pxgstrf_super_bnd_dfs.c
	pxgstrf_synch.c
	qrnzcnt.c
	sp_coletree.c
	sp_colorder.c
	sp_ienv.c
	superlu_timer.c
	util.c
	# xerbla.c
)

SET(SZAUX
	slamch.c
)

SET(DZAUX
	dlamch.c
)

SET(SLUSRC
	psgssv.c
	psgssvx.c
	psgstrf.c
	psgstrf_bmod1D.c
	psgstrf_bmod1D_mv2.c
	psgstrf_bmod2D.c
	psgstrf_bmod2D_mv2.c
	psgstrf_column_bmod.c
	psgstrf_column_dfs.c
	psgstrf_copy_to_ucol.c
	psgstrf_factor_snode.c
	psgstrf_init.c
	psgstrf_panel_bmod.c
	psgstrf_panel_dfs.c
	psgstrf_pivotL.c
	psgstrf_snode_bmod.c
	psgstrf_snode_dfs.c
	psgstrf_thread.c
	psgstrf_thread_finalize.c
	psgstrf_thread_init.c
	psmemory.c
	psutil.c
	sgscon.c
	sgsequ.c
	sgsrfs.c
	sgstrs.c
	slacon.c
	slangs.c
	slaqgs.c
	smatgen.c
	smyblas2.c
	spivotgrowth.c
	sreadhb.c
	sreadmt.c
	sreadrb.c
	ssp_blas2.c
	ssp_blas3.c
)

SET(DLUSRC
	dgscon.c
	dgsequ.c
	dgsrfs.c
	dgstrs.c
	dlacon.c
	dlangs.c
	dlaqgs.c
	dmatgen.c
	dmyblas2.c
	dpivotgrowth.c
	dreadhb.c
	dreadmt.c
	dreadrb.c
	dsp_blas2.c
	dsp_blas3.c
	pdgssv.c
	pdgssvx.c
	pdgstrf.c
	pdgstrf_bmod1D.c
	pdgstrf_bmod1D_mv2.c
	pdgstrf_bmod2D.c
	pdgstrf_bmod2D_mv2.c
	pdgstrf_column_bmod.c
	pdgstrf_column_dfs.c
	pdgstrf_copy_to_ucol.c
	pdgstrf_factor_snode.c
	pdgstrf_init.c
	pdgstrf_panel_bmod.c
	pdgstrf_panel_dfs.c
	pdgstrf_pivotL.c
	pdgstrf_snode_bmod.c
	pdgstrf_snode_dfs.c
	pdgstrf_thread.c
	pdgstrf_thread_finalize.c
	pdgstrf_thread_init.c
	pdmemory.c
	pdutil.c
)

SET(CLUSRC
	cgscon.c
	cgsequ.c
	cgsrfs.c
	cgstrs.c
	clacon.c
	clangs.c
	claqgs.c
	cmatgen.c
	cmyblas2.c
	cpivotgrowth.c
	creadhb.c
	creadmt.c
	creadrb.c
	csp_blas2.c
	csp_blas3.c
	icmax1.c
	pcgssv.c
	pcgssvx.c
	pcgstrf.c
	pcgstrf_bmod1D.c
	pcgstrf_bmod1D_mv2.c
	pcgstrf_bmod2D.c
	pcgstrf_bmod2D_mv2.c
	pcgstrf_column_bmod.c
	pcgstrf_column_dfs.c
	pcgstrf_copy_to_ucol.c
	pcgstrf_factor_snode.c
	pcgstrf_init.c
	pcgstrf_panel_bmod.c
	pcgstrf_panel_dfs.c
	pcgstrf_pivotL.c
	pcgstrf_snode_bmod.c
	pcgstrf_snode_dfs.c
	pcgstrf_thread.c
	pcgstrf_thread_finalize.c
	pcgstrf_thread_init.c
	pcmemory.c
	pcutil.c
	scomplex.c
	scsum1.c
)

SET(ZLUSRC
	dcomplex.c
	dzsum1.c
	izmax1.c
	pzgssv.c
	pzgssvx.c
	pzgstrf.c
	pzgstrf_bmod1D.c
	pzgstrf_bmod1D_mv2.c
	pzgstrf_bmod2D.c
	pzgstrf_bmod2D_mv2.c
	pzgstrf_column_bmod.c
	pzgstrf_column_dfs.c
	pzgstrf_copy_to_ucol.c
	pzgstrf_factor_snode.c
	pzgstrf_init.c
	pzgstrf_panel_bmod.c
	pzgstrf_panel_dfs.c
	pzgstrf_pivotL.c
	pzgstrf_snode_bmod.c
	pzgstrf_snode_dfs.c
	pzgstrf_thread.c
	pzgstrf_thread_finalize.c
	pzgstrf_thread_init.c
	pzmemory.c
	pzutil.c
	zgscon.c
	zgsequ.c
	zgsrfs.c
	zgstrs.c
	zlacon.c
	zlangs.c
	zlaqgs.c
	zmatgen.c
	zmyblas2.c
	zpivotgrowth.c
	zreadhb.c
	zreadmt.c
	zreadrb.c
	zsp_blas2.c
	zsp_blas3.c
)


############### Specific file compilation options ###############

set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(slamch.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(dlamch.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(dclock.c PROPERTIES COMPILE_FLAGS -O0)

# TODO
## Do not optimize this on Cray; related to 'volatile' variable.
# await.c
# : await.c
# $(CC) -c $(NOOPTS) $< $(VERBOSE)

######### Sources and headers groups #########
SET(LUSRC "${ALLAUX}") # all sources that will be compiled must go here
SET(LUHEAD "${ALLHEAD}") # all headers that will be included must go here

if(enable_single)
  LIST(APPEND LUSRC "${SLUSRC}")
  LIST(APPEND LUHEAD "${SLUHEAD}")
endif(enable_single)

if(enable_double)
  LIST(APPEND LUSRC "${DLUSRC}" )
  LIST(APPEND LUHEAD "${DLUHEAD}")
endif(enable_double)

if(enable_complex)
  LIST(APPEND LUSRC "${CLUSRC}" )
  LIST(APPEND LUHEAD "${CLUHEAD}")
endif(enable_complex)

if(enable_complex16)
  LIST(APPEND LUSRC "${ZLUSRC}" )
  LIST(APPEND LUHEAD "${ZLUHEAD}")
endif(enable_complex16)


if(enable_complex16 OR enable_double)
  LIST(APPEND LUSRC "${DZAUX}")
endif(enable_complex16 OR enable_double)

if(enable_complex OR enable_single)
  LIST(APPEND LUSRC "${SZAUX}")
endif(enable_complex OR enable_single)


############### Build library ############### 
add_library(superlu_mt ${LUSRC} ${LUHEAD})

include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(superlu_mt ${BLAS_LIBRARIES} ${MULTITHREADING_LIBRARIES})

if (NOT WIN32)
	target_link_libraries(superlu_mt m)
endif()

set_target_properties(superlu_mt PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${VERSION_MAJOR}
  COMPILE_DEFINITIONS "${SUPERLU_DEFS}"
  )

if (CMAKE_COMPILER_IS_GNUC OR CMAKE_COMPILER_IS_GNUCXX)
	include(GNUInstallDirs)
endif()

############### Install library ###############
install(TARGETS superlu_mt
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES ${LUHEAD}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
