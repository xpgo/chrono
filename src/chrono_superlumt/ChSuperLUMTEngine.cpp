// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni, Radu Serban
// =============================================================================
// Interfacing to the Pardiso Sparse Direct Solver from the Intel® MKL Library.
// =============================================================================

#include "chrono_superlumt/ChSuperLUMTEngine.h"
#include <omp.h>
#include <algorithm>


namespace chrono {


	ChSuperLUMTEngine::ChSuperLUMTEngine()
	{

		ResetSolver();

		panel_size = sp_ienv(1);
		relax = sp_ienv(2);

		dCreate_Dense_Matrix(&m_sol_Super, 0, 0, nullptr, 0, SLU_DN, SLU_D, SLU_GE);
		dCreate_Dense_Matrix(&m_rhs_Super, 0, 0, nullptr, 0, SLU_DN, SLU_D, SLU_GE);
		dCreate_CompCol_Matrix(&m_mat_Super, 0, 0, 0, nullptr, nullptr, nullptr, SLU_NC, SLU_D, SLU_GE);

		StatAlloc(m_n, superlumt_options.nprocs, panel_size, relax, &Gstat);

	}

	ChSuperLUMTEngine::~ChSuperLUMTEngine()
	{
		//Destroy_CompCol_Matrix(&m_mat_Super); this would destroy also the CSR arrays
		Destroy_SuperMatrix_Store(&m_mat_Super);
		Destroy_SuperMatrix_Store(&m_rhs_Super);
		Destroy_SuperMatrix_Store(&m_sol_Super);
		if (lwork == 0) {
			if (L.nrow==m_n && L.ncol == m_n && m_n!=0) // checks if the solver has ever been called
			{
				Destroy_SuperNode_Matrix(&L);
				Destroy_CompCol_Matrix(&U);
			}
		}
	}

	void ChSuperLUMTEngine::SetMatrix(ChSparseMatrix& Z)
	{

		SetMatrix(Z.GetNumRows(), Z.GetCSR_ValueArray(), Z.GetCSR_TrailingIndexArray(), Z.GetCSR_LeadingIndexArray());
		//TODO: change here the matrix type (symmetric, etc) according to what Z says
	}

	void ChSuperLUMTEngine::SetMatrix(int pb_size, double* values, int* rowIndex, int* colIndex)
	{
		m_n = pb_size;
		m_mat_Super.nrow = m_n;
		m_mat_Super.ncol = m_n;
		auto m_mat_Super_store = static_cast<NCformat*>(m_mat_Super.Store);
		m_mat_Super_store->nnz = colIndex[pb_size];
		m_mat_Super_store->nzval = values;
		m_mat_Super_store->rowind = rowIndex;
		m_mat_Super_store->colptr = colIndex;

		R.resize(m_mat_Super.nrow);
		C.resize(m_mat_Super.ncol);

		perm_r.resize(m_n);
		perm_c.resize(m_n);

		get_perm_c(MMD_ATA, &m_mat_Super, perm_c.data());

		// Internally computed
		etree.resize(m_n);
		colcnt_h.resize(m_n);
		part_super_h.resize(m_n);

		superlumt_options.etree = etree.data();
		superlumt_options.colcnt_h = colcnt_h.data();
		superlumt_options.part_super_h = part_super_h.data();

		// update lda values also for rhs and sol because they could be
		// not initialize if a factorization is called before
		// initializing the sol and rhs
		static_cast<DNformat*>(m_rhs_Super.Store)->lda = m_n;
		static_cast<DNformat*>(m_sol_Super.Store)->lda = m_n;
	}

	void ChSuperLUMTEngine::SetSolutionVector(ChMatrix<>& x)
	{
		assert(x.GetRows() == m_n);

		SetSolutionVector(x.GetAddress());
	}

	void ChSuperLUMTEngine::SetSolutionVector(double* x)
	{
		m_sol_Super.nrow = m_n;
		m_sol_Super.ncol = 1;
		auto m_sol_Super_store = static_cast<DNformat*>(m_sol_Super.Store);
		//m_sol_Super_store->lda = m_n; // done during the SuperLU_MT call
		m_sol_Super_store->nzval = x;

		ferr.resize(m_nrhs);
	}

	void ChSuperLUMTEngine::SetRhsVector(ChMatrix<>& b)
	{
		assert(b.GetRows() == m_n);

		// that should deal multiple rhs
		if (b.GetColumns()>1)
			b.MatrTranspose(); // ChMatrix is row-major

		SetRhsVector(b.GetAddress(), b.GetColumns());
	}

	void ChSuperLUMTEngine::SetRhsVector(double* b, int nrhs)
	{
		m_nrhs = nrhs;

		m_rhs_Super.nrow = m_n;
		m_rhs_Super.ncol = nrhs;
		auto m_rhs_Super_store = static_cast<DNformat*>(m_rhs_Super.Store);
		//m_rhs_Super_store->lda = m_n; // done during the SuperLU_MT call
		m_rhs_Super_store->nzval = b;

		berr.resize(m_nrhs);
	}

	void ChSuperLUMTEngine::SetProblem(ChSparseMatrix& Z, ChMatrix<>& b, ChMatrix<>& x)
	{
		SetMatrix(Z);
		SetRhsVector(b);
        SetSolutionVector(x);
    }

	int ChSuperLUMTEngine::SuperLUMTCall(phase_t phase, int verbose)
	{
		auto m_rhs_Super_ncol_bkp = m_rhs_Super.ncol;
		auto m_sol_Super_ncol_bkp = m_sol_Super.ncol;

        // Set the number of right-hand side
        // if put to 0 then factorize without solve
        m_rhs_Super.ncol = (phase == phase_t::ANALYSIS_NUMFACTORIZATION) ? 0 : m_nrhs;
		m_sol_Super.ncol = (phase == phase_t::ANALYSIS_NUMFACTORIZATION) ? 0 : m_nrhs;
		superlumt_options.fact = (phase == phase_t::SOLVE) ? FACTORED : EQUILIBRATE; /* Indicate the factored form of m_mat_Super is supplied. */
		//TODO: superlumt_options.refact


		// call to SuperLU_MT modified routine
		pdgssvx_mod(superlumt_options.nprocs, &superlumt_options, &m_mat_Super, perm_c.data(), perm_r.data(), &equed, R.data(), C.data(),
			&L, &U, &m_rhs_Super, &m_sol_Super, &rpg, &m_rcond, ferr.data(), berr.data(),
			&superlu_memusage, &info);


		// restore previous values of columns of rhs and sol
		// otherwise if someone set rhs and sol before factorization they will be corrupted
		m_rhs_Super.ncol = m_rhs_Super_ncol_bkp;
		m_sol_Super.ncol = m_sol_Super_ncol_bkp;

		if (verbose>1)
		{
			if (phase == phase_t::ANALYSIS_NUMFACTORIZATION || phase == phase_t::COMPLETE)
			{
				if (info == 0 || info == m_n + 1)
				{
					printf("Recip. pivot growth = %e\n", rpg);
					if (rcond_evaluation) printf("Recip. condition number = %e\n", m_rcond);
					printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
					for (auto i = 0; i < m_nrhs; ++i) {
						printf(IFMT "%16e%16e\n", i + 1, ferr[i], berr[i]);
					}

					auto Lnnz = static_cast<SCPformat*>(L.Store)->nnz;
					auto Unnz = static_cast<NCPformat*>(U.Store)->nnz;
					printf("No of nonzeros in factor L = " IFMT "\n", Lnnz);
					printf("No of nonzeros in factor U = " IFMT "\n", Unnz);
					printf("No of nonzeros in L+U = " IFMT "\n", Lnnz + Unnz - m_n);
					printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
						superlu_memusage.for_lu / 1e6, superlu_memusage.total_needed / 1e6,
						superlu_memusage.expansions);

				}
				else if (info > 0 && lwork == -1)
					printf("** Estimated memory: %d bytes\n", info - m_n);
			}

			if (phase == phase_t::SOLVE || phase == phase_t::COMPLETE)
			{
				if (info == 0 || info == m_n + 1)
					printf("Triangular solve: dgssvx() returns info %d\n", info);
				else if (info > 0 && lwork == -1)
					printf("** Estimated memory: " IFMT " bytes\n", info - m_n);
			}
			
			printf("time for [EQUIL] phase: %f\n", Gstat.utime[EQUIL]);
			printf("time for [FACT] phase: %f\n", Gstat.utime[FACT]);
			printf("time for [RCOND] phase: %f\n", Gstat.utime[RCOND]);
			printf("time for [SOLVE] phase: %f\n", Gstat.utime[SOLVE]);
			printf("time for [REFINE] phase: %f\n", Gstat.utime[REFINE]);

		}

		return info;
	}

	void ChSuperLUMTEngine::ResetSolver()
	{
		superlumt_options.nprocs = omp_get_num_procs();
		superlumt_options.fact = EQUILIBRATE;
		superlumt_options.trans = NOTRANS;
		superlumt_options.refact = NO;
		superlumt_options.panel_size = panel_size;
		superlumt_options.relax = relax;
		superlumt_options.usepr = NO;
		superlumt_options.drop_tol = 0.0;
		superlumt_options.diag_pivot_thresh = 1.0;
		superlumt_options.SymmetricMode = NO;
		superlumt_options.PrintStat = NO;
		superlumt_options.perm_c = perm_c.data();
		superlumt_options.perm_r = perm_r.data();
		superlumt_options.work = work;
		superlumt_options.lwork = lwork;
	}

	void ChSuperLUMTEngine::SetNumProcs(int nprocs_in)
	{
		StatFree(&Gstat);
		superlumt_options.nprocs = std::max(nprocs_in,omp_get_num_procs());
		StatAlloc(m_n, superlumt_options.nprocs, panel_size, relax, &Gstat);
	}

	void ChSuperLUMTEngine::GetResidual(ChMatrix<>& res) const
	{
		assert(res.GetRows() >= m_n);
		GetResidual(res.GetAddress());
	}

	void ChSuperLUMTEngine::GetResidual(double* res) const
	{
		auto m_mat_Super_store = static_cast<NCformat*>(m_mat_Super.Store);
		auto m_a = static_cast<double*>(m_mat_Super_store->nzval);
		auto m_ia = static_cast<int*>(m_mat_Super_store->rowind);
		auto m_ja = static_cast<int*>(m_mat_Super_store->colptr);

		auto m_b = static_cast<double*>(static_cast<DNformat*>(m_rhs_Super.Store)->nzval);
		auto m_x = static_cast<double*>(static_cast<DNformat*>(m_sol_Super.Store)->nzval);


		for (auto col_sel = 0; col_sel<m_n; ++col_sel)
		{
			for (auto row_sel = m_ja[col_sel]; row_sel<m_ja[col_sel+1]; ++row_sel)
			{
				res[m_ia[col_sel]] = m_b[col_sel] - m_a[row_sel] * m_x[col_sel];
			}
		}
	}

	double ChSuperLUMTEngine::GetResidualNorm() const
	{
		std::vector<double> res(m_n);
		GetResidual(res.data());
		double norm = 0;
		for (auto i = 0; i < m_n; i++) {
			norm += res[i] * res[i];
		}

		// also BLAS function can be used here
		return std::sqrt(norm);
	}

	void ChSuperLUMTEngine::pdgssvx_mod(int_t nprocs, superlumt_options_t* superlumt_options, SuperMatrix* A, int_t* perm_c, int_t* perm_r, equed_t* equed, double* R, double* C, SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperMatrix* X, double* recip_pivot_growth, double* rcond, double* ferr, double* berr, superlu_memusage_t* superlu_memusage, int_t* info)
	{
		DNformat  *Bstore = static_cast<DNformat*>(B->Store);
		DNformat  *Xstore = static_cast<DNformat*>(X->Store);
		double    *Bval = static_cast<double*>(Bstore->nzval);
		double    *Xval = static_cast<double*>(Xstore->nzval);
		int_t     ldb = Bstore->lda;
		int_t     ldx = Bstore->lda;
		int_t     nrhs = B->ncol;
		SuperMatrix *AA; /* A in NC format used by the factorization routine.*/
		SuperMatrix AC; /* Matrix postmultiplied by Pc */
		bool       colequ, rowequ;
		trans_t   trant;
		int_t     info1;
		double amax, anorm, bignum, smlnum, colcnd, rowcnd, rcmax, rcmin;
		int_t     panel_size = superlumt_options->panel_size;
		double    t0;      /* temporary time */
		

		/* External functions */
		//extern double dlangs(char *, SuperMatrix *);
		//extern double dlamch_(char *);

		superlumt_options->perm_c = perm_c;
		superlumt_options->perm_r = perm_r;


		*info = 0;
		auto dofact = (superlumt_options->fact == DOFACT);
		auto equil = (superlumt_options->fact == EQUILIBRATE);
		auto notran = (superlumt_options->trans == NOTRANS);

		if (dofact || equil) {
			*equed = NOEQUIL;
			rowequ = FALSE;
			colequ = FALSE;
		}
		else {
			rowequ = (*equed == ROW) || (*equed == BOTH);
			colequ = (*equed == COL) || (*equed == BOTH);
		}
		smlnum = dlamch_("Safe minimum");
		bignum = 1. / smlnum;

		/* ------------------------------------------------------------
		Test the input parameters.
		------------------------------------------------------------*/
		if (nprocs <= 0) *info = -1;
		else if ((!dofact && !equil && superlumt_options->fact != FACTORED)
			|| (!notran && (superlumt_options->trans != TRANS) && (superlumt_options->trans != CONJ))
			|| (superlumt_options->refact != YES && superlumt_options->refact != NO)
			|| (superlumt_options->usepr != YES && superlumt_options->usepr != NO)
			|| superlumt_options->lwork < -1)
			*info = -2;
		else if (A->nrow != A->ncol || A->nrow < 0 ||
			    (A->Stype != SLU_NC && A->Stype != SLU_NR) ||
			     A->Dtype != SLU_D || A->Mtype != SLU_GE)
			*info = -3;
		else if (superlumt_options->fact == FACTORED && !(rowequ || colequ || *equed == NOEQUIL))
			*info = -6;
		else {
			if (rowequ) {
				rcmin = bignum;
				rcmax = 0.;
				for (auto j = 0; j < A->nrow; ++j) {
					rcmin = SUPERLU_MIN(rcmin, R[j]);
					rcmax = SUPERLU_MAX(rcmax, R[j]);
				}
				if (rcmin <= 0.) *info = -7;
				else if (A->nrow > 0)
					rowcnd = SUPERLU_MAX(rcmin, smlnum) / SUPERLU_MIN(rcmax, bignum);
				else rowcnd = 1.;
			}
			if (colequ && *info == 0) {
				rcmin = bignum;
				rcmax = 0.;
				for (auto j = 0; j < A->nrow; ++j) {
					rcmin = SUPERLU_MIN(rcmin, C[j]);
					rcmax = SUPERLU_MAX(rcmax, C[j]);
				}
				if (rcmin <= 0.) *info = -8;
				else if (A->nrow > 0)
					colcnd = SUPERLU_MAX(rcmin, smlnum) / SUPERLU_MIN(rcmax, bignum);
				else colcnd = 1.;
			}
			if (*info == 0) {
				if (B->ncol < 0 || Bstore->lda < SUPERLU_MAX(0, A->nrow) ||
					B->Stype != SLU_DN || B->Dtype != SLU_D ||
					B->Mtype != SLU_GE)
					*info = -11;
				else if (X->ncol < 0 || Xstore->lda < SUPERLU_MAX(0, A->nrow) ||
					B->ncol != X->ncol || X->Stype != SLU_DN ||
					X->Dtype != SLU_D || X->Mtype != SLU_GE)
					*info = -12;
			}
		}

		if (*info != 0) {
			int_t i = -(*info);
			xerbla_("pdgssvx", &i);
			return;
		}


		/* ------------------------------------------------------------
		Allocate storage and initialize statistics variables.
		------------------------------------------------------------*/
		StatInit(m_n, superlumt_options->nprocs, &Gstat);

		/* ------------------------------------------------------------
		Convert A to NC format when necessary.
		------------------------------------------------------------*/
		if (A->Stype == SLU_NR) {
			NRformat *AstoreNR = static_cast<NRformat*>(A->Store);
			AA = static_cast<SuperMatrix *>(SUPERLU_MALLOC(sizeof(SuperMatrix)));
			dCreate_CompCol_Matrix(AA, A->ncol, A->nrow, AstoreNR->nnz,
				static_cast<double*>(AstoreNR->nzval), AstoreNR->colind, AstoreNR->rowptr,
				SLU_NC, A->Dtype, A->Mtype);
			if (notran) { /* Reverse the transpose argument. */
				trant = TRANS;
				notran = 0;
			}
			else {
				trant = NOTRANS;
				notran = 1;
			}
		}
		else { /* A->Stype == NC */
			trant = superlumt_options->trans;
			AA = A;
		}

		/* ------------------------------------------------------------
		Diagonal scaling to equilibrate the matrix.
		------------------------------------------------------------*/
		if (equil) {
			t0 = SuperLU_timer_();
			/* Compute row and column scalings to equilibrate the matrix A. */
			dgsequ(AA, R, C, &rowcnd, &colcnd, &amax, &info1);

			if (info1 == 0) {
				/* Equilibrate matrix A. */
				dlaqgs(AA, R, C, rowcnd, colcnd, amax, equed);
				rowequ = (*equed == ROW) || (*equed == BOTH);
				colequ = (*equed == COL) || (*equed == BOTH);
			}
            Gstat.utime[EQUIL] = SuperLU_timer_() - t0;
		}

		/* ------------------------------------------------------------
		Scale the right hand side.
		------------------------------------------------------------*/
		if (notran) {
			if (rowequ) {
				for (auto j = 0; j < nrhs; ++j)
					for (auto i = 0; i < A->nrow; ++i) {
						Bval[i + j*ldb] *= R[i];
					}
			}
		}
		else if (colequ) {
			for (auto j = 0; j < nrhs; ++j)
				for (auto i = 0; i < A->nrow; ++i) {
					Bval[i + j*ldb] *= C[i];
				}
		}

		/* ------------------------------------------------------------
		Perform the LU factorization.
		------------------------------------------------------------*/
		if (dofact || equil) {
            /* Obtain column etree, the column count (colcnt_h) and supernode
            partition (part_super_h) for the Householder matrix. */
            t0 = SuperLU_timer_();
			sp_colorder(AA, perm_c, superlumt_options, &AC);
            Gstat.utime[ETREE] = SuperLU_timer_() - t0;

#if (PRNTlevel>=2)
			printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n",
				superlumt_options->relax, panel_size, sp_ienv(3), sp_ienv(4));
			fflush(stdout);
#endif

			/* Compute the LU factorization of A*Pc. */
			t0 = SuperLU_timer_();
			pdgstrf(superlumt_options, &AC, perm_r, L, U, &Gstat, info);
            Gstat.utime[FACT] = SuperLU_timer_() - t0;

			flops_t flopcnt = 0;
			for (auto i = 0; i < nprocs; ++i) flopcnt += Gstat.procstat[i].fcops;
            Gstat.ops[FACT] = flopcnt;

			if (superlumt_options->lwork == -1) {
				superlu_memusage->total_needed = static_cast<float_t>(*info - A->ncol);
				return;
			}
		}

		// if factorization errors occurred...
		if (*info > 0) { 
			if (*info <= A->ncol) {
				/* Compute the reciprocal pivot growth factor of the leading
				rank-deficient *info columns of A. */
				*recip_pivot_growth = dPivotGrowth(*info, AA, perm_c, L, U);
			}
		}
		// if no errors in factorization...
		else {

			/* ------------------------------------------------------------
			Compute the reciprocal pivot growth factor *recip_pivot_growth.
			------------------------------------------------------------*/
			*recip_pivot_growth = dPivotGrowth(A->ncol, AA, perm_c, L, U);

			if (rcond_evaluation)
			{
				/* ------------------------------------------------------------
				Estimate the reciprocal of the condition number of A.
				------------------------------------------------------------*/
				char      norm[1];
				t0 = SuperLU_timer_();
				norm[0] = notran ? '1' : 'I';
				anorm = dlangs(norm, AA);
				dgscon(norm, L, U, anorm, rcond, info);
                Gstat.utime[RCOND] = SuperLU_timer_() - t0;
			}

			// if there is actually a rhs
			if (B->ncol>0)
			{
				/* ------------------------------------------------------------
				Compute the solution matrix X.
				------------------------------------------------------------*/
				for (auto j = 0; j < nrhs; j++)    /* Save a copy of the right hand sides */
					for (auto i = 0; i < B->nrow; i++)
						Xval[i + j*ldx] = Bval[i + j*ldb];

				t0 = SuperLU_timer_();
				dgstrs(trant, L, U, perm_r, perm_c, X, &Gstat, info);
				Gstat.utime[SOLVE] = SuperLU_timer_() - t0;
				Gstat.ops[SOLVE] = Gstat.ops[TRISOLVE];

				if (iterative_refinement)
				{
					/* ------------------------------------------------------------
					Use iterative refinement to improve the computed solution and
					compute error bounds and backward error estimates for it.
					------------------------------------------------------------*/
					t0 = SuperLU_timer_();
					dgsrfs(trant, AA, L, U, perm_r, perm_c, *equed, R, C, B, X, ferr, berr, &Gstat, info);
					Gstat.utime[REFINE] = SuperLU_timer_() - t0;
				}


				/* ------------------------------------------------------------
				Transform the solution matrix X to a solution of the original
				system.
				------------------------------------------------------------*/
				if (notran) {
					if (colequ) {
						for (auto j = 0; j < nrhs; ++j)
							for (auto i = 0; i < A->nrow; ++i) {
								Xval[i + j*ldx] *= C[i];
							}
					}
				}
				else if (rowequ) {
					for (auto j = 0; j < nrhs; ++j)
						for (auto i = 0; i < A->nrow; ++i) {
							Xval[i + j*ldx] *= R[i];
						}
				}

				if (rcond_evaluation)
				{
					/* Set INFO = A->ncol+1 if the matrix is singular to
					working precision.*/
					if (*rcond < dlamch_("E"))
						*info = A->ncol + 1;
				}
			} // end if (rhs has non-zero dimension)
		}
		superlu_dQuerySpace(nprocs, L, U, panel_size, superlu_memusage);

        /* ------------------------------------------------------------
		Deallocate storage after factorization.
		------------------------------------------------------------*/
		if (dofact || equil) {
			Destroy_CompCol_Permuted(&AC);
		}
		if (A->Stype == SLU_NR) {
			Destroy_SuperMatrix_Store(AA);
			SUPERLU_FREE(AA);
		}

		/* ------------------------------------------------------------
		Print timings, then deallocate statistic variables.
		------------------------------------------------------------*/
#ifdef PROFILE
		{
			SCPformat *Lstore = (SCPformat *)L->Store;
			ParallelProfile(n, Lstore->nsuper + 1, Gstat.num_panels, nprocs, &Gstat);
		}
		printf("colcnt_h[0] %lld\n", superlumt_options->colcnt_h[0]);
		PrintStat(&Gstat);
#endif



	}


	
}  // end namespace chrono
