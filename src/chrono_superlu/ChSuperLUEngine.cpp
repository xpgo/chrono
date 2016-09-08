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

#include <algorithm>

#include "chrono_superlu/ChSuperLUEngine.h"

namespace chrono {
	enum phase_t {
		COMPLETE = 13,
		ANALYSIS_NUMFACTORIZATION = 12,
		SOLVE = 33,
	};

	ChSuperLUEngine::ChSuperLUEngine()
	{
		lwork = 0;
		panel_size = sp_ienv(1);
		relax = sp_ienv(2);

		ResetSolver();

		// it is fundamental to leave m_sol_Super.ncol=0, and m_rhs_Super.ncol=0
		// in order to avoid checking of dgssvx routine
		dCreate_Dense_Matrix(&m_sol_Super, 0, 0, nullptr, 0, SLU_DN, SLU_D, SLU_GE);
		dCreate_Dense_Matrix(&m_rhs_Super, 0, 0, nullptr, 0, SLU_DN, SLU_D, SLU_GE);
		dCreate_CompCol_Matrix(&m_mat_Super, 0, 0, 0, nullptr, nullptr, nullptr, SLU_NC, SLU_D, SLU_GE);
	}

	ChSuperLUEngine::~ChSuperLUEngine()
	{
		
		//Destroy_CompCol_Matrix(&m_mat_Super); this would destroy also the CSR arrays
		Destroy_SuperMatrix_Store(&m_mat_Super);
		Destroy_SuperMatrix_Store(&m_rhs_Super);
		Destroy_SuperMatrix_Store(&m_sol_Super);
		if (lwork == 0) {
			if (L.nrow==m_n && L.ncol == m_n) // checks if the solver has ever been called
			{
				Destroy_SuperNode_Matrix(&L);
				Destroy_CompCol_Matrix(&U);
			}
		}
	}

	void ChSuperLUEngine::SetMatrix(ChSparseMatrix& Z)
	{

		SetMatrix(Z.GetNumRows(), Z.GetCSR_ValueArray(), Z.GetCSR_TrailingIndexArray(), Z.GetCSR_LeadingIndexArray());
		//TODO: change here the matrix type (symmetric, etc) according to what Z says
	}

	void ChSuperLUEngine::SetMatrix(int pb_size, double* values, int* rowIndex, int* colIndex)
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


		get_perm_c(1, &m_mat_Super, perm_c.data());

		superlumt_options.perm_c = perm_c.data();
		superlumt_options.perm_r = perm_r.data();

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

	void ChSuperLUEngine::SetSolutionVector(ChMatrix<>& x)
	{
		assert(x.GetRows() == m_n);

		SetSolutionVector(x.GetAddress());
	}

	void ChSuperLUEngine::SetSolutionVector(double* x)
	{
		m_sol_Super.nrow = m_n;
		m_sol_Super.ncol = 1;
		auto m_sol_Super_store = static_cast<DNformat*>(m_sol_Super.Store);
		//m_sol_Super_store->lda = m_n; // done during the SuperLU call
		m_sol_Super_store->nzval = x;

		ferr.resize(m_nrhs);
	}

	void ChSuperLUEngine::SetRhsVector(ChMatrix<>& b)
	{
		assert(b.GetRows() == m_n);

		// that should deal multiple rhs
		if (b.GetColumns()>1)
			b.MatrTranspose(); // ChMatrix is row-major

		SetRhsVector(b.GetAddress(), b.GetColumns());
	}

	void ChSuperLUEngine::SetRhsVector(double* b, int nrhs)
	{
		m_nrhs = nrhs;

		m_rhs_Super.nrow = m_n;
		m_rhs_Super.ncol = nrhs;
		auto m_rhs_Super_store = static_cast<DNformat*>(m_rhs_Super.Store);
		//m_rhs_Super_store->lda = m_n; // done during the SuperLU call
		m_rhs_Super_store->nzval = b;

		berr.resize(m_nrhs);
	}

	void ChSuperLUEngine::SetProblem(ChSparseMatrix& Z, ChMatrix<>& b, ChMatrix<>& x)
	{
		SetMatrix(Z);
		SetRhsVector(b);
		SetSolutionVector(x);
	}

	int ChSuperLUEngine::SuperLUCall(int phase, int verbose)
	{
		auto m_rhs_Super_ncol_bkp = m_rhs_Super.ncol;
		auto m_sol_Super_ncol_bkp = m_sol_Super.ncol;

		// Set the number of right-hand side
		// if put to 0 then factorize without solve
		m_rhs_Super.ncol = (phase == ANALYSIS_NUMFACTORIZATION) ? 0 : m_nrhs;
		m_sol_Super.ncol = (phase == ANALYSIS_NUMFACTORIZATION) ? 0 : m_nrhs;
		superlumt_options.fact = (phase == SOLVE) ? FACTORED : DOFACT; /* Indicate the factored form of m_mat_Super is supplied. */
		//TODO: superlumt_options.refact


		// call to SuperLU
		pdgssvx(nprocs, &superlumt_options, &m_mat_Super, perm_c.data(), perm_r.data(), &equed, R.data(), C.data(),
			&L, &U, &m_rhs_Super, &m_sol_Super, &rpg, &rcond, ferr.data(), berr.data(),
			&superlu_memusage, &info);

		// restore previous values of columns of rhs and sol
		// otherwise if someone set rhs and sol before factorization they will be corrupted
		m_rhs_Super.ncol = m_rhs_Super_ncol_bkp;
		m_sol_Super.ncol = m_sol_Super_ncol_bkp;

		if (verbose>1 && (phase == ANALYSIS_NUMFACTORIZATION || phase == COMPLETE))
		{
			if (info == 0 || info == m_n + 1)
			{
				printf("Recip. pivot growth = %e\n", rpg);
				printf("Recip. condition number = %e\n", rcond);
				printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
				for (i = 0; i < m_nrhs; ++i) {
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
			{
				printf("** Estimated memory: %d bytes\n", info - m_n);
			}

		}

		if (verbose>1 && (phase == SOLVE || phase == COMPLETE))
		{
			if (info == 0 || info == m_n + 1)
			{
				printf("Triangular solve: dgssvx() returns info %d\n", info);
			}
			else if (info > 0 && lwork == -1)
			{
				printf("** Estimated memory: " IFMT " bytes\n", info - m_n);
			}
		}
		

		return info;
	}

	void ChSuperLUEngine::ResetSolver()
	{
		superlumt_options.nprocs = nprocs;
		superlumt_options.fact = fact;
		superlumt_options.trans = trans;
		superlumt_options.refact = refact;
		superlumt_options.panel_size = panel_size;
		superlumt_options.relax = relax;
		superlumt_options.usepr = usepr;
		superlumt_options.drop_tol = 0.0;
		superlumt_options.diag_pivot_thresh = 1.0;
		superlumt_options.SymmetricMode = NO;
		superlumt_options.PrintStat = NO;
		superlumt_options.perm_c = perm_c.data();
		superlumt_options.perm_r = perm_r.data();
		superlumt_options.work = work;
		superlumt_options.lwork = lwork;
	}

	void ChSuperLUEngine::GetResidual(ChMatrix<>& res) const
	{
		assert(res.GetRows() >= m_n);
		GetResidual(res.GetAddress());
	}

	void ChSuperLUEngine::GetResidual(double* res) const
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

	double ChSuperLUEngine::GetResidualNorm() const
	{
		std::vector<double> res(m_n);
		GetResidual(res.data());
		double norm = 0;
		for (int i = 0; i < m_n; i++) {
			norm += res[i] * res[i];
		}

		// also BLAS function can be used
		return std::sqrt(norm);
	}
}  // end namespace chrono
