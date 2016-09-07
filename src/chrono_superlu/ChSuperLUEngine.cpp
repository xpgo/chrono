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
#include <cmath>

#include "chrono_superlu/ChSuperLUEngine.h"

namespace chrono {
	enum phase_t {
		COMPLETE = 13,
		ANALYSIS_NUMFACTORIZATION = 12,
		SOLVE = 33,
	};

	ChSuperLUEngine::ChSuperLUEngine(int pb_size): m_n(pb_size)
	{
		/* Defaults */
		lwork = 0;

		set_default_options(&options);

		options.Equil = YES;
		options.DiagPivotThresh = 1.0;
		options.Trans = NOTRANS;
	}

	ChSuperLUEngine::~ChSuperLUEngine()
	{
		StatFree(&stat);
		//Destroy_CompCol_Matrix(&m_mat_Super); this would destroy also the CSR arrays
		Destroy_SuperMatrix_Store(&m_mat_Super);
		Destroy_SuperMatrix_Store(&m_rhs_Super);
		Destroy_SuperMatrix_Store(&m_sol_Super);
		if (lwork == 0) {
			Destroy_SuperNode_Matrix(&L);
			Destroy_CompCol_Matrix(&U);
		}
	}

	void ChSuperLUEngine::SetMatrix(ChSparseMatrix& Z)
	{

		// this function is identical to SetMatrix(int pb_size, double* values, int* rowIndex, int* colIndex)
		// but it is here because the type of the matrix can be deduced in this case from Z
		m_n = Z.GetNumRows();
		dCreate_CompCol_Matrix(&m_mat_Super, m_n, m_n, Z.GetNNZ(), Z.GetCSR_ValueArray(), Z.GetCSR_TrailingIndexArray(), Z.GetCSR_LeadingIndexArray(), SLU_NC, SLU_D, SLU_GE);

		etree.resize(m_n);
		perm_r.resize(m_n);
		perm_c.resize(m_n);

		R.resize(m_mat_Super.nrow);
		C.resize(m_mat_Super.ncol);
	}

	void ChSuperLUEngine::SetMatrix(int pb_size, double* values, int* rowIndex, int* colIndex)
	{
		m_n = pb_size;
		m_a = values;
		m_ia = rowIndex;
		m_ja = colIndex;
		dCreate_CompCol_Matrix(&m_mat_Super, m_n, m_n, m_ja[m_n], m_a, m_ia, m_ja, SLU_NC, SLU_D, SLU_GE);
		auto Astore = static_cast<double*>(m_mat_Super.Store);

		etree.resize(m_n);
		perm_r.resize(m_n);
		perm_c.resize(m_n);

		R.resize(m_mat_Super.nrow);
		C.resize(m_mat_Super.ncol);
	}

	void ChSuperLUEngine::SetSolutionVector(ChMatrix<>& x)
	{
		SetSolutionVector(x.GetAddress());
	}

	void ChSuperLUEngine::SetSolutionVector(double* x)
	{
		m_x = x;
		dCreate_Dense_Matrix(&m_sol_Super, m_n, m_nrhs, m_x, m_n, SLU_DN, SLU_D, SLU_GE);

		ferr.resize(m_nrhs);
	}

	void ChSuperLUEngine::SetRhsVector(ChMatrix<>& b)
	{
		// that should deal multiple rhs
		if (b.GetColumns()>1)
			b.MatrTranspose(); // ChMatrix is row-major

		SetRhsVector(b.GetAddress(), b.GetColumns());
	}

	void ChSuperLUEngine::SetRhsVector(double* b, int nrhs)
	{
		m_b = b;
		m_nrhs = nrhs;
		dCreate_Dense_Matrix(&m_rhs_Super, m_n, m_nrhs, m_b, m_n, SLU_DN, SLU_D, SLU_GE);

		berr.resize(m_nrhs);
	}

	void ChSuperLUEngine::SetProblem(ChSparseMatrix& Z, ChMatrix<>& b, ChMatrix<>& x)
	{
		SetMatrix(Z);
		SetRhsVector(b);
		SetSolutionVector(x);
	}

	int ChSuperLUEngine::SuperLUCall(int phase, bool verbose)
	{

		/* Initialize the statistics variables. */
		StatInit(&stat);

		// Set the number of right-hand side
		// if put to 0 then factorize without solve
		m_rhs_Super.ncol = (phase == ANALYSIS_NUMFACTORIZATION) ? 0 : m_nrhs;
		options.Fact = (phase == SOLVE) ? FACTORED : DOFACT; /* Indicate the factored form of m_mat_Super is supplied. */

		
		// call to SuperLU
		dgssvx(&options, &m_mat_Super, perm_c.data(), perm_r.data(), etree.data(), equed, R.data(), C.data(),
			&L, &U, work, lwork, &m_rhs_Super, &m_sol_Super, &rpg, &rcond, ferr.data(), berr.data(),
			&Glu, &mem_usage, &stat, &info);


		if (phase != ANALYSIS_NUMFACTORIZATION)
		{
			/* This is how you could access the solution matrix. */
			double* sol = static_cast<double*>(static_cast<DNformat*>(m_sol_Super.Store)->nzval);

			for (auto row_sel = 0; row_sel<m_n; ++row_sel)
			{
				m_x[row_sel]= sol[row_sel];
			}
		}

		if (phase == ANALYSIS_NUMFACTORIZATION || phase == COMPLETE)
		{
			if (verbose && (info == 0 || info == m_n + 1))
			{
				printf("LU factorization: dgssvx() returns info %d\n", info);

				if (options.PivotGrowth) printf("Recip. pivot growth = %e\n", rpg);
				if (options.ConditionNumber)
					printf("Recip. condition number = %e\n", rcond);
				auto Lnnz = static_cast<SCformat*>(L.Store)->nnz;
				auto Unnz = static_cast<NCformat*>(U.Store)->nnz;
				printf("No of nonzeros in factor L = %d\n", Lnnz);
				printf("No of nonzeros in factor U = %d\n", Unnz);
				printf("No of nonzeros in L+U = %d\n", Lnnz + Unnz - m_n);
				printf("FILL ratio = %.1f\n", static_cast<float>(Lnnz + Unnz - m_n) / nnz);

				printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
					mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6);
			}
			else if (verbose && (info > 0 && lwork == -1))
			{
				printf("** Estimated memory: %d bytes\n", info - m_n);
			}

			if (options.PrintStat && phase == ANALYSIS_NUMFACTORIZATION) StatPrint(&stat);
		}

		if (phase == SOLVE || phase == COMPLETE)
		{
			if (verbose && (info == 0 || info == m_n + 1))
			{
				printf("Triangular solve: dgssvx() returns info %d\n", info);

				if (options.IterRefine)
				{
					printf("Iterative Refinement:\n");
					printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
					for (i = 0; i < m_nrhs; ++i)
						printf("%8d%8d%16e%16e\n", i + 1, stat.RefineSteps, ferr[i], berr[i]);
				}
			}
			else if (verbose && (info > 0 && lwork == -1))
			{
				printf("** Estimated memory: %d bytes\n", info - m_n);
			}

			if (options.PrintStat) StatPrint(&stat);
		}
		

		StatFree(&stat);

		return info;
	}

}  // end namespace chrono
