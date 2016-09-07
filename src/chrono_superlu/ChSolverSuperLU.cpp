// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by values_vect BSD-style license that can be found
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

#include "ChSolverSuperLU.h"

namespace chrono {


	template <typename Matrix>
	bool ChSolverSuperLU<Matrix>::Setup(ChSystemDescriptor& sysd)
	{
		m_timer_setup_assembly.start();

		// Set the lock on the matrix sparsity pattern (if enabled).
		m_mat.SetSparsityPatternLock(m_lock);

		// Calculate problem size at first call.
		if (m_solver_call == 0)
		{
			m_n = sysd.CountActiveVariables() + sysd.CountActiveConstraints();
		}

		// If an NNZ value for the underlying matrix was specified, perform an initial resizing, *before*
		// values_vect call to ChSystemDescriptor::ConvertToMatrixForm(), to allow for possible size optimizations.
		// Otherwise, do this only at the first call, using the default sparsity fill-in.
		if (m_nnz != 0)
		{
			m_mat.Reset(m_n, m_n, m_nnz);
		}
		else if (m_solver_call == 0)
		{
			m_mat.Reset(m_n, m_n, static_cast<int>(m_n * (m_n * SPM_DEF_FULLNESS)));
		}

		// Assemble the matrix.
		sysd.ConvertToMatrixForm(&m_mat, nullptr);
		m_n = m_mat.GetNumRows();


		// Allow the matrix to be compressed.
		bool change = m_mat.Compress();

		m_timer_setup_assembly.stop();

		/* ONLY PERFORM THE LU DECOMPOSITION */
		m_engine.SetMatrix(m_mat);

		m_timer_setup_superlu.start();
		m_engine.SuperLUCall(12);
		m_timer_setup_superlu.stop();


		return true;
	}

	template <typename Matrix>
	double ChSolverSuperLU<Matrix>::Solve(ChSystemDescriptor& sysd)
	{
		// Assemble the problem right-hand side vector.
		m_timer_solve_assembly.start();
		sysd.ConvertToMatrixForm(nullptr, &m_rhs);
		m_sol.Resize(m_rhs.GetRows(), 1);
		m_engine.SetRhsVector(m_rhs);
		m_engine.SetSolutionVector(m_sol);
		m_timer_solve_assembly.stop();


		/* ------------------------------------------------------------
		NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF m_mat_Super.
		------------------------------------------------------------*/

		m_timer_solve_superlu.start();
		m_engine.SuperLUCall(33);
		m_timer_solve_superlu.stop();
		m_solver_call++;


		// Scatter solution vector to the system descriptor.
		m_timer_solve_assembly.start();
		sysd.FromVectorToUnknowns(m_sol);
		m_timer_solve_assembly.stop();

		return 0.0f;
	}

}  // end namespace chrono
