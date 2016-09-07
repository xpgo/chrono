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

#ifndef CHSOLVERSUPERLU_H
#define CHSOLVERSUPERLU_H

#include "chrono/solver/ChSolver.h"
#include "chrono/solver/ChSystemDescriptor.h"
#include "chrono/core/ChSparseMatrix.h"
#include "chrono/core/ChMatrixDynamic.h"
#include "chrono/core/ChTimer.h"

#include "chrono/core/ChCSR3Matrix.h"
#include "chrono_superlu/ChSuperLUEngine.h"

namespace chrono {

/// @addtogroup superlu_module
/// @{

/// Class that wraps the SuperLU solver.
/// It can solve linear systems, but not VI and complementarity problems.
template <typename Matrix = ChCSR3Matrix>
class ChSolverSuperLU : public ChSolver {

public:
	ChSolverSuperLU(){}

	~ChSolverSuperLU(){}

	Matrix& GetMatrix() { return m_mat; }

	/// Enable/disable locking the sparsity pattern (default: false).
	/// If on_off is set to true, then the sparsity pattern of the problem matrix is assumed
	/// to be unchanged from call to call.
	void SetSparsityPatternLock(bool val) { m_lock = val; }

	/// Set the number of non-zero entries in the problem matrix.
	void SetMatrixNNZ(int nnz) { m_nnz = nnz; }

	/// Reset timers for internal phases in Solve and Setup.
	void ResetTimers() {
		m_timer_setup_assembly.reset();
		m_timer_setup_superlu.reset();
		m_timer_solve_assembly.reset();
		m_timer_solve_superlu.reset();
	}

	/// Get cumulative time for assembly operations in Solve phase.
	double GetTimeSolveAssembly() const { return m_timer_solve_assembly(); }
	/// Get cumulative time for Pardiso calls in Solve phase.
	double GetTimeSolveSuperLU() const { return m_timer_solve_superlu(); }
	/// Get cumulative time for assembly operations in Setup phase.
	double GetTimeSetupAssembly() const { return m_timer_setup_assembly(); }
	/// Get cumulative time for Pardiso calls in Setup phase.
	double GetTimeSetupSuperLU() const { return m_timer_setup_superlu(); }

	/// Indicate whether or not the Solve() phase requires an up-to-date problem matrix.
	/// As typical of direct solvers, the Pardiso solver only requires the matrix for its Setup() phase.
	virtual bool SolveRequiresMatrix() const override { return false; }

	/// Solve using the SuperLU solver.
	/// It uses the matrix factorization obtained at the last call to Setup().
	virtual double Solve(ChSystemDescriptor& sysd) override {
		// Assemble the problem right-hand side vector.
		m_timer_solve_assembly.start();
		sysd.ConvertToMatrixForm(nullptr, &m_rhs);
		m_sol.Resize(m_rhs.GetRows(), 1);
		m_engine.SetRhsVector(m_rhs);
		m_engine.SetSolutionVector(m_sol);
		m_timer_solve_assembly.stop();

		// Solve the system
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

	/// Perform the solver setup operations.
	/// For the SuperLU solver, this means assembling and factorizing the system matrix.
	/// Returns true if successful and false otherwise.
	virtual bool Setup(ChSystemDescriptor& sysd) override {
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

		m_engine.SetMatrix(m_mat);
		// Performs factorization (LU decomposition)
		m_timer_setup_superlu.start();
		m_engine.SuperLUCall(12);
		m_timer_setup_superlu.stop();


		return true;
	}

	/// Method to allow serialization of transient data to archives.
	virtual void ArchiveOUT(ChArchiveOut& marchive) override {
		// version number
		marchive.VersionWrite(1);
		// serialize parent class
		ChSolver::ArchiveOUT(marchive);
		// serialize all member data:
		marchive << CHNVP(m_lock);
	}

	/// Method to allow de serialization of transient data from archives.
	virtual void ArchiveIN(ChArchiveIn& marchive) override {
		// version number
		int version = marchive.VersionRead();
		// deserialize parent class
		ChSolver::ArchiveIN(marchive);
		// stream in all member data:
		marchive >> CHNVP(m_lock);
	}

private:
	Matrix m_mat{1, 1};                 ///< problem matrix
    ChMatrixDynamic<double> m_rhs;      ///< right-hand side vector
	ChMatrixDynamic<double> m_sol;      ///< solution vector
    int m_n = 0;                      ///< problem size
    int m_nnz = 0;                      ///< user-supplied estimate of NNZ
    int m_solver_call = 0;              ///< counter for calls to Solve

    bool m_lock = false;                ///< is the matrix sparsity pattern locked?

	ChSuperLUEngine m_engine;

    ChTimer<> m_timer_setup_assembly;  ///< timer for matrix assembly
    ChTimer<> m_timer_setup_superlu;   ///< timer for factorization
    ChTimer<> m_timer_solve_assembly;  ///< timer for RHS assembly
    ChTimer<> m_timer_solve_superlu;   ///< timer for solution


};

	/// @} superlu_module

}  // end namespace chrono

#endif
