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
// Interfacing to the SuperLU_MT solver.
// =============================================================================

#ifndef CHSUPERLUMTENGINE_H
#define CHSUPERLUMTENGINE_H

#include "chrono_superlumt/ChApiSuperLUMT.h"
#include "chrono/core/ChSparseMatrix.h"
#include "chrono/core/ChMatrixDynamic.h"
#include <slu_mt_ddefs.h>


namespace chrono {

/// @addtogroup superlumt_module
/// @{

/// Interface class to SuperLU_MT solver.
/// This class wraps the C interface of the solver in order to fit Chrono data structures.
/// This class can still be called by the end-user in order to solve linear systems.
/// See demo_SUPERLU_Engine for the related demo.
class ChApiSuperLUMT ChSuperLUMTEngine {
  public:
    ChSuperLUMTEngine();
    ~ChSuperLUMTEngine();

    /// Set problem dimension.
    void SetProblemSize(int pb_size) { m_n = pb_size; }

    /// Set the problem matrix.
    /// This will also update the problem dimension as well as the matrix symmetry type.
    void SetMatrix(ChSparseMatrix& Z);

    /// Set directly the CSR matrix arrays.
    /// Note that it is implied that the matrix symmetry type is GENERAL.
    void SetMatrix(int pb_size, double* values, int* rowIndex, int* colIndex);

    /// Set the solution vector.
    /// Note that it is the caller's responsibility to provide an array of appropriate size.
    void SetSolutionVector(ChMatrix<>& x);
    void SetSolutionVector(double* x);

    /// Set the right-hand side vector.
    /// Note that it is the caller's responsibility to ensure that the size is appropriate.
    void SetRhsVector(ChMatrix<>& b);
    void SetRhsVector(double* b, int nrhs = 1);

    /// Set the matrix, as well as the right-hand side and solution arrays.
    void SetProblem(ChSparseMatrix& Z, ChMatrix<>& b, ChMatrix<>& x);

    /// Solver routine.
    int SuperLUMTCall(int phase, int verbose = 0);

    /// Reinitializes the solver to default values.
    void ResetSolver();

	// Auxiliary functions
	/// Returns the Options vector
	superlumt_options_t& GetOptions() { return superlumt_options; }


    // Output functions
    /// Calculate and return the problem residual res=b-Ax.
    /// Note that it is the caller's responsibility to provide an array of appropriate size.
    void GetResidual(ChMatrix<>& res) const;
    void GetResidual(double* res) const;

    /// Calculate and return the L2-norm of the problem residual, ||b-Ax||.
    double GetResidualNorm() const;

  private:

	// Problem properties
	int m_n = 0;     ///< (square) matrix size
	int m_nrhs = 1;  ///< number of rhs vectors

	/* SuperLU_MT data */
	int& ldx = m_n;

	SuperMatrix    m_mat_Super, m_rhs_Super, m_sol_Super;
	std::vector<int> perm_c; /* column permutation vector */
	std::vector<int> perm_r; /* row permutations from partial pivoting */
	std::vector<int> etree;

	int            lwork = 0; // allocate space internally by system malloc (don't use 'work' variable)
	int            i = 0, nnz = 0;

	std::vector<double> R = {0.0}, C = { 0.0 };
	std::vector<double> ferr = {0.0}, berr={0.0};
	double         rpg, rcond;

	// internally used and never directly modified by user
	SuperMatrix    L, U;
	void*          work = nullptr; // don't used to allocate space
	int            info = 0;

	// SuperLU_MT datas (different from serial SuperLU)
	superlu_memusage_t    superlu_memusage;
	superlumt_options_t superlumt_options;
	int         nprocs = 4;
	fact_t      fact = EQUILIBRATE;
	trans_t     trans = NOTRANS;
	yes_no_t    refact = NO, usepr = NO;
	equed_t     equed = NOEQUIL;
	int relax, panel_size;
	std::vector<int> colcnt_h;
	std::vector<int> part_super_h;



    // SuperLU_MT solver settings
};



	/// @} superlumt_module

}  // end of namespace chrono

#endif