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
// Interfacing to the SuperLU solver.
// =============================================================================

#ifndef CHSUPERLUENGINE_H
#define CHSUPERLUENGINE_H

#include "chrono_superlu/ChApiSuperLU.h"
#include "chrono/core/ChSparseMatrix.h"
#include "chrono/core/ChMatrixDynamic.h"
#include "slu_ddefs.h"


namespace chrono {

/// @addtogroup superlu_module
/// @{

/// Interface class to SuperLU solver.
/// This class wraps the C interface of the solver in order to fit Chrono data structures.
/// This class can still be called by the end-user in order to solve linear systems.
/// See demo_SUPERLU_Engine for the related demo.
class ChApiSuperLU ChSuperLUEngine {
  public:
    ChSuperLUEngine(int pb_size = 0);
    ~ChSuperLUEngine();

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
    int SuperLUCall(int phase, bool verbose = false);

    /// Reinitializes the solver to default values.
    void ResetSolver();

    // Output functions

    /// Calculate and return the problem residual res=b-Ax.
    /// Note that it is the caller's responsibility to provide an array of appropriate size.
    void GetResidual(ChMatrix<>& res) const;
    void GetResidual(double* res) const;

    /// Calculate and return the L2-norm of the problem residual, ||b-Ax||.
    double GetResidualNorm() const;


  private:
	// Data

	// Matrix in CSR3 format.
	// Note that ChSuperLUEngine does not own this data.
	double* m_a = nullptr;    ///< pointer to the CSR array of non-zero elements of the A
	int* m_ia = nullptr;  ///< pointer to the CSR array of row indices
	int* m_ja = nullptr;  ///< pointer to the CSR array of columns indices

	// Right-hand side and solution arrays.
	// Note that ChSuperLUEngine does not own this data.
	double* m_b = nullptr;  ///< rhs vector
	double* m_x = nullptr;  ///< solution vector

				// Problem properties
	int m_n = 0;     ///< (square) matrix size
	int m_nrhs = 1;  ///< number of rhs vectors

	/* SuperLU data */
	int& ldx = m_n;

	char           equed[1];
	SuperMatrix    m_mat_Super, L, U;
	SuperMatrix    m_rhs_Super, m_sol_Super;
	GlobalLU_t	   Glu; /* facilitate multiple factorizations with
						SamePattern_SameRowPerm                  */
	double         *values_vect;
	int            *rowIndex_vect, *colIndex_ptr_vect;
	std::vector<int> perm_c; /* column permutation vector */
	std::vector<int> perm_r; /* row permutations from partial pivoting */
	std::vector<int> etree;
	void*          work = nullptr;
	int            info = 0;
	int            lwork = 0; // allocate space internally by system malloc (don't use 'work' variable)
	int            i =0, nnz=0;
	//std::vector<double> rhsb, rhsx;
	std::vector<double> R, C;
	std::vector<double> ferr, berr;
	double         rpg, rcond;
	mem_usage_t    mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;

    // SuperLU solver settings
};



	/// @} superlu_module

}  // end of namespace chrono

#endif