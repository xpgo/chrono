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
// Authors: Radu Serban
// =============================================================================

#ifndef CH_MAP_MATRIX_H
#define CH_MAP_MATRIX_H

#include <vector>
#include <unordered_map>

#include "chrono/core/ChSparseMatrix.h"
#include "chrono/core/ChMatrixDynamic.h"

namespace chrono {

// -----------------------
// SPARSE MAP MATRIX CLASS
// -----------------------

/// This class defines a sparse matrix, implemented using unordered_maps for each row.
class ChApi ChMapMatrix : public ChSparseMatrix {
  public:
    /// Create a sparse matrix with given dimensions.
    ChMapMatrix(int nrows = 1, int ncols = 1);

    /// Create a sparse matrix from a given dense matrix.
    ChMapMatrix(const ChMatrix<>& mat);

    /// Copy constructor.
    ChMapMatrix(const ChMapMatrix& other);

    /// Destructor.
    ~ChMapMatrix();

    /// Resize this matrix.
    virtual bool Resize(int nrows, int ncols, int nonzeros = 0) override;

    /// Reset to null matrix and (if needed) changes the size.
    virtual void Reset(int nrows, int ncols, int nonzeros = 0) override;

    /// Set/update the specified matrix element.
    virtual void SetElement(int row, int col, double elem, bool overwrite = true) override;

    /// Get the element at the specified location.
    virtual double GetElement(int row, int col) const override;

    /// Return the row index array in the CSR representation of this matrix.
    virtual int* GetCSR_LeadingIndexArray() const override;

    /// Return the column index array in the CSR representation of this matrix.
    virtual int* GetCSR_TrailingIndexArray() const override;

    /// Return the array of matrix values in the CSR representation of this matrix.
    virtual double* GetCSR_ValueArray() const override;

    /// Convert to dense matrix.
    void ConvertToDense(ChMatrixDynamic<double>& mat);

    /// Convert to CSR arrays.
	template <class index_vector_t, class value_vector_t>
	void ConvertToCSR(index_vector_t& ia, index_vector_t& ja, value_vector_t& a) const
	{
		ia.resize(m_num_rows + 1);
		ja.resize(m_nnz);
		a.resize(m_nnz);

		// Loop over all rows, accumulating the number of non-zero elements.
		ia[0] = 0;
		int nnz = 0;
		for (int ir = 0; ir < m_num_rows; ir++) {
			// Update row index array.
			ia[ir + 1] = ia[ir] + m_rows[ir].m_nnz;

			// Copy the keys for the row data into a vector and sort them.
			std::vector<int> col_idx;
			col_idx.reserve(m_rows[ir].m_nnz);
			for (auto& it : m_rows[ir].m_data) {
				col_idx.push_back(it.first);
			}
			std::sort(col_idx.begin(), col_idx.end());

			// Extract the non-zero elements in the row, in ascending order of their keys.
			// Note that we need not test the return iterator from find here (key guaranteed to exist).
			// Update the column index and value arrays.
			for (auto ic : col_idx) {
				ja[nnz] = ic;
				a[nnz] = m_rows[ir].m_data.find(ic)->second;
				nnz++;
			}
		}

		m_CSR_current = true;
	}

    /// Method to allow serializing transient data into in ASCII stream (e.g., a file) as a
    /// Matlab sparse matrix format; each row in file has three elements: {row, column, value}.
    /// Note: the row and column indexes start from 1.
    void StreamOUTsparseMatlabFormat(ChStreamOutAscii& mstream);

    /// Write first few rows and columns to the console.
    /// Method to allow serializing transient data into in ASCII format.
    void StreamOUT(ChStreamOutAscii& mstream);

  private:
    struct MatrixRow {
        MatrixRow() : m_nnz(0) {}
		int m_nnz;                               ///< number of non-zero elements in row
        std::unordered_map<int, double> m_data;  ///< column - value pairs in row
    };

    std::vector<MatrixRow> m_rows;  ///< vector of data structures for each row

    // Storage for CSR arrays. These are allocated only if needed.
    mutable std::vector<int> m_ia;    ///< CSR row index array
    mutable std::vector<int> m_ja;    ///< CSR column index array
    mutable std::vector<double> m_a;  ///< CSR values array

    mutable bool m_CSR_current;  ///< flag indicating whether or not the CSR arrays are up-to-date
};

}  // end namespace chrono

#endif
