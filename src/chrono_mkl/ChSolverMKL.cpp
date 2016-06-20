#include "chrono_mkl/ChSolverMKL.h"

namespace chrono {
// Register into the object factory, to enable run-time
// dynamic creation and persistence
ChClassRegister<ChSolverMKL> a_registration_ChSolverMKL;

/** \brief It calls Intel MKL Pardiso Sparse Direct Solver.
*
*	On the first call the routine resets the matrix and vectors involved to the size of the current problem.
*	Every time the system stored in \c sysd is converted in a CSR3 standard matrix that will be later used by Pardiso.
*	If the sparsity pattern lock is turned on then the matrix will try, as long as possible, to preserve not only the
*   arrays dimensions, but it also keeps column and row indexes through different calls.
*/

double ChSolverMKL::Solve(ChSystemDescriptor& sysd) {
    if (!manual_factorization)
        Factorize(sysd);

    timer_buildmat.start();
    sysd.ConvertToMatrixForm(nullptr, &rhs);
    timer_buildmat.stop();

    mkl_engine.SetProblem(matCSR3, rhs, sol);
    timer_solve.start();
    int pardiso_message_phase33 = mkl_engine.PardisoCall(33, 0);
    timer_solve.stop();

    solve_phase_calls++;

    if (pardiso_message_phase33) {
        GetLog() << "Pardiso solve+refine error code = " << pardiso_message_phase33 << "\n";
        GetLog() << "Matrix verification code = " << matCSR3.VerifyMatrix() << "\n";
        GetLog() << "Matrix MKL verification code = " << matCSR3.VerifyMatrixByMKL() << "\n";
    }

    if (verbose) {
        // Get residual;
        mkl_engine.GetResidual(res);
        double res_norm = mkl_engine.GetResidualNorm(res);

        GetLog() << "Pardiso solve phase call " << solve_phase_calls << "  |residual| = " << res_norm << "\n";
    }

    // Replicate the changes to vvariables and vconstraint into SystemDescriptor
    sysd.FromVectorToUnknowns(sol);

    return 0.0f;
}

double ChSolverMKL::Factorize(ChSystemDescriptor& sysd) {
    // Initial resizing;
	timer_buildmat.start();
    if (factorize_phase_calls == 0)
    {
        // not mandatory, but it speeds up the first build of the matrix, guessing its sparsity; needs to stay BEFORE ConvertToMatrixForm()
        n = sysd.CountActiveVariables() + sysd.CountActiveConstraints();
        nnz ? matCSR3.Reset(n, n, nnz) : matCSR3.Reset(n, n, static_cast<int>(n*(n*SPM_DEF_FULLNESS)));
        sol.Resize(n, 1); // ConvertToMatrixForm() takes care of eventually resizing matCSR3 and rhs, but not sol; this can be done also AFTER CTMF()
        res.Resize(n, 1);
    }
    else
        if (nnz)
            matCSR3.Reset(n, n, nnz);

    // Build matrix and rhs;
    // in case the matrix changes size (rows) then RowIndexLockBroken is turned on
    // in case the matrix has to insert a new element in the matrix then ColIndexLockBroken is turned on
    // in case the matrix inserts LESS elements than the previous cycle:
    //     - if sparsity_pattern_lock is ON the matrix will remain with unuseful filling zeros; they have to be
    //     eventually removed with a Purge();
    //     - if sparsity_pattern_lock is ON the matrix will have some uninitialized element flagged with -1 in the
    //     colIndex;
    //          those have to be removed by Compress().
    // manual_factorization ? sysd.ConvertToMatrixForm(&matCSR3, nullptr) : sysd.ConvertToMatrixForm(&matCSR3, &rhs);
    sysd.ConvertToMatrixForm(&matCSR3, nullptr);

    // Set up the locks;
    // after the first build the sparsity pattern is eventually locked; needs to stay AFTER ConvertToMatrixForm()
    if (factorize_phase_calls == 0) {
        matCSR3.SetRowIndexLock(sparsity_pattern_lock);
        matCSR3.SetColIndexLock(sparsity_pattern_lock);
    }

    // remember: the matrix is constructed with broken locks; so for solver_call==0 they are broken!
    if (!sparsity_pattern_lock || matCSR3.IsRowIndexLockBroken() || matCSR3.IsColIndexLockBroken()) {
        // breaking the row_index_block means that the size has changed so every dimension has to be reset
        if (matCSR3.IsRowIndexLockBroken()) {
            n = matCSR3.GetRows();
            sol.Resize(n, 1);
            res.Resize(n, 1);
        }

        // if sparsity is not locked OR the sparsity_lock is broken (like in the first cycle!); the matrix must be
        // recompressed
        matCSR3.Compress();

        // the permutation vector is based on the sparsity of the matrix;
        // if ColIndexLockBroken/RowIndexLockBroken are on then the permutation must be updated
        if (use_perm)
            mkl_engine.UsePermutationVector(true);
    }

    // the sparsity of rhs must be updated at every cycle (am I wrong?)
    if (use_rhs_sparsity && !use_perm)
        mkl_engine.UsePartialSolution(2);
	
	timer_buildmat.stop();

    // Solve with Pardiso Sparse Direct Solver
    // the problem size must be updated also in the Engine: this is done by SetProblem() itself.
    mkl_engine.SetProblem(matCSR3, rhs, sol);
    timer_factorize.start();
    int pardiso_message_phase12 = mkl_engine.PardisoCall(12, 0);
    timer_factorize.stop();

    factorize_phase_calls++;

    if (pardiso_message_phase12) {
        GetLog() << "Pardiso analyze+reorder+factorize error code = " << pardiso_message_phase12 << "\n";
        GetLog() << "Matrix verification code = " << matCSR3.VerifyMatrix() << "\n";
        GetLog() << "Matrix MKL verification code = " << matCSR3.VerifyMatrixByMKL() << "\n";
    }

    return pardiso_message_phase12;
}

}  // end namespace chrono
