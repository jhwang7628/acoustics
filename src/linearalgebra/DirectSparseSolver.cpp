#include "DirectSparseSolver.h"
#include <stdlib.h>
#include <stdio.h>

void DirectSparseSolver::load_matrix(const PardisoMatrix<double>* A)
{
    _INTEGER_t error;

    if ( m_allocated )
    {
        error = dss_delete(m_handle, m_opt);
        if ( error != MKL_DSS_SUCCESS )
        {
            fprintf(stderr, "[dss_delete] Solver returned error code %d\n", error);
            exit(1);
        }
        m_allocated = false;
    }

    error = dss_create(m_handle, m_opt);
    if ( error != MKL_DSS_SUCCESS )
    {
        fprintf(stderr, "[dss_create] Solver returned error code %d\n", error);
        exit(1);
    }
    m_allocated = true;
    //// define the structure
    int sym = A->is_symmetric() ? MKL_DSS_SYMMETRIC : MKL_DSS_NON_SYMMETRIC;

    int nr = A->num_rows();
    int nc = A->num_cols();
    int nnz = A->num_nonzeros();
    error = dss_define_structure(m_handle, sym, 
            A->row_indices(), nr, nc, 
            A->col_indices(), nnz);
    if ( error != MKL_DSS_SUCCESS )
    {
        fprintf(stderr, "[dss_define_structure] Solver returned error code %d\n", error);
        exit(1);
    }

    //// reorder the matrix
    error = dss_reorder(m_handle, m_opt, 0);
    if ( error != MKL_DSS_SUCCESS )
    {
        fprintf(stderr, "[dss_reorder] Solver returned error code %d\n", error);
        exit(1);
    }

    int type = A->is_positive_definite() ? MKL_DSS_POSITIVE_DEFINITE : MKL_DSS_INDEFINITE;
    //// factorize the matrix
    error = dss_factor_real(m_handle, type, A->data());
    if ( error != MKL_DSS_SUCCESS )
    {
        fprintf(stderr, "[dss_factor_real] Solver returned error code %d\n", error);
        exit(1);
    }
}

void DirectSparseSolver::solve(const double* b, double* x)
{
    const int NRHS = 1;
    _INTEGER_t error = dss_solve_real(m_handle, m_opt, b, NRHS, x);
    if ( error != MKL_DSS_SUCCESS )
    {
        fprintf(stderr, "[dss_solve_real] Solver returned error code %d\n", error);
        exit(1);
    }
}

