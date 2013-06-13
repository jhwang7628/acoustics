#include "PardisoSolver.h"
#include <stdio.h>
#include <mkl.h>

#include "utils/nano_timer.h"

/*!
 * solve for Ax = b
 */
void PardisoSolver::solve(const PardisoMatrix<double>& A, 
                          const double* b, double* x)
{
    int phase;

    int idum, nrhs = 1;
    double ddum; /* Double dummy */
    int error;

    m_size = A.num_rows();
    m_type = A.is_symmetric() ? (A.is_positive_definite() ? 2 : -2) : 11;

    double st;
    if ( m_doSymbFact )
    {
        st = GetNanoTimed();
        phase = 11;
        PARDISO(m_pt, &m_maxFact, &m_num, &m_type, &phase, &m_size,
                const_cast<double *>(A.data()), 
                const_cast<int *>(A.row_indices()),
                const_cast<int *>(A.col_indices()),
                &idum, &nrhs, m_parm, &m_msgLvl, 
                &ddum, &ddum, &error);

        if ( error != 0 )
        {
            fprintf(stderr, "ERROR during symbolic factorization: %d\n", error);
            exit(1);
        }
        m_doSymbFact = false;
        printf("Number of nonzeros in factors = %d\n", m_parm[17]);
        printf("Number of factorization MFLOPS = %d\n", m_parm[18]);
        std::cout << "TIME-I: " << GetNanoTimed() - st << std::endl;
    }

    st = GetNanoTimed();
    phase = 23;
    PARDISO(m_pt, &m_maxFact, &m_num, &m_type, &phase, &m_size,
            const_cast<double *>(A.data()), 
            const_cast<int *>(A.row_indices()),
            const_cast<int *>(A.col_indices()),
            &idum, &nrhs, m_parm, &m_msgLvl,
            const_cast<double *>(b), x, &error);
    std::cout << "TIME-II: " << GetNanoTimed() - st << std::endl;
    if ( m_parm[3] % 10 ) // CG/CGS method
        printf("Number of iterations = %d\n", m_parm[19]);
    if ( error != 0 )
    {
        fprintf(stderr, "ERROR during solution: %d\n", error);
        exit(1);
    }


    /*
    // release L and U matrix
    phase = 0;
    PARDISO(m_pt, &m_maxFact, &m_num, &m_type, &phase, &m_size,
            const_cast<double *>(A.data()), 
            const_cast<int *>(A.row_indices()),
            const_cast<int *>(A.col_indices()),
            &idum, &nrhs, m_parm, &m_msgLvl,
            &ddum, &ddum, &error);
    if ( error != 0 )
    {
        fprintf(stderr, "ERROR during release internal memory for L and U: %d\n", error);
        exit(1);
    }
    */
}

void PardisoSolver::release_memory()
{
    int idum, nrhs = 1, phase = -1, error;
    double ddum;
    PARDISO(m_pt, &m_maxFact, &m_num, &m_type, &phase, &m_size,
            &ddum, &idum, &idum,
            &idum, &nrhs, m_parm, &m_msgLvl,
            &ddum, &ddum, &error);
    if ( error != 0 )
    {
        fprintf(stderr, "ERROR during release internal memory: %d\n", error);
        exit(1);
    }
}
