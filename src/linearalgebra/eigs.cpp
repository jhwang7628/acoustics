#include "eigs.h"
#include <arlsmat.h>
#include <arlssym.h>
#include "utils/macros.h"
#include "logging/logging.h"

void std_symm_eigs(const PardisoMatrix<double>& mat,
        int nev, double* eval, double* evec, char* op)
{
    MSG_ASSERT(mat.is_symmetric(), "Given matrix should be symmetric");

    const double* data = mat.data();
    const int*    coli = mat.col_indices();
    const int*    rowi = mat.row_indices();
    const int     nonz = mat.num_nonzeros();
    const int     nrow = mat.num_rows() + 1;

    int * idxcol = new int[nonz];
    int * ptrrow = new int[nrow];

    memcpy(idxcol, coli, nonz*sizeof(int));
    memcpy(ptrrow, rowi, nrow*sizeof(int));

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(idxcol)
#endif
    for(int i = 0;i < nonz;++ i) -- idxcol[i];

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ptrrow)
#endif
    for(int i = 0;i < nrow;++ i) -- ptrrow[i];

    //// create the symmetric matrix
    ARluSymMatrix<double> K(nrow-1, nonz, data, idxcol, ptrrow);
    ARluSymStdEig<double> solver(nev, K, 0, op);
    solver.FindEigenvectors();

    LOGGING_INFO("---------------- Eigen Solver ---------------");
    LOGGING_INFO("max #iter:  %d", solver.GetMaxit());
    LOGGING_INFO("dimension:  %d", solver.GetN());
    LOGGING_INFO("tolenrence: %f", solver.GetTol());

    if ( solver.ConvergedEigenvalues() != nev )
    {
        LOGGING_ERROR("# of converged eigenvalues is not %d", nev);
        delete []idxcol;
        delete []ptrrow;
        exit(1);
    }
    LOGGING_INFO("---------------------------------------------");

    // copy back the results
    memcpy(eval, solver.RawEigenvalues(), nev*sizeof(double));
    memcpy(evec, solver.RawEigenvectors(), nev*(nrow-1)*sizeof(double));

    // release memory
    delete []idxcol;
    delete []ptrrow;
}
