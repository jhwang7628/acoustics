#ifndef LINEARALGEBRA_EIGS_H
#   define LINEARALGEBRA_EIGS_H

#include "PardisoMatrix.hpp"

/*
 * compute the \nev standard eigenvalues for the given symmetric 
 * matrix. The result eigenvalues are put into \eval, and eigenvectors
 * are put into \evec
 */
void std_symm_eigs(const PardisoMatrix<double>& mat, 
        int nev, double* eval, double* evec, char* op = "LM");

#endif
