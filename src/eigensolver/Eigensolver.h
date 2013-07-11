//////////////////////////////////////////////////////////////////////
// Eigensolver.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef ARNOLDI_EIGENSOLVER_H
#define ARNOLDI_EIGENSOLVER_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/SPARSE_MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <superlu-interface/SuperLU_Interface.h>

#include <TYPES.h>

#include "areigs.h"

//////////////////////////////////////////////////////////////////////
// Eigensolver class
//
// Solve generalized eigenvalue problems via a shift-and-invert
// method.  Matrix inversion is accomplished via the SuperLU library
//////////////////////////////////////////////////////////////////////
class Eigensolver : public ArpackDS {
    public:
        // Set up the solver for the generalized eigenvalue problem
        //  K * u = l * M * u
        Eigensolver( const SPARSE_MATRIX &K, const SPARSE_MATRIX &M,
                     double shift );

        // Destructor
        virtual ~Eigensolver();

    protected:
        // Applies the matrix (K - sigma_ * M)^{-1} * M to the input vector
        // and places the result in opx.  Also puts M * x in Mx
        virtual void times_OP1(double* x, double* opx);
        virtual void times_OP2(double* x, double* opx, double* Mx);
        virtual void times_M  (double* x, double* Mx);

    private:
        const SPARSE_MATRIX                 &_K;
        SPARSE_MATRIX::SparseColumnMatrix    _M;

        // Stores K - sigma_ * M
        SPARSE_MATRIX                        _inverseMatrix;

        SuperLU_Solver                      *_linearSolver;

};

#endif
