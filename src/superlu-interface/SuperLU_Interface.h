//////////////////////////////////////////////////////////////////////
// SuperLU_Interface.h: Interface for the SuperLU_Interface class
//
//////////////////////////////////////////////////////////////////////

#ifndef SUPERLU_INTERFACE_H
#define SUPERLU_INTERFACE_H

#include <TYPES.h>

#include <linearalgebra/SPARSE_MATRIX.h>
#include <linearalgebra/VECTOR.h>

//////////////////////////////////////////////////////////////////////
// SuperLU includes
//////////////////////////////////////////////////////////////////////
#include <slu_ddefs.h>

//////////////////////////////////////////////////////////////////////
// SuperLU_Solver class
//
// Simple C++ interface for the SuperLU library
//////////////////////////////////////////////////////////////////////
class SuperLU_Solver {
    public:
        // Initialize the solver and build factors for the given
        // sparse matrix
        SuperLU_Solver( const SPARSE_MATRIX &A );
        SuperLU_Solver( const SPARSE_MATRIX::SparseColumnMatrix &A );

        // Destructor
        virtual ~SuperLU_Solver();

        // Solves the system with the given right hand side
        void solve( const REAL *rhs, REAL *solution );
        void solve( const VECTOR &rhs, VECTOR &solution );

    protected:

    private:
        // Initialize from a matrix in compressed column format
        void init( const SPARSE_MATRIX::SparseColumnMatrix &A );

    private:
        // System matrix in SuperLU format
        SuperMatrix          _A;

        // We also have to explicitly store the data fields for this matrix
        REAL                *_nzval;
        int                 *_rowindices;
        int                 *_colptr;

        // Factor matrices in SuperLU format
        SuperMatrix          _L, _U;

        // Other SuperLU data fields
        superlu_options_t    _options;
        SuperLUStat_t        _stat;

        // Permutation vectors
        int                 *_columnPermutation;
        int                 *_rowPermutation;

};

#endif
