// MATRIX.h: interface for the MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include <TYPES.h>

#include <map>
#include <iostream>
#include <cstdio>

#include <linearalgebra/Matrix3.hpp>
#include <linearalgebra/VECTOR.h>

#include <utils/IO.h>
#include <utils/trace.h>

#include <memory.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// An arbitrary dimension matrix class
//////////////////////////////////////////////////////////////////////
class MATRIX {

    public:
        MATRIX();
        MATRIX(int rows, int cols);
        MATRIX(int rows, int cols, REAL* data);
        MATRIX(int rows, int cols, REAL* data, int lda);
        MATRIX(const char* filename);
        MATRIX(const MATRIX& m);
        // matrix with "vec" along the diagonal
        MATRIX(VECTOR& vec);

        // make a copy of a 3x3 matrix
        virtual ~MATRIX();

        inline REAL& operator()(int row, int col) {
            return _matrix[row * _cols + col];
        };

        REAL operator()(int row, int col) const {
            return _matrix[row * _cols + col];
        };

        int& rows() { return _rows; };
        int& cols() { return _cols; };

        int rows() const { return _rows; };
        int cols() const { return _cols; };

        // return a pointer to the beginning of a row
        REAL* row(int index) { return &_matrix[index * _cols]; };
        const REAL* row(int index) const { return &_matrix[index * _cols]; };

        // wipe the whole matrix
        void clear();

        // Copy without reallocating
        void copyInplace( const MATRIX &A );

        // write the matrix to a binary file
        // everything is always written as a double
        void write(const char* filename);

        // read from a binary file
        // everything is always read in as a double, then
        // converted if necessary
        void read(const char* filename);

        // resize the matrix and wipe to zero
        void resizeAndWipe(int rows, int cols);

        // overload operators
        MATRIX& operator=(const MATRIX m);
        MATRIX& operator-=(const MATRIX& m);
        MATRIX& operator+=(const MATRIX& m);
        MATRIX& operator*=(const REAL& alpha);

        // return the matrix diagonal
        VECTOR diagonal();

        // return the transpose of the current matrix
        MATRIX transpose();

        // raw data pointer
        REAL* data() { return _matrix; };
        const REAL *data() const { return _matrix; }

        // 
        REAL sum()
        {
            REAL sum = 0.0;
            for( int i = 0; i < _rows*_cols; i++ )
            {
                sum += _matrix[i];
            }
            return sum;
        }

        // stomp the current matrix with the given matrix 
        // starting at row number "row". It is your responsibility
        // to ensure that you don't fall off the end of this matrix.
        void setSubmatrix(MATRIX& matrix, int row);

        // Stomp the current matrix with the given vector,
        // starting at (row,col). So the first elemenet of the vector will
        // be copied into (row,col), the next into (row+1, col), etc.
        void setVector( VECTOR& vector, int row, int col );

        void copyRowFrom( MATRIX& src, int srcRow, int row );
        void copyRowFrom( VECTOR& src, int srcRow, int row );

        // BLAS axpy operation: B += alpha * A, where B is this matrix
        //
        // Note that axpy actually applies to vectors, but in this
        // case we can just treat the matrix as a vector and multiply
        // all its elements by alpha
        void axpy(REAL alpha, MATRIX& A);

        // Try doing a parallelized version of the above, in which the
        // matrix is divided in to row blocks
        void parallelAxpy( REAL alpha, MATRIX &A );
        void parallelAxpy( MATRIX &A );
        void parallelCopy( MATRIX &A );
        void parallelCopyAdd( MATRIX &A1, MATRIX &A2, MATRIX &A3 );

        // same as axpy above, but matrix contents are stomped as well
        void clearingAxpy(REAL alpha, MATRIX& A);

        // BLAS gemm operation: C += alpha * A * B
        // where C is this matrix
        void gemm(REAL alpha, MATRIX& A, MATRIX& B);
        void gemm(MATRIX& A, MATRIX& B) { gemm(1.0f, A, B); };

        // same as gemm above, but matrix contents are stomped as well
        void clearingGemm(REAL alpha, MATRIX& A, MATRIX& B,
                bool transposeA = false, bool tranposeB = false);
        void clearingGemm(MATRIX& A, MATRIX& B,
                bool transposeA = false, bool transposeB = false)
        {
            clearingGemm(1.0f, A, B, transposeA, transposeB);
        }

        // solve the linear system Ax = b, return x in the passed in b
        void solve(VECTOR& b);
        void solveSPD(VECTOR& b); // Same thing for symmetric positive-definite systems

        // Note: This version involves taking tranposes (eg. dynamic
        // allocation).
        void solveSPD(MATRIX &B);

        // Inverts this matrix in place
        void invert();

        // Returns the inverse of this matrix in A
        void inverse( MATRIX &A );

        // multiply in place
        // this * x = y
        // Do *NOT* let x == y!
        void multiplyInplace(VECTOR& x, VECTOR& y,
                REAL alpha = 1.0, bool add = false);
        void multiplyInplace(REAL *x, REAL *y,
                REAL alpha = 1.0, bool add = false);

        void subMatrixMultiplyInplace( VECTOR& x, VECTOR& prod, int subRows,
                int subCols, bool transpose );

        // Assumes matrix is symmetric
        void uppertriMultiplyInplace( VECTOR& x, VECTOR& prod );

        // solve for eigenvalues
        void eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors);

        // Returns the Frobenius-norm of the difference between this and B
        REAL differenceFrobeniusSq( MATRIX& B );

        REAL frobeniusNorm() const;

        // Computes the solution to the least squares system A * x = b
        // where A is this matrix.
        // Results are placed in the leading part of b
        //
        // NOTE: This changes the contents of the matrix (see MKL
        //       reference manual)
        int leastSquaresFullRank( VECTOR &b, bool transpose = false );

        //////////////////////////////////////////////////////////////////////
        // Static routines
        //////////////////////////////////////////////////////////////////////

        // Least squares via truncated singular value decomposition
        //
        // Note: both matrices are overwritten by this call
        static void LeastSquares_TSVD( MATRIX &A, MATRIX &B,
                VECTOR &singularValues,
                int &rank, REAL tolerance );

        //////////////////////////////////////////////////////////////////////
        // Static routines for handling array-based matrices.
        // Use with care - nothing here does any bound checking.
        // We assume row-major format in all cases.
        //////////////////////////////////////////////////////////////////////

        // Zero out a matrix
        static inline void clear( REAL *A, int rows, int cols )
        {
            memset( (void *)A, 0, rows * cols * sizeof( REAL ) );
        }

        // Zero a matrix with a given leading dimension
        static inline void clear( REAL *A, int rows, int cols, int lda )
        {
            for ( int row_idx = 0; row_idx < rows; row_idx++ )
            {
                clear( A + row_idx * lda, 1, cols );
            }
        }

        // Generate a diagonal matrix
        static inline void diagonal( REAL *A, REAL *D, int rows )
        {
            for ( int i = 0; i < rows; i++ )
            {
                access( A, rows, rows, i, i ) = D[ i ];
            }
        }

        // Copy one matrix in to another (A <- B)
        static void copy( REAL *A, const REAL *B, int rows, int cols );

        // Copies a row (or part of it) from one matrix to another
        static void copyRow( const REAL *A, REAL *B,
                int rowA, int rowB,
                int nColsA, int nColsB,
                int nCopyCols = -1 );

        // Copies a set of rows r0, r1, r2, ... from A to rows 0, 1, 2, ... of B
        static void copyRows( const REAL *A, REAL *B,
                const IntArray &rowsA, int nCols );

        // Copies rows 0, 1, 2, of A in to r0, r1, r2, of B
        static void scatterRows( const REAL *A, REAL *B,
                const IntArray &rowsB, int nCols );

        // Access element from a matrix
        static inline REAL &access( REAL *A, int rows, int cols, int i, int j )
        {
            return A[ i * cols + j ];
        }
        static inline REAL access( const REAL *A, int rows, int cols, int i, int j )
        {
            return A[ i * cols + j ];
        }

        // Add one matrix to another (A := A + alpha * B)
        static void axpy( REAL *A, const REAL *B, int rows, int cols,
                REAL alpha = 1.0,
                int incB = 1, int incA = 1 );

        // Matrix-matrix multiplication (C := beta * C + alpha * A * B)
        static void gemm( const REAL *A, const REAL *B, REAL *C,
                int rowsA, int colsA, int rowsB, int colsB,
                bool transposeA, bool transposeB,
                REAL alpha = 1.0, REAL beta = 0.0 );

        // Version of the above requiring leading dimensions for everything
        static void gemm( const REAL *A, const REAL *B, REAL *C,
                int rowsA, int colsA, int rowsB, int colsB,
                bool transposeA, bool transposeB,
                int lda_A, int lda_B, int lda_C,
                REAL alpha = 1.0, REAL beta = 0.0 );

        // Matrix-vector multiplication (C := beta * C + alpha * A * b)
        static void gemv( const REAL *A, const REAL *b, REAL *c,
                int rowsA, int colsA,
                bool transposeA,
                REAL alpha = 1.0, REAL beta = 0.0 );

        // Symmetrix matrix-matrix update (C := beta * C + alpha * A * A')
        // Updates the lower-triangular part of C
        //
        // If trans == true then C := beta * C + alpha * A' * A
        static void syrk( const REAL *A, REAL *C,
                int rowsC, int k, /* A is either n x k or k x n */
                bool transpose = false,
                REAL alpha = 1.0, REAL beta = 0.0 );

        // Version of the above which explicitly requires leaging dimensions
        static void syrk( const REAL *A, REAL *C,
                int rowsC, int k, /* A is either n x k or k x n */
                int lda_A, int lda_C,
                bool transpose = false,
                REAL alpha = 1.0, REAL beta = 0.0 );

        // Get transpose (A = B^T)
        static void transpose( REAL *A, const REAL *B, int rows, int cols );

        static void transposeBLAS( REAL *A, const REAL *B, int rows, int cols );

        // In place transpose for square matrices
        static void transpose( REAL *A, int rows );

        // Compute eigenvalues/vectors.  We require everything here,
        // including workspaces, etc.
        // workspace should have size (7 * rows)
        // vectorWorkspace should have size (rows * rows)
        static void eigensystem( REAL *A, int rows,
                REAL *eigenvalues, REAL *eigenvectors,
                REAL *workspace, REAL *vectorWorkspace );

        // Performs QR factorization of the given system.  Other calls
        // will be needed following this to work the factorization,
        // since this is just a wrapper for the basic LAPACK call.
        //
        // No pivoting is used here
        static int qr( REAL *A, int nRows, int nCols,
                REAL *extraData, REAL *workspace,
                int workSize );

        // Extracts the leading columns of a QR factor computed
        // via the qr function.  Parameters should be the same.
        static int extractQRfactor( REAL *A, int nRows, int nCols,
                REAL *extraData, REAL *workspace,
                int workSize );

        // Computes the Cholesky factor of the matrix stored in A.  A is
        // overwritten.
        // Computes L such that A = L * L'
        static int cholesky( REAL *A, int rows );

        // Given a triangular matrix, and another input matrix, solve
        // the associated triangular system.
        //
        // By default, solves L * X = B, where L is non-unit lower-triangular
        static void triangularSolve( const REAL *L, REAL *B,
                int rowsB, int colsB,
                bool leftSide = true,
                bool lower = true,
                bool transpose = false,
                bool unitDiag = false,
                REAL alpha = 1.0 );

        // Same as the above, but for a vector
        static void triangularSolve( const REAL *L, REAL *b,
                int N,
                bool lower = true,
                bool transpose = false,
                bool unitDiag = false,
                REAL alpha = 1.0 );

        // Full rank least squares solve: A * x = b
        //
        // Overwrites A, and overwrites b with the solution
        static int leastSquaresFullRank( REAL *A, REAL *b, int nRows, int nCols,
                int nRHS = 1, bool transpose = false );

        // Least squares for an ill-conditioned matrix
        static int leastSquaresSVD( REAL *A, REAL *b, REAL *singularValues,
                int nRows, int nCols,
                REAL tolerance, int &rank, int nRHS = 1 );

        // Symmetric positive definite solver
        static int solveSPD( REAL *A, REAL *b, int nRows, int lda_A );

        // 2x2 eigensolver
        static void eigensystem2x2( REAL *A, REAL *eigenvalues, REAL *eigenvectors );

        // 3x3 low precision eigensolver
        static void eigensystem3x3( REAL *A, REAL *eigenvalues, REAL *eigenvectors );

        // 3x3 symmetric eigensolver
        static void symmetricEigensystem( const Matrix3d &M, Vector3d &eigenvalues,
                                          Matrix3d &eigenvectors );

        // Scales a matrix.  A *= alpha
        static void scale( REAL *A, int rows, int cols, REAL alpha );

    protected:
        int _rows;
        int _cols;

        REAL* _matrix;
};

// overloaded operators
VECTOR operator*(MATRIX& A, VECTOR& x);
MATRIX operator*(MATRIX& A, REAL alpha);
MATRIX operator*(MATRIX& A, MATRIX& B);
ostream& operator<<(ostream &out, MATRIX& matrix);

// multiply by the transpose of A
VECTOR operator^(MATRIX& A, VECTOR& x);
MATRIX operator^(MATRIX& A, MATRIX& B);

#endif
