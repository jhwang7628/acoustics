// SPARSE_MATRIX.h: interface for the SPARSE_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <TYPES.h>

#include "MATRIX.h"
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cstdio>

#include <memory.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// A sparse matrix class based on maps
//////////////////////////////////////////////////////////////////////
class SPARSE_MATRIX {

    public:
        SPARSE_MATRIX();
        SPARSE_MATRIX(MATRIX& matrix);
        SPARSE_MATRIX(int rows, int cols);

        // see if an entry exists in the matrix
        bool exists(int col, int row);

        // get the reference to an entry -
        // note that if an entry doesn't exist, it will be created and
        // set to zero.
        // 
        // to check if an entry already exists, use the exists() function
        REAL& operator()(int row, int col);
        REAL operator()(int row, int col) const;

        int& rows() { return _rows; };
        int& cols() { return _cols; };
        int rows() const { return _rows; }
        int cols() const { return _cols; }
        void resize(int rows, int cols) { _rows = rows; _cols = cols; };
        SPARSE_MATRIX& operator*=(const REAL& alpha);
        SPARSE_MATRIX& operator+=(const SPARSE_MATRIX& A);
        SPARSE_MATRIX& operator-=(const SPARSE_MATRIX& A);

        // Set the matrix to zero. Note this will *NOT* stomp the underlying map!
        // It will instead set all current entries to zero so that we are not
        // forced to reallocate the sparsity structure again
        void clear();

        // Clear which stomps the map structure as well
        void clearFull();

        // write the sparse matrix to Matlab format
        void writeToMatlab(string filename = string("mysparse.m"),
                string varName = string("A"));

        // write to a simple binary format
        void writeToBinary(string filename);

        bool readFromBinary( string filename );

        typedef map<pair<int,int>, REAL>::iterator iterator;
        typedef map<pair<int,int>, REAL>::const_iterator const_iterator;

        // direct access to the matrix
        map<pair<int,int>, REAL>& matrix() { return _matrix; };
        const map<pair<int,int>, REAL>& matrix() const { return _matrix; };

        iterator begin() { return _matrix.begin(); }
        iterator end() { return _matrix.end(); }

        const_iterator begin() const { return _matrix.begin(); }
        const_iterator end() const { return _matrix.end(); }

        // All the matrix entries repackaged into vectors,
        // ie _matrix(rows[x], cols[x]) = values[x];
        void entries(vector<int>& rows, vector<int>& cols, vector<REAL>& values);
        void entries(vector<int>& rows, vector<int>& cols, vector<const REAL*>& values);

        // count of non-zero elements per row -- PetSc needs this to
        // efficiently preallocate
        vector<int> nonZerosPerRow();

        // BLAS axpy operation: B += alpha * A, where B is this matrix
        //
        // Note that axpy actually applies to vectors, but in this
        // case we can just treat the matrix as a vector and multiply
        // all its elements by alpha
        void axpy(REAL alpha, SPARSE_MATRIX& A);

        // gemv operation
        //
        // y = alpha * A * x
        // (Don't let y == x)
        void gemv(REAL alpha, VECTOR &x, VECTOR &y);
        void gemv(VECTOR &x, VECTOR &y);

        // copy A into the current object
        void copies(SPARSE_MATRIX& A);
        void equals(SPARSE_MATRIX& A) { copies(A); };

        // copy A into the current object, where A(0,0) starts copying 
        // into (*this)(row,col)
        void subcopies(MATRIX& A, int row, int col);

        // add A into the current object, where A(0,0) starts copying 
        // into (*this)(row,col)
        void add(MATRIX& A, int row, int col);

        // subtract A from the current object, where A(0,0) starts copying 
        // into (*this)(row,col)
        void subtract(MATRIX& A, int row, int col);

        // print a specific row
        void printRow(int row);

        // convert to a full matrix
        MATRIX full();

        // erase a row/column
        void eraseRowColumn(int rowCol);

        virtual int size() const { return _matrix.size(); };

        REAL norm1();
        REAL normInf();
        REAL normFrob();

        // Returns the sum of all entries. Mainly for debugging.
        REAL sum() const;

        // Performs out = basis^T * this * basis
        // workMat should be same size as basis, but transposed
        void projectInplace( MATRIX& basis, MATRIX& workMat, MATRIX& out );

        // Multiplies the given dense matrix on the left by this, and places
        // the result in out.
        // Sets out = this * A
        void denseMultiply( MATRIX &A, MATRIX &out );

        // Converts this to a dense matrix
        MATRIX dense();

        // Writes a text representation of this matrix
        void writeText( std::ostream &os ) const
        {
            os << setprecision(12);
            os << "A = SPARSE_MATRIX( " << rows() << ", " << cols() << " );" << endl << endl;
            for ( map<pair<int, int>, REAL>::const_iterator iter = _matrix.begin();
                  iter != _matrix.end(); iter++ ) {
                os << "A( " << iter->first.first << ", " << iter->first.second
                   << " ) = " << iter->second << ";" << endl;
            }
        }

        // Sparse column matrix representation
        // This representation follows from the one used by CHOLMOD,
        // except that we always represent the full matrix, not just
        // the lower or upper triangular part.
        struct SparseColumnMatrix {
            SparseColumnMatrix()
                : _nrow( 0 ),
                  _ncol( 0 ),
                  _nzmax( 0 ),
                  _p( NULL ),
                  _i( NULL ),
                  _x( NULL )
            {
            }

            SparseColumnMatrix( SparseColumnMatrix &M )
                : _nrow( 0 ),
                  _ncol( 0 ),
                  _nzmax( 0 ),
                  _p( NULL ),
                  _i( NULL ),
                  _x( NULL )
            {
                copy( M );
            }

            virtual ~SparseColumnMatrix()
            {
                clear();
            }

            SparseColumnMatrix &operator=( SparseColumnMatrix &M )
            {
                copy( M );
                return *this;
            }

            void clear()
            {
                if ( _p != NULL ) {
                    free( _p );
                }

                if ( _i != NULL ) {
                    free( _i );
                }

                if ( _x != NULL ) {
                    free( _x );
                } 

                _p = NULL;
                _i = NULL;
                _x = NULL;

                _nrow = 0;
                _ncol = 0;
                _nzmax = 0;
            }

            size_t usage() const
            {
                size_t totalUsage = 0;

                if ( _p != NULL ) {
                    totalUsage += (_ncol + 1) * sizeof( int );
                }

                if ( _i != NULL ) {
                    totalUsage += _nzmax * sizeof( int );
                }

                if ( _x != NULL ) {
                    totalUsage += _nzmax * sizeof( REAL );
                }

                return totalUsage;
            }

            // Counts non-zeros in the lower triangular part of the matrix
            size_t symmetricNonZeros() const
            {
                size_t totalNZ = 0;

                for ( int col_idx = 0; col_idx < (int)_ncol; col_idx++ ) {
                    for ( int row_ptr = _p[ col_idx ]; row_ptr < _p[ col_idx + 1 ];
                            row_ptr++ )
                    {
                        if ( _i[ row_ptr ] >= col_idx ) {
                            totalNZ++;
                        }
                    }
                }

                return totalNZ;
            }

            REAL usageMB() const
            {
                return (REAL)usage() / 1024.0 / 1024.0;
            }

            // Assumes size bounds have been set
            void allocate()
            {
                _p = (int *)malloc( (_ncol + 1) * sizeof( int ) );
                _i = (int *)malloc( _nzmax * sizeof( int ) );
                _x = (REAL *)malloc( _nzmax * sizeof( REAL ) );

                memset( (void *)_p, 0, sizeof( int ) * ( _ncol + 1 ) );
                memset( (void *)_i, 0, sizeof( int ) * _nzmax );
                memset( (void *)_x, 0, sizeof( REAL ) * _nzmax );
            }

            void copy( SparseColumnMatrix &M )
            {
                clear();

                _nrow = M._nrow;
                _ncol = M._ncol;
                _nzmax = M._nzmax;

                // Copy the fields from M
                if ( M._x != NULL )
                {
                    _p = (int *)malloc( (_ncol + 1) * sizeof( int ) );
                    _i = (int *)malloc( _nzmax * sizeof( int ) );
                    _x = (REAL *)malloc( _nzmax * sizeof( REAL ) );

                    memcpy( (void *)_p, (void *)M._p, (_ncol + 1) * sizeof( int ) );
                    memcpy( (void *)_i, (void *)M._i, _nzmax * sizeof( int ) );
                    memcpy( (void *)_x, (void *)M._x, _nzmax * sizeof( REAL ) );
                }
            }

            size_t           _nrow;   // Number of rows
            size_t           _ncol;   // Number of columns
            size_t           _nzmax;  // Number of nonzeros

            int             *_p;      // [0..ncol] Column pointers
            int             *_i;      // [0..nzmax-1] Row indices
            REAL            *_x;      // Actual matrix data
        };

        // Copies the matrix in to sparse column form
        void constructSparseColumnCopy( SparseColumnMatrix &M ) const;

        // Multiplies the matrix with the given vector
        static int matrixMultiply( const SparseColumnMatrix &A,
                                   const REAL *g, REAL *b,
                                   REAL alpha,
                                   bool clear = true, bool transpose = false,
                                   int ldaG = -1, int ldaB = -1 );

        // Multiplies a sub matrix of the given sparse column matrix
        // with the given matrix G, placing the result in B
        // Eg. this forms B = A * G (or B = A' * G)
        static int subMatrixMultiply( const SparseColumnMatrix &A,
                                      const REAL *G, REAL *B,
                                      int rowStart, int columnStart,
                                      int nRows, int nCols, int nColsG,
                                      REAL alpha,
                                      bool clear = true, bool transpose = false,
                                      int ldaG = -1, int ldaB = -1 );

        // Same as the above, but with vector inputs g and b
        static void subMatrixMultiply( const SparseColumnMatrix &A,
                                       const REAL *g, REAL *b,
                                       int rowStart, int columnStart,
                                       int nRows, int nCols,
                                       REAL alpha,
                                       bool clear = true, bool transpose = false,
                                       int ldaG = -1, int ldaB = -1 );

        // Version in which all rows after the given start rows in the
        // provided column set are used
        static void subMatrixMultiply( const SparseColumnMatrix &A,
                                       const REAL *g, REAL *b,
                                       int rowStart, int columnStart,
                                       const IntArray &startRows,
                                       REAL alpha,
                                       bool clear = true, bool transpose = false,
                                       int ldaG = -1 );

        // Version of the above with an extra workspace for faster vectorized
        // operations
        static void subMatrixMultiply( const SparseColumnMatrix &A,
                                       const REAL *g, REAL *b,
                                       int rowStart, int columnStart,
                                       const IntArray &startRows,
                                       REAL alpha,
                                       REAL *workspace,
                                       bool clear = true, bool transpose = false,
                                       int ldaG = -1 );

        // Same as the above, but uses a specific row set
        static int subMatrixMultiply( const SparseColumnMatrix &A,
                                      const REAL *G, REAL *B,
                                      int rowStart, int columnStart,
                                      int nRows, int nCols, int nColsG,
                                      const IntArray &inverseRowList,
                                      REAL alpha,
                                      bool clear = true, bool transpose = false,
                                      int ldaG = -1, int ldaB = -1 );

        // Similar to the above, but forms B = G' * A (or B = G' * A')
        static void subMatrixLeftMultiply( const SparseColumnMatrix &A,
                                           const REAL *G, REAL *B,
                                           int rowStart, int columnStart,
                                           int nRows, int nCols, int nRowsG,
                                           bool clear = true, bool transpose = false,
                                           int ldaG = -1, int ldaB = -1 );

        // Converts a sparse column matrix to it's graph laplacian
        // NOTE: Assumes that the matrix is symmetric
        static void convertToGraphLaplacian( SparseColumnMatrix &A,
                                             bool normalize );

    protected:
        int _rows;
        int _cols;

        // a dud REAL to pass back if the index is out of bounds
        REAL _dud;

        map<pair<int,int>, REAL> _matrix;

    public:

        // Creates a random sparse matrix, with random sparsity, with random
        // values in [0,1]
        static SPARSE_MATRIX* makeRandom( int rows, int cols, int nnz );

        // Writes a sparse column matrix to disk
        static void writeToBinary( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                   const char *filename,
                                   bool newVersion = false );

        // Reads a sparse column matrix from disk
        static void readFromBinary( SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const char *filename );

        // Computes Gaussian elimination fill-in for a matrix stored in
        // sparse column format.
        //
        // Returns the number of non-zeros in the factor matrix
        static long int computeFillIn( const SPARSE_MATRIX::SparseColumnMatrix &A );

};

ostream& operator<<(ostream &out, const SPARSE_MATRIX& matrix);
VECTOR operator*(const SPARSE_MATRIX& A, const VECTOR& x);
VECTOR operator*(const SPARSE_MATRIX::SparseColumnMatrix &A, const VECTOR &x);
SPARSE_MATRIX operator*(const SPARSE_MATRIX& A, const REAL& alpha);
MATRIX operator*(const SPARSE_MATRIX& A, const MATRIX& B);

#endif
