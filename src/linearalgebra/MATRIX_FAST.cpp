// MATRIX.h: interface for the MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include "MATRIX.h"

#ifdef USE_MKL
#include <mkl_lapack.h>
#ifdef USE_LAPACKE
#include <mkl_lapacke.h>
#endif
#include <mkl_types.h>
#include <mkl_cblas.h>
#endif
#ifdef USING_ATLAS
extern "C" {
#include <clapack.h>
#include <cblas.h>
}
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif



#include <utils/trace.h>

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
MATRIX::MATRIX() :
    _rows(0), _cols(0)
{
    _matrix = NULL;
}

MATRIX::MATRIX(int rows, int cols) :
    _rows(rows), _cols(cols)
{
    _matrix = new REAL[_rows * _cols];
    clear();
}

MATRIX::MATRIX(int rows, int cols, REAL* data) :
    _rows(rows), _cols(cols)
{
    _matrix = new REAL[_rows * _cols];
    memcpy(_matrix, data, _rows * _cols * sizeof(REAL));
}

MATRIX::MATRIX(int rows, int cols, REAL* data, int lda) :
    _rows(rows), _cols(cols)
{
    _matrix = new REAL[_rows * _cols];

    for ( int row_idx = 0; row_idx < rows; row_idx++ )
    {
        MATRIX::copyRow( data, _matrix, row_idx, row_idx, lda, cols, cols );
    }
}

MATRIX::MATRIX(const char* filename)
{
    _matrix = NULL;
    read(filename);
}

MATRIX::MATRIX(const MATRIX& m)
{
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new REAL[_rows * _cols];
    memcpy(_matrix, m._matrix, _rows * _cols * sizeof(REAL));
}

MATRIX::MATRIX(VECTOR& vec)
{
    _cols = vec.size();
    _rows = vec.size();

    _matrix = new REAL[_rows * _cols];
    clear();

    for (int x = 0; x < vec.size(); x++)
        (*this)(x,x) = vec(x);
}

MATRIX::~MATRIX()
{
    delete[] _matrix;
}

//////////////////////////////////////////////////////////////////////
// wipe the whole matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::clear()
{
    memset(_matrix, 0, _rows * _cols * sizeof(REAL));
}

//////////////////////////////////////////////////////////////////////
// Copy without reallocating
//////////////////////////////////////////////////////////////////////
void MATRIX::copyInplace( const MATRIX &A )
{
    TRACE_ASSERT( _rows == A.rows() && _cols == A.cols() );

    MATRIX::copy( this->_matrix, A._matrix, _rows, _cols );
}

//////////////////////////////////////////////////////////////////////
// resize and wipe the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::resizeAndWipe(int rows, int cols)
{
    if (_rows != rows || _cols != cols)
        delete[] _matrix;
    _rows = rows;
    _cols = cols;

    _matrix = new REAL[_rows * _cols];
    clear();
}

//////////////////////////////////////////////////////////////////////
// write the matrix to a file
//////////////////////////////////////////////////////////////////////
void MATRIX::write(const char* filename)
{
    FILE* file;
    file = fopen(filename, "wb");

    // write dimensions
    fwrite((void*)&_rows, sizeof(int), 1, file);
    fwrite((void*)&_cols, sizeof(int), 1, file);

    // always write out as a double
    if (sizeof(REAL) != sizeof(double))
    {
        double* matrixDouble = new double[_rows * _cols];
        for (int x = 0; x < _rows * _cols; x++)
            matrixDouble[x] = _matrix[x];

        fwrite((void*)matrixDouble, sizeof(double), _rows * _cols, file);
        delete[] matrixDouble;
        fclose(file);
    }
    else
        fwrite((void*)_matrix, sizeof(REAL), _rows * _cols, file);
    fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
void MATRIX::read(const char* filename)
{
    FILE* file;
    file = fopen(filename, "rb");
    if (file == NULL)
    {
        cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
        return;
    }

    // read dimensions
    fread((void*)&_rows, sizeof(int), 1, file);
    fread((void*)&_cols, sizeof(int), 1, file);

    // read data
    if( _matrix != NULL ) delete[] _matrix;

    _matrix = new REAL[_rows * _cols];

    // always read in a double
    if (sizeof(REAL) != sizeof(double))
    {
        double* matrixDouble = new double[_rows * _cols];
        fread((void*)matrixDouble, sizeof(double), _rows * _cols, file);
        for (int x = 0; x < _rows * _cols; x++)
            _matrix[x] = matrixDouble[x];
        delete[] matrixDouble;
    }
    else 
        fread((void*)_matrix, sizeof(REAL), _rows * _cols, file);
    fclose(file);
}

//////////////////////////////////////////////////////////////////////
// return transpose of current matrix
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::transpose()
{
    MATRIX toReturn(_cols, _rows);

    for (int y = 0; y < _cols; y++)
        for (int x = 0; x < _rows; x++)
            toReturn(y,x) = (*this)(x,y);

    return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(MATRIX& A, VECTOR& x) 
{
    VECTOR y(A.rows());

#ifdef SINGLE_PRECISION
    cblas_sgemv(
            CblasRowMajor, 
            CblasNoTrans, 
            A.rows(), A.cols(),
            1.0f, A.data(), A.cols(), x.data(), 1, 
            0.0f, y.data(), 1);
#else
    cblas_dgemv(
            CblasRowMajor, 
            CblasNoTrans, 
            A.rows(), A.cols(),
            1.0, A.data(), A.cols(), x.data(), 1, 
            0.0, y.data(), 1);
#endif

    return y;
}

//////////////////////////////////////////////////////////////////////
// Scale matrix
//////////////////////////////////////////////////////////////////////
MATRIX operator*(MATRIX& A, REAL alpha) 
{
    MATRIX y(A.rows(), A.cols());

    for (int i = 0; i < A.rows(); i++)
        for (int j = 0; j < A.cols(); j++)
            y(i,j) = A(i, j) * alpha;

    return y;
}

//////////////////////////////////////////////////////////////////////
// scale matrix by a constant
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator*=(const REAL& alpha)
{
#ifdef SINGLE_PRECISION
    cblas_sscal(_cols * _rows, alpha, _matrix, 1);
#else
    cblas_dscal(_cols * _rows, alpha, _matrix, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
// Matrix-Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(MATRIX& A, MATRIX& B) 
{
    MATRIX y(A.rows(), B.cols());

#ifdef SINGLE_PRECISION
    cblas_sgemm(
            CblasRowMajor,
            CblasNoTrans,
            CblasNoTrans,
            A.rows(), 
            B.cols(),
            A.cols(), 
            1.0f,
            A.data(), A.cols(),
            B.data(), B.cols(),
            0.0f,
            y.data(), y.cols());
#else
    cblas_dgemm(
            CblasRowMajor,
            CblasNoTrans,
            CblasNoTrans,
            A.rows(), 
            B.cols(),
            A.cols(), 
            1.0,
            A.data(), A.cols(),
            B.data(), B.cols(),
            0.0,
            y.data(), y.cols());
#endif
    return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply where A is transposed
//////////////////////////////////////////////////////////////////////
VECTOR operator^(MATRIX& A, VECTOR& x) 
{
    VECTOR y(A.cols());

#ifdef SINGLE_PRECISION
    cblas_sgemv(
            CblasRowMajor, 
            CblasTrans, 
            A.rows(), A.cols(), 
            1.0f, A.data(), A.cols(), x.data(), 1, 
            0.0f, y.data(), 1);
#else
    cblas_dgemv(
            CblasRowMajor, 
            CblasTrans, 
            A.rows(), A.cols(), 
            1.0, A.data(), A.cols(), x.data(), 1, 
            0.0, y.data(), 1);
#endif
    return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix^T -Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator^(MATRIX& A, MATRIX& B) 
{
    MATRIX y(A.cols(), B.cols());

#ifdef SINGLE_PRECISION
    cblas_sgemm(
            CblasRowMajor,
            CblasTrans,
            CblasNoTrans,
            A.cols(), 
            B.cols(), 
            A.rows(),
            1.0f,
            A.data(), A.cols(),
            B.data(), B.cols(),
            0.0f,
            y.data(), y.cols());
#else
    cblas_dgemm(
            CblasRowMajor,
            CblasTrans,
            CblasNoTrans,
            A.cols(), 
            B.cols(), 
            A.rows(),
            1.0,
            A.data(), A.cols(),
            B.data(), B.cols(),
            0.0,
            y.data(), y.cols());
#endif
    return y;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, MATRIX& matrix)
{
    out << "[" << endl; 
    for (int row = 0; row < matrix.rows(); row++)
    {
        for (int col = 0; col < matrix.cols(); col++)
            out << matrix(row, col) << " ";
        out << endl;
    }
    out << "]" << endl;
    return out;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator=(const MATRIX m)
{
    if (m._cols != _cols || m._rows != _rows)
    {
        delete[] _matrix;
        _cols = m._cols;
        _rows = m._rows;

        _matrix = new REAL[_rows * _cols];
    }
    memcpy(_matrix, m._matrix, _rows * _cols * sizeof(REAL));
    return *this;
}

//////////////////////////////////////////////////////////////////////
// self minus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator-=(const MATRIX& m)
{
    if (m._cols != _cols || m._rows != _rows)
    {
        delete[] _matrix;
        _cols = m._cols;
        _rows = m._rows;

        _matrix = new REAL[_rows * _cols];
    }
#ifdef SINGLE_PRECISION
    cblas_saxpy(_cols * _rows, -1.0f, m._matrix, 1, _matrix, 1);
#else
    cblas_daxpy(_cols * _rows, -1.0, m._matrix, 1, _matrix, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
// self plus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator+=(const MATRIX& m)
{
    if (m._cols != _cols || m._rows != _rows)
    {
        delete[] _matrix;
        _cols = m._cols;
        _rows = m._rows;

        _matrix = new REAL[_rows * _cols];
    }
#ifdef SINGLE_PRECISION
    cblas_saxpy(_cols * _rows, 1.0f, m._matrix, 1, _matrix, 1);
#else
    cblas_daxpy(_cols * _rows, 1.0, m._matrix, 1, _matrix, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
// Return the matrix diagonal
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::diagonal()
{
    int minDim = (_rows > _cols) ? _cols : _rows;
    VECTOR diag(minDim);
    for (int x = 0; x < minDim; x++)
        diag(x) = (*this)(x,x);

    return diag;
}

//////////////////////////////////////////////////////////////////////
// stomp the current matrix with the given matrix starting at "row". 
// It is your responsibility to ensure that you don't fall off the 
// end of this matrix.
//////////////////////////////////////////////////////////////////////
void MATRIX::setSubmatrix(MATRIX& matrix, int row)
{
    int totalSize = matrix.rows() * matrix.cols();
    int index = row * _cols;

    for (int x = 0; x < totalSize; x++, index++)
        _matrix[index] = matrix._matrix[x];
}

//////////////////////////////////////////////////////////////////////
// Stomp the current matrix with the given vector,
// starting at (row,col). So the first elemenet of the vector will
// be copied into (row,col), the next into (row+1, col), etc.
//////////////////////////////////////////////////////////////////////
void MATRIX::setVector( VECTOR& vector, int row, int col )
{
    for( int j = 0; j < vector.size(); j++ ) {
        (*this)(row+j, col) = vector(j);
    }
}

//////////////////////////////////////////////////////////////////////
// This assumes row-major storage
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRowFrom( MATRIX& src, int srcRow, int row )
{
    memcpy( &(*this)(row,0), &src(srcRow,0), sizeof(REAL)*_cols );
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy(REAL alpha, MATRIX& A)
{
#ifdef SINGLE_PRECISION
    cblas_saxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#else
    cblas_daxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Try doing a parallelized version of the above, in which the
// matrix is divided in to row blocks
//////////////////////////////////////////////////////////////////////
void MATRIX::parallelAxpy( REAL alpha, MATRIX &A )
{
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( int row_idx = 0; row_idx < _rows; row_idx += 1 )
    {
        for ( int col_idx = 0; col_idx < _cols; col_idx++ )
        {
            this->operator()( row_idx, col_idx ) += alpha * A( row_idx, col_idx );
        }
#if 0
        REAL                    *output = _matrix + row_idx * _cols;
        REAL                    *input = A._matrix + row_idx * _cols;

        int                      nRows = min( _rows - row_idx, BLOCK_SIZE );

        axpy( output, input, BLOCK_SIZE, _cols, alpha );
#endif
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::parallelAxpy( MATRIX &A )
{
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( int row_idx = 0; row_idx < _rows; row_idx += 1 )
    {
        for ( int col_idx = 0; col_idx < _cols; col_idx++ )
        {
            this->operator()( row_idx, col_idx ) += A( row_idx, col_idx );
        }
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::parallelCopy( MATRIX &A )
{
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( int row_idx = 0; row_idx < _rows; row_idx += 1 )
    {
        for ( int col_idx = 0; col_idx < _cols; col_idx++ )
        {
            this->operator()( row_idx, col_idx ) = A( row_idx, col_idx );
        }
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::parallelCopyAdd( MATRIX &A1, MATRIX &A2, MATRIX &A3 )
{
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( int row_idx = 0; row_idx < _rows; row_idx += 1 )
    {
        for ( int col_idx = 0; col_idx < _cols; col_idx++ )
        {
            this->operator()( row_idx, col_idx ) = A1( row_idx, col_idx )
                + A2( row_idx, col_idx )
                + A3( row_idx, col_idx );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B = alpha * A, where B is this matrix, and 
// current contents of B are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingAxpy(REAL alpha, MATRIX& A)
{
    memset(_matrix, 0, _rows * _cols * sizeof(REAL));
#ifdef SINGLE_PRECISION
    cblas_saxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#else
    cblas_daxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C += alpha * A * B where C is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm(REAL alpha, MATRIX& A, MATRIX& B)
{
#ifdef SINGLE_PRECISION
    cblas_sgemm(
            CblasRowMajor,
            CblasNoTrans,
            CblasNoTrans,
            A.rows(), 
            B.cols(),
            A.cols(), 
            1.0f,
            A.data(), A.cols(),
            B.data(), B.cols(),
            0.0f,
            _matrix, _cols);
#else
    cblas_dgemm(
            CblasRowMajor,
            CblasNoTrans,
            CblasNoTrans,
            A.rows(), 
            B.cols(),
            A.cols(), 
            1.0,
            A.data(), A.cols(),
            B.data(), B.cols(),
            0.0,
            _matrix, _cols);
#endif
}

/*
 * Untested -- don't uncomment until a test case comes up
 * 
//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(VECTOR& x)
{
if (x.size() != _rows)
cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dimensions do not match!" << endl;

VECTOR y(x.size());
for (int j = 0; j < _cols; j++)
for (int i = 0; i < _rows; i++)
y(j) += (*this)(i,j) * x(j);

return y;
}
 */

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C = alpha * A * B where C is this matrix and
// current contents of C are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingGemm(REAL alpha, MATRIX& A, MATRIX& B,
        bool transposeA, bool transposeB)
{
#ifdef SINGLE_PRECISION
    cblas_sgemm(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            transposeB ? CblasTrans : CblasNoTrans,
            transposeA ? A.cols() : A.rows(), 
            transposeB ? B.rows() : B.cols(),
            transposeA ? A.rows() : A.cols(), 
            alpha,
            A.data(),
            A.cols(),
            //transposeA ? A.rows() : A.cols(),
            B.data(),
            B.cols(),
            //transposeB ? B.rows() : B.cols(),
            0.0f,
            _matrix,
            _cols); /// Should be rows if tranposed B?
    //transposeB ? _rows : _cols); /// Should be rows if tranposed B?
#else
    cblas_dgemm(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            transposeB ? CblasTrans : CblasNoTrans,
            transposeA ? A.cols() : A.rows(), 
            transposeB ? B.rows() : B.cols(),
            transposeA ? A.rows() : A.cols(), 
            alpha,
            A.data(),
            A.cols(),
            //transposeA ? A.rows() : A.cols(),
            B.data(),
            B.cols(),
            //transposeB ? B.rows() : B.cols(),
            0.0,
            _matrix,
            _cols); /// Should be rows if tranposed B?
    //transposeB ? _rows : _cols); /// Should be rows if tranposed B?
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//
// Do *NOT* let x == y!
//////////////////////////////////////////////////////////////////////
void MATRIX::multiplyInplace(VECTOR& x, VECTOR& y, REAL alpha, bool add) 
{
#ifdef SINGLE_PRECISION
    cblas_sgemv(
            CblasRowMajor, 
            CblasNoTrans,
            _rows, _cols,
            //1.0f, _matrix, _cols, x.data(), 1, 
            alpha, _matrix, _cols, x.data(), 1, 
            add ? 1.0f : 0.0f,
            y.data(), 1);
#else
    cblas_dgemv(
            CblasRowMajor, 
            CblasNoTrans,
            _rows, _cols,
            //1.0, _matrix, _cols, x.data(), 1, 
            alpha, _matrix, _cols, x.data(), 1, 
            add ? 1.0 : 0.0,
            y.data(), 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//
// Do *NOT* let x == y!
//////////////////////////////////////////////////////////////////////
void MATRIX::multiplyInplace(REAL *x, REAL *y, REAL alpha, bool add)
{
#ifdef SINGLE_PRECISION
    cblas_sgemv(
            CblasRowMajor, 
            CblasNoTrans,
            _rows, _cols,
            //1.0f, _matrix, _cols, x, 1, 
            alpha, _matrix, _cols, x, 1, 
            add ? 1.0f : 0.0f,
            y, 1);
#else
    cblas_dgemv(
            CblasRowMajor, 
            CblasNoTrans,
            _rows, _cols,
            //1.0, _matrix, _cols, x, 1, 
            alpha, _matrix, _cols, x, 1, 
            add ? 1.0 : 0.0,
            y, 1);
#endif
}

void MATRIX::subMatrixMultiplyInplace( VECTOR& x, VECTOR& prod, int subRows, int subCols, bool transpose )
{
#ifdef SINGLE_PRECISION
    cblas_sgemv(
            CblasRowMajor, 
            transpose ? CblasTrans : CblasNoTrans,
            subRows, subCols,
            1.0f, _matrix, _cols, x.data(), 1, 
            0.0f, prod.data(), 1);
#else
    cblas_dgemv(
            CblasRowMajor, 
            transpose ? CblasTrans : CblasNoTrans,
            subRows, subCols,
            1.0, _matrix, _cols, x.data(), 1, 
            0.0, prod.data(), 1);
#endif
}

void MATRIX::uppertriMultiplyInplace( VECTOR& x, VECTOR& prod )
{
#ifdef SINGLE_PRECISION
    cblas_ssymv(
            CblasRowMajor, 
            CblasUpper, _rows,
            1.0, _matrix, _cols, x.data(), 1,
            0.0, prod.data(), 1 );
#else
    cblas_dsymv(
            CblasRowMajor, 
            CblasUpper, _rows,
            1.0, _matrix, _cols, x.data(), 1,
            0.0, prod.data(), 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
// Here, A is a general invertible system.
//////////////////////////////////////////////////////////////////////
void MATRIX::solve(VECTOR& b)
{
    char uplo = 'U';
    int nrhs = 1;
    int info = 0;
    int R = _rows;
    int ipiv[_rows];

#ifdef USE_MKL  
    // call the general LU solver from MKL
#ifdef SINGLE_PRECISION
    sgesv(&R, &nrhs, _matrix, &R, ipiv, b.data(), &R, &info);
#else
    dgesv(&R, &nrhs, _matrix, &R, ipiv, b.data(), &R, &info);
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
// Here, A is a general invertible system.
//////////////////////////////////////////////////////////////////////
void MATRIX::solveSPD(VECTOR& b)
{
    char uplo = 'U';
    int nrhs = 1;
    int info = 0;
    int R = _rows;

#ifdef USE_MKL  
    // call the SPD solver from MKL
#ifdef SINGLE_PRECISION
    sposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
#else
    dposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
#endif
#endif
#ifdef USING_ATLAS
    // call the general LU solver from ATLAS
#ifdef SINGLE_PRECISION
    clapack_sposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#else
    clapack_dposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system AX = B, return x in the passed in b
// Here, A is a general invertible system.
//
// Note: This version involves taking tranposes (eg. dynamic
// allocation).
//////////////////////////////////////////////////////////////////////
void MATRIX::solveSPD(MATRIX &B)
{
    char uplo = 'U';
    int nrhs = B.cols();
    int info = 0;
    int R = _rows;

    // Create column major ordering of data
    REAL *rhsData = new REAL[ B.rows() * B.cols() ];

    int idx = 0;
    for ( int j = 0; j < B.cols(); j++ )
    {
        for ( int i = 0; i < B.rows(); i++ )
        {
            rhsData[ idx ] = B( i, j );
            idx++;
        }
    }

#ifdef USE_MKL  
    // call the SPD solver from MKL
#ifdef SINGLE_PRECISION
    sposv(&uplo, &R, &nrhs, _matrix, &R, rhsData, &R, &info);
#else
    dposv(&uplo, &R, &nrhs, _matrix, &R, rhsData, &R, &info);
#endif
#endif
#ifdef USING_ATLAS
    // call the general LU solver from ATLAS
#ifdef SINGLE_PRECISION
    clapack_sposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#else
    clapack_dposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#endif
#endif

    if ( info < 0 )
    {
        cout << "ERROR: dposv failed!  Invalid argument" << endl;
        abort();
    }
    else if ( info > 0 )
    {
        cout << "ERROR: dposv failed!  Matrix not positive definite" << endl;
        abort();
    }

    // Copy back in to B
    idx = 0;
    for ( int j = 0; j < B.cols(); j++ )
    {
        for ( int i = 0; i < B.rows(); i++ )
        {
            B( i, j ) = rhsData[ idx ];
            idx++;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Inverts this matrix in place
//////////////////////////////////////////////////////////////////////
void MATRIX::invert()
{
    // basic error checking
    if (_rows != _cols)
    {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Matrix must be square to be inverted! " << endl;
        return;
    }

    int R = _rows;
    int ipiv[_rows];
    int info = 0;
    REAL work_query;
    int l_work = -1;

#ifdef USE_MKL
#ifdef SINGLE_PRECISION
    sgetrf( &R, &R, _matrix, &R, ipiv, &info );

    // Query for workspace size
    sgetri( &R, _matrix, &R, ipiv, &work_query, &l_work, &info );

    l_work = (int)work_query;
    REAL work[l_work];

    sgetri( &R, _matrix, &R, ipiv, work, &l_work, &info );
#else
    dgetrf( &R, &R, _matrix, &R, ipiv, &info );

    // Query for workspace size
    dgetri( &R, _matrix, &R, ipiv, &work_query, &l_work, &info );

    l_work = (int)work_query;
    REAL work[l_work];

    dgetri( &R, _matrix, &R, ipiv, work, &l_work, &info );
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// Returns the inverse of this matrix in A
//////////////////////////////////////////////////////////////////////
void MATRIX::inverse( MATRIX &A )
{
    // basic error checking
    if (_rows != _cols)
    {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Matrix must be square to be inverted! " << endl;
        return;
    }

    A = *this;

    int R = _rows;
    int ipiv[_rows];
    int info = 0;
    REAL work_query;
    int l_work = -1;

#ifdef USE_MKL
#ifdef SINGLE_PRECISION
    sgetrf( &R, &R, A._matrix, &R, ipiv, &info );

    // Query for workspace size
    sgetri( &R, A._matrix, &R, ipiv, &work_query, &l_work, &info );

    l_work = (int)work_query;
    REAL work[l_work];

    sgetri( &R, A._matrix, &R, ipiv, work, &l_work, &info );
#else
    dgetrf( &R, &R, A._matrix, &R, ipiv, &info );

    // Query for workspace size
    dgetri( &R, A._matrix, &R, ipiv, &work_query, &l_work, &info );

    l_work = (int)work_query;
    REAL work[l_work];

    dgetri( &R, A._matrix, &R, ipiv, work, &l_work, &info );
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// solve for the eigensystem of the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
    // basic error checking
    if (_rows != _cols)
    {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Matrix must be square to get eigenvalues! " << endl;
        return;
    }

    // resize result space
    eigenvalues.resizeAndWipe(_rows);
    eigenvectors.resizeAndWipe(_rows, _rows);

    int rowsize = _rows;
    int worksize = 5 * _rows;
    REAL* work = new REAL[worksize];
    REAL* valuesREAL = new REAL[2 * _rows];
    REAL* valuesImag = valuesREAL + _rows;
    REAL* vectors = new REAL[_rows * _rows];

    // the actual LAPACK call
    int error;
#ifdef SINGLE_PRECISION
    sgeev("V","N", &rowsize, _matrix, &rowsize, 
            valuesREAL, valuesImag, 
            vectors, &rowsize, NULL, &rowsize,
            work, &worksize, &error);
#else
    dgeev("V","N", &rowsize, _matrix, &rowsize, 
            valuesREAL, valuesImag, 
            vectors, &rowsize, NULL, &rowsize,
            work, &worksize, &error);
#endif

    if (error != 0)
    {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " eigenvalue solver bombed!" << endl;
    }

    // copy out results
    for (int x = 0; x < _rows; x++)
        eigenvalues(x) = valuesREAL[x];

    for (int x = 0; x < _rows; x++)
        for (int y = 0; y < _rows; y++)
            eigenvectors(x,y) = vectors[x + y * _cols];

    // cleanup
    delete[] work;
    delete[] valuesREAL;
    delete[] vectors;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL MATRIX::differenceFrobeniusSq( MATRIX& B )
{
    // Not MKL optimized
    REAL frobSq = 0.0;
    for( int i = 0; i < _rows; i++ ) {
        for( int j = 0; j < _cols; j++ ) {
            REAL diff = (*this)(i,j) - B(i,j);
            frobSq += diff*diff;
        }
    }
    return frobSq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL MATRIX::frobeniusNorm() const
{
#ifdef SINGLE_PRECISION
    return cblas_snrm2( _rows * _cols, _matrix, 1 );
#else
    return cblas_dnrm2( _rows * _cols, _matrix, 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// Computes the solution to the least squares system A * x = b
// where A is this matrix.
// Results are placed in the leading part of b
//
// NOTE: This changes the contents of the matrix (see MKL
//       reference manual)
//////////////////////////////////////////////////////////////////////
int MATRIX::leastSquaresFullRank( VECTOR &b, bool transpose )
{
    int                        info;

#ifdef USE_LAPACKE
#ifdef SINGLE_PRECISION
    info = LAPACKE_sgels( LAPACK_ROW_MAJOR,
            transpose ? 'T' : 'N',
            _rows, _cols, 1,
            _matrix, _cols,
            b.data(), 1 );
#else
    info = LAPACKE_dgels( LAPACK_ROW_MAJOR,
            transpose ? 'T' : 'N',
            _rows, _cols, 1,
            _matrix, _cols,
            b.data(), 1 );
#endif
#else
    cerr << "NOT IMPLEMENTED!" << endl;
    abort();
#endif

    return info;
}

//////////////////////////////////////////////////////////////////////
// Least squares via truncated singular value decomposition
//
// Note: both matrices are overwritten by this call
//////////////////////////////////////////////////////////////////////
void MATRIX::LeastSquares_TSVD( MATRIX &A, MATRIX &B,
        VECTOR &singularValues,
        int &rank, REAL tolerance )
{
    singularValues.resizeAndWipe( min( A.rows(), A.cols() ) );

    leastSquaresSVD( A.data(), B.data(), singularValues.data(),
            A.rows(), A.cols(), tolerance, rank, B.cols() );
}

//////////////////////////////////////////////////////////////////////
// Static routines for handling array-based matrices.
// Use with case - nothing here does any bound checking
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Copy one matrix in to another (A <- B)
//////////////////////////////////////////////////////////////////////
void MATRIX::copy( REAL *A, const REAL *B, int rows, int cols )
{
    memcpy( A, B, rows * cols * sizeof( REAL ) );
}

//////////////////////////////////////////////////////////////////////
// Copies a row (or part of it) from one matrix to another
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRow( const REAL *A, REAL *B,
        int rowA, int rowB,
        int nColsA, int nColsB,
        int nCopyCols )
{
    if ( nCopyCols < 0 )
    {
        nCopyCols = nColsA;
    }

    memcpy( (void *)( B + rowB * nColsB ), (void *)( A + rowA * nColsA ),
            nCopyCols * sizeof( REAL ) );
}

//////////////////////////////////////////////////////////////////////
// Copies a set of rows r0, r1, r2, ... from A to rows 0, 1, 2, ... of B
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRows( const REAL *A, REAL *B,
        const IntArray &rowsA, int nCols )
{
    for ( int row_idx = 0; row_idx < rowsA.size(); row_idx++ )
    {
        MATRIX::copyRow( A, B, rowsA[ row_idx ], row_idx, nCols, nCols );
    }
}

//////////////////////////////////////////////////////////////////////
// Copies rows 0, 1, 2, of A in to r0, r1, r2, of B
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterRows( const REAL *A, REAL *B,
        const IntArray &rowsB, int nCols )
{
    for ( int row_idx = 0; row_idx < rowsB.size(); row_idx++ )
    {
        MATRIX::copyRow( A, B, row_idx, rowsB[ row_idx ], nCols, nCols );
    }
}

//////////////////////////////////////////////////////////////////////
// Add one matrix to another (A = A + alpha * B
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy( REAL *A, const REAL *B, int rows, int cols,
        REAL alpha, int incA, int incB )
{
#ifdef SINGLE_PRECISION
    cblas_saxpy(rows * cols, alpha, B, incB, A, incA);
#else
    cblas_daxpy(rows * cols, alpha, B, incB, A, incA);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-matrix multiplication (C = beta * C + alpha * A * B)
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm( const REAL *A, const REAL *B, REAL *C,
        int rowsA, int colsA, int rowsB, int colsB,
        bool transposeA, bool transposeB,
        REAL alpha, REAL beta )
{
    int requiredCols = transposeB ? rowsB : colsB;

#ifdef SINGLE_PRECISION
    cblas_sgemm(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            transposeB ? CblasTrans : CblasNoTrans,
            transposeA ? colsA : rowsA, 
            transposeB ? rowsB : colsB,
            transposeA ? rowsA : colsA, 
            alpha,
            A, colsA,
            B, colsB,
            beta,
            C, requiredCols);
#else
    cblas_dgemm(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            transposeB ? CblasTrans : CblasNoTrans,
            transposeA ? colsA : rowsA, 
            transposeB ? rowsB : colsB,
            transposeA ? rowsA : colsA, 
            alpha,
            A, colsA,
            B, colsB,
            beta,
            C, requiredCols);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-matrix multiplication (C = beta * C + alpha * A * B)
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm( const REAL *A, const REAL *B, REAL *C,
        int rowsA, int colsA, int rowsB, int colsB,
        bool transposeA, bool transposeB,
        int lda_A, int lda_B, int lda_C,
        REAL alpha, REAL beta )
{
    int requiredCols = transposeB ? rowsB : colsB;

#ifdef SINGLE_PRECISION
    cblas_sgemm(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            transposeB ? CblasTrans : CblasNoTrans,
            transposeA ? colsA : rowsA, 
            transposeB ? rowsB : colsB,
            transposeA ? rowsA : colsA, 
            alpha,
            A, lda_A,
            B, lda_B,
            beta,
            C, lda_C);
#else
    cblas_dgemm(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            transposeB ? CblasTrans : CblasNoTrans,
            transposeA ? colsA : rowsA, 
            transposeB ? rowsB : colsB,
            transposeA ? rowsA : colsA, 
            alpha,
            A, lda_A,
            B, lda_B,
            beta,
            C, lda_C);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiplication (C := beta * C + alpha * A * b)
//////////////////////////////////////////////////////////////////////
void MATRIX::gemv( const REAL *A, const REAL *b, REAL *c,
        int rowsA, int colsA,
        bool transposeA,
        REAL alpha, REAL beta )
{
#ifdef SINGLE_PRECISION
    cblas_sgemv(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            rowsA,
            colsA,
            alpha,
            A, colsA,
            b, 1,
            beta,
            c, 1);
#else
    cblas_dgemv(
            CblasRowMajor,
            transposeA ? CblasTrans : CblasNoTrans,
            rowsA,
            colsA,
            alpha,
            A, colsA,
            b, 1,
            beta,
            c, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Symmetrix matrix-matrix update (C := beta * C + alpha * A * A')
// Updates the lower-triangular part of C
//
// If trans == true then C := beta * C + alpha * A' * A
//////////////////////////////////////////////////////////////////////
void MATRIX::syrk( const REAL *A, REAL *C,
        int rowsC, int k, /* A is either n x k or k x n */
        bool transpose, REAL alpha, REAL beta )
{
    int requiredCols = transpose ? rowsC : k;

#ifdef SINGLE_PRECISION
    cblas_ssyrk(
            CblasRowMajor,
            CblasLower, /* Work with lower triangular part of C */
            transpose ? CblasTrans: CblasNoTrans,
            rowsC, k,
            alpha,
            A, requiredCols,
            beta,
            C, rowsC);
#else
    cblas_dsyrk(
            CblasRowMajor,
            CblasLower, /* Work with lower triangular part of C */
            transpose ? CblasTrans: CblasNoTrans,
            rowsC, k,
            alpha,
            A, requiredCols,
            beta,
            C, rowsC);
#endif
}

//////////////////////////////////////////////////////////////////////
// Symmetrix matrix-matrix update (C := beta * C + alpha * A * A')
// Updates the lower-triangular part of C
//
// If trans == true then C := beta * C + alpha * A' * A
//////////////////////////////////////////////////////////////////////
void MATRIX::syrk( const REAL *A, REAL *C,
        int rowsC, int k, /* A is either n x k or k x n */
        int lda_A, int lda_C,
        bool transpose, REAL alpha, REAL beta )
{
    int requiredCols = transpose ? rowsC : k;

#ifdef SINGLE_PRECISION
    cblas_ssyrk(
            CblasRowMajor,
            CblasLower, /* Work with lower triangular part of C */
            transpose ? CblasTrans: CblasNoTrans,
            rowsC, k,
            alpha,
            A, lda_A,
            beta,
            C, lda_C);
#else
    cblas_dsyrk(
            CblasRowMajor,
            CblasLower, /* Work with lower triangular part of C */
            transpose ? CblasTrans: CblasNoTrans,
            rowsC, k,
            alpha,
            A, lda_A,
            beta,
            C, lda_C);
#endif
}

//////////////////////////////////////////////////////////////////////
// Get transpose (B = A^T)
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( REAL *A, const REAL *B, int rows, int cols )
{
    for (int y = 0; y < cols; y++)
        for (int x = 0; x < rows; x++)
            access( A, cols, rows, y, x ) = access( B, rows, cols, x, y );
}

//////////////////////////////////////////////////////////////////////
// Get transpose (B = A^T)
//////////////////////////////////////////////////////////////////////
void MATRIX::transposeBLAS( REAL *A, const REAL *B, int rows, int cols )
{
#ifdef SINGLE_PRECISION
    for ( int x = 0; x < rows; x++ )
    {
        cblas_scopy( cols, B + x * cols, 1, A + x, rows );
    }
#else
    for ( int x = 0; x < rows; x++ )
    {
        cblas_dcopy( cols, B + x * cols, 1, A + x, rows );
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// In place transpose for square matrices
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( REAL *A, int rows )
{
    for ( int i = 1; i < rows; i++ )
    {
        for ( int j = 0; j < i; j++ )
        {
            REAL temp = access( A, rows, rows, i, j );
            access( A, rows, rows, i, j ) = access( A, rows, rows, j, i );
            access( A, rows, rows, j, i ) = temp;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Compute eigenvalues/vectors.  We require everything here,
// including workspaces, etc.
// workspace should have size (7 * rows)
// vectorWorkspace should have size (rows * rows)
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem( REAL *A, int rows,
        REAL *eigenvalues, REAL *eigenvectors,
        REAL *workspace, REAL *vectorWorkspace )
{
    int rowsize = rows;
    int worksize = 5 * rows;

    REAL *work = workspace;
    REAL *valuesREAL = workspace + worksize;
    REAL *valuesImag = valuesREAL + rows;
    REAL *vectors = vectorWorkspace;

    // the actual LAPACK call
    int error;
#ifdef SINGLE_PRECISION
    sgeev("V","N", &rowsize, A, &rowsize, 
            valuesREAL, valuesImag, 
            vectors, &rowsize, NULL, &rowsize,
            work, &worksize, &error);
#else
    dgeev("V","N", &rowsize, A, &rowsize, 
            valuesREAL, valuesImag, 
            vectors, &rowsize, NULL, &rowsize,
            work, &worksize, &error);
#endif

    if (error != 0)
    {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " eigenvalue solver bombed!" << endl;
    }

    // copy out results
    for (int x = 0; x < rows; x++)
        eigenvalues[x] = valuesREAL[x];

    for (int x = 0; x < rows; x++)
        for (int y = 0; y < rows; y++)
            MATRIX::access( eigenvectors, rows, rows, x, y ) = vectors[x + y * rows];
}

//////////////////////////////////////////////////////////////////////
// Computes the Cholesky factor of the matrix stored in A.  A is
// overwritten
//////////////////////////////////////////////////////////////////////
int MATRIX::cholesky( REAL *A, int rows )
{
    int info = 0;

#if 0
#ifdef SINGLE_PRECISION
    spotrf( "U", &rows, A, &rows, &info );
#else
    dpotrf( "U", &rows, A, &rows, &info );
#endif
#endif
#ifdef USE_LAPACKE
#ifdef SINGLE_PRECISION
    info = LAPACKE_spotrf( LAPACK_ROW_MAJOR, 'L', rows, A, rows );
#else
    info = LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'L', rows, A, rows );
#endif
#else
#ifdef SINGLE_PRECISION
    spotrf( "U", &rows, A, &rows, &info );
#else
    dpotrf( "U", &rows, A, &rows, &info );
#endif
#endif

    return info;
}

//////////////////////////////////////////////////////////////////////
// Given a triangular matrix, and another input matrix, solve
// the associated triangular system.
//
// By default, solves L * X = B, where L is non-unit lower-triangular
//////////////////////////////////////////////////////////////////////
void MATRIX::triangularSolve( const REAL *L, REAL *B,
        int rowsB, int colsB,
        bool leftSide,
        bool lower,
        bool transpose,
        bool unitDiag,
        REAL alpha )
{
    CBLAS_SIDE side = leftSide ? CblasLeft : CblasRight;
    CBLAS_UPLO uplo = lower ? CblasLower : CblasUpper;
    CBLAS_TRANSPOSE transL = transpose ? CblasTrans : CblasNoTrans;
    CBLAS_DIAG diag = unitDiag ? CblasUnit : CblasNonUnit;

#ifdef SINGLE_PRECISION
    cblas_strsm(
            CblasRowMajor,
            side,
            uplo,
            transL,
            diag,
            rowsB, colsB,
            alpha,
            L,
            leftSide ? rowsB : colsB,
            B,
            colsB);
#else
    cblas_dtrsm(
            CblasRowMajor,
            side,
            uplo,
            transL,
            diag,
            rowsB, colsB,
            alpha,
            L,
            leftSide ? rowsB : colsB,
            B,
            colsB);
#endif
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but for a vector
//////////////////////////////////////////////////////////////////////
void MATRIX::triangularSolve( const REAL *L, REAL *b,
        int N,
        bool lower,
        bool transpose,
        bool unitDiag,
        REAL alpha )
{
    CBLAS_UPLO uplo = lower ? CblasLower : CblasUpper;
    CBLAS_TRANSPOSE transL = transpose ? CblasTrans : CblasNoTrans;
    CBLAS_DIAG diag = unitDiag ? CblasUnit : CblasNonUnit;

#ifdef SINGLE_PRECISION
    cblas_strsv(
            CblasRowMajor,
            uplo,
            transL,
            diag,
            N,
            L, N,
            b, 1);
#else
    cblas_dtrsv(
            CblasRowMajor,
            uplo,
            transL,
            diag,
            N,
            L, N,
            b, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Full rank least squares solve: A * x = b
//
// Overwrites A, and overwrites b with the solution
//////////////////////////////////////////////////////////////////////
int MATRIX::leastSquaresFullRank( REAL *A, REAL *b, int nRows, int nCols,
        int nRHS, bool transpose )
{
    int                        info;

#ifdef USE_LAPACKE
#ifdef SINGLE_PRECISION
    info = LAPACKE_sgels( LAPACK_ROW_MAJOR,
            transpose ? 'T' : 'N',
            nRows, nCols, nRHS,
            A, nCols,
            b, nRHS );
#else
    info = LAPACKE_dgels( LAPACK_ROW_MAJOR,
            transpose ? 'T' : 'N',
            nRows, nCols, nRHS,
            A, nCols,
            b, nRHS );
#endif
#else
    cerr << "NOT IMPLEMENTED!!!" << endl;
    abort();
#endif

    return info;
}

//////////////////////////////////////////////////////////////////////
// Least squares for an ill-conditioned matrix
//////////////////////////////////////////////////////////////////////
int MATRIX::leastSquaresSVD( REAL *A, REAL *b, REAL *singularValues,
        int nRows, int nCols,
        REAL tolerance, int &rank,
        int nRHS )
{
    int                        info;

#ifdef USE_LAPACKE
#ifdef SINGLE_PRECISION
    info = LAPACKE_sgelss( LAPACK_ROW_MAJOR,
            nRows, nCols, nRHS,
            A, nCols,
            b, nRHS,
            singularValues,
            tolerance, &rank );
#else
    info = LAPACKE_dgelss( LAPACK_ROW_MAJOR,
            nRows, nCols, nRHS,
            A, nCols,
            b, nRHS,
            singularValues,
            tolerance, &rank );
#endif
#else
    cerr << "NOT IMPLEMENTED!!!" << endl;
    abort();
#endif

    return info;
}

//////////////////////////////////////////////////////////////////////
// Symmetric positive definite solver
//////////////////////////////////////////////////////////////////////
int MATRIX::solveSPD( REAL *A, REAL *b, int nRows, int lda_A )
{
    int                        info;

#ifdef USE_LAPACKE
#ifdef SINGLE_PRECISION
    info = LAPACKE_sposv( LAPACK_ROW_MAJOR, 'L',
            nRows, 1,
            A, lda_A,
            b, 1 );
#else
    info = LAPACKE_dposv( LAPACK_ROW_MAJOR, 'L',
            nRows, 1,
            A, lda_A,
            b, 1 );
#endif
#else
    cerr << "NOT IMPLEMENTED!!!" << endl;
    abort();
#endif

    return info;
}

//////////////////////////////////////////////////////////////////////
// 2x2 eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem2x2( REAL *matrix, REAL *eigenvalues, REAL *eigenvectors )
{
    REAL A = access( matrix, 2, 2, 0, 0);
    REAL B = access( matrix, 2, 2, 0, 1);
    REAL C = access( matrix, 2, 2, 1, 0);
    REAL D = access( matrix, 2, 2, 1, 1);

    if (B * C <= 0.1e-20) {
        eigenvalues[0] = A; 
        access( eigenvectors, 2, 2, 0,0) = 1; 
        access( eigenvectors, 2, 2, 1,0) = 0;
        eigenvalues[1] = D; 
        access( eigenvectors, 2, 2, 0,1) = 0; 
        access( eigenvectors, 2, 2, 1,1) = 1;
        return;
    }

    REAL tr = A + D;
    REAL det = A * D - B * C;
    REAL S = sqrt(tr * tr * 0.25 - det );
    eigenvalues[0] = tr * 0.5 + S;
    eigenvalues[1] = tr * 0.5 - S;

    REAL temp = (A - D) * (A - D) * 0.25 + B * C;
    REAL SS = (temp < 0.0) ? 0.0 : sqrt(temp);
    if (A - D < 0.0) {
        access( eigenvectors, 2, 2, 0,0) = C;
        access( eigenvectors, 2, 2, 1,0) = -(A - D) * 0.5 + SS;
        access( eigenvectors, 2, 2, 0,1) = (A - D) * 0.5 - SS;
        access( eigenvectors, 2, 2, 1,1) = B;
    } 
    else {
        access( eigenvectors, 2, 2, 0,1) = C;
        access( eigenvectors, 2, 2, 1,1) = -(A - D) * 0.5 - SS;
        access( eigenvectors, 2, 2, 0,0) = (A - D) * 0.5 + SS;
        access( eigenvectors, 2, 2, 1,0) = B;
    }

    REAL n1 = sqrt(access( eigenvectors, 2, 2, 0,0) * access( eigenvectors, 2, 2, 0,0) +
            access( eigenvectors, 2, 2, 1,0) * access( eigenvectors, 2, 2, 1,0));
    REAL inv = 1.0 / n1;
    access( eigenvectors, 2, 2, 0,0) *= inv; 
    access( eigenvectors, 2, 2, 1,0) *= inv;
    REAL n2 = sqrt(access( eigenvectors, 2, 2, 0,1) * access( eigenvectors, 2, 2, 0,1) +
            access( eigenvectors, 2, 2, 1,1) * access( eigenvectors, 2, 2, 1,1));
    inv = 1.0 / n2;
    access( eigenvectors, 2, 2, 0,1) *= inv; 
    access( eigenvectors, 2, 2, 1,1) *= inv;
}

//////////////////////////////////////////////////////////////////////
// 3x3 low precision eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem3x3( REAL *A, REAL *eigenvalues, REAL *eigenvectors )
{
    register int i, j, k, l;

    //float TOL = 1e-8f;
    //int MAX_SWEEPS = 500;
    float TOL = 1e-3f;
    int MAX_SWEEPS = 5;
    unsigned int n = 3;

    float a[9];
    float d[3];
    float v[9];
    i = 0;
    for (int x = 0; x < 3; x++)
        for (int y = 0; y < 3; y++, i++)
        {
            a[i] = access( A, 3, 3, x, y);
            v[i] = (x == y) ? 1.0 : 0.0;
        }


    float onorm, dnorm;
    float b, dma, q, t, c, s;
    float atemp, vtemp, dtemp;

    // Set v to the identity matrix, set the vector d to contain the
    // diagonal elements of the matrix a
    d[0] = a[0];
    d[1] = a[4];
    d[2] = a[8];

    for (l = 1; l <= MAX_SWEEPS; l++)
    {
        // Set dnorm to be the maximum norm of the diagonal elements, set
        // onorm to the maximum norm of the off-diagonal elements

        dnorm = (float)fabs(d[0]) + (float)fabs(d[1]) + (float)fabs(d[2]);
        onorm = (float)fabs(a[1]) + (float)fabs(a[2]) + (float)fabs(a[5]);
        // Normal end point of this algorithm.
        if((onorm/dnorm) <= TOL)
            goto Exit_now;

        for (j = 1; j < static_cast<int>(n); j++)
        {
            for (i = 0; i <= j - 1; i++)
            {

                b = a[n*i+j];
                if(fabs(b) > 0.0f)
                {
                    dma = d[j] - d[i];
                    if((fabs(dma) + fabs(b)) <= fabs(dma))
                        t = b / dma;
                    else
                    {
                        q = 0.5f * dma / b;
                        t = 1.0f/((float)fabs(q) + (float)sqrt(1.0f+q*q));
                        if (q < 0.0)
                            t = -t;
                    }

                    c = 1.0f/(float)sqrt(t*t + 1.0f);
                    s = t * c;
                    a[n*i+j] = 0.0f;

                    for (k = 0; k <= i-1; k++)
                    {
                        atemp = c * a[n*k+i] - s * a[n*k+j];
                        a[n*k+j] = s * a[n*k+i] + c * a[n*k+j];
                        a[n*k+i] = atemp;
                    }

                    for (k = i+1; k <= j-1; k++)
                    {
                        atemp = c * a[n*i+k] - s * a[n*k+j];
                        a[n*k+j] = s * a[n*i+k] + c * a[n*k+j];
                        a[n*i+k] = atemp;
                    }

                    for (k = j+1; k < static_cast<int>(n); k++)
                    {
                        atemp = c * a[n*i+k] - s * a[n*j+k];
                        a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
                        a[n*i+k] = atemp;
                    }

                    for (k = 0; k < static_cast<int>(n); k++)
                    {
                        vtemp = c * v[n*k+i] - s * v[n*k+j];
                        v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
                        v[n*k+i] = vtemp;
                    }

                    dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
                    d[j] = s*s*d[i] + c*c*d[j] + 2.0f*c*s*b;
                    d[i] = dtemp;
                } /* end if */
            } /* end for i */
        } /* end for j */
    } /* end for l */

Exit_now:
    for (int x = 0; x < 3; x++)
        eigenvalues[x] = d[x];

    i = 0;
    for (int x = 0; x < 3; x++)
        for (int y = 0; y < 3; y++, i++)
            access( eigenvectors, 3, 3, x, y) = v[i];

    return;
}

//////////////////////////////////////////////////////////////////////
// 3x3 symmetric eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::symmetricEigensystem( const Matrix3<REAL> &M, Vector3d &eigenvalues,
                                   Matrix3d &eigenvectors )
{
    // Storage for the original matrix
    REAL                       A[ 9 ];

    // Storage for the output matrix
    REAL                       Z[ 9 ];
    REAL                       evOut[ 3 ];

    int                        info;
    int                        numEigenvalues;

    int                        isuppz[ 6 ];

    // Copy the input matrix in row-major format
    A[ 0 ] = M.cols[ 0 ][ 0 ];
    A[ 1 ] = M.cols[ 1 ][ 0 ];
    A[ 2 ] = M.cols[ 2 ][ 0 ];
    A[ 3 ] = M.cols[ 0 ][ 1 ];
    A[ 4 ] = M.cols[ 1 ][ 1 ];
    A[ 5 ] = M.cols[ 2 ][ 1 ];
    A[ 6 ] = M.cols[ 0 ][ 2 ];
    A[ 7 ] = M.cols[ 1 ][ 2 ];
    A[ 8 ] = M.cols[ 2 ][ 2 ];

#ifdef SINGLE_PRECISION
    info = LAPACKE_ssyevr( LAPACK_ROW_MAJOR,
                           'V', /* evals and evecs */
                           'A', /* all eigenvalues */
                           'L', /* lower triangular part is stored */
                           3, /* Matrix dimension */
                           A, 3, /* System and leading dimension */
                           0.0, 0.0, 0, 0, /* Unused parameters */
                           1e-14, /* Absolute tolerance */
                           &numEigenvalues,
                           evOut,
                           Z, 3, /* eigenvectors and leading dimension */
                           isuppz );
#else
    info = LAPACKE_dsyevr( LAPACK_ROW_MAJOR,
                         'V', /* evals and evecs */
                         'A', /* all eigenvalues */
                         'L', /* lower triangular part is stored */
                         3, /* Matrix dimension */
                         A, 3, /* System and leading dimension */
                         0.0, 0.0, 0, 0, /* Unused parameters */
                         1e-14, /* Absolute tolerance */
                         &numEigenvalues,
                         evOut,
                         Z, 3, /* eigenvectors and leading dimension */
                         isuppz );
#endif

    TRACE_ASSERT( info == 0, "Eigensolver failed" );

    // Copy the output
    eigenvectors.cols[ 0 ][ 0 ] = Z[ 0 ];
    eigenvectors.cols[ 1 ][ 0 ] = Z[ 1 ];
    eigenvectors.cols[ 2 ][ 0 ] = Z[ 2 ];
    eigenvectors.cols[ 0 ][ 1 ] = Z[ 3 ];
    eigenvectors.cols[ 1 ][ 1 ] = Z[ 4 ];
    eigenvectors.cols[ 2 ][ 1 ] = Z[ 5 ];
    eigenvectors.cols[ 0 ][ 2 ] = Z[ 6 ];
    eigenvectors.cols[ 1 ][ 2 ] = Z[ 7 ];
    eigenvectors.cols[ 2 ][ 2 ] = Z[ 8 ];

    eigenvalues[ 0 ] = evOut[ 0 ];
    eigenvalues[ 1 ] = evOut[ 1 ];
    eigenvalues[ 2 ] = evOut[ 2 ];
}

//////////////////////////////////////////////////////////////////////
// Scales a matrix.  A *= alpha
//////////////////////////////////////////////////////////////////////
void MATRIX::scale( REAL *A, int rows, int cols, REAL alpha )
{
#ifdef SINGLE_PRECISION
    cblas_sscal(cols * rows, alpha, A, 1);
#else
    cblas_dscal(cols * rows, alpha, A, 1);
#endif
}
