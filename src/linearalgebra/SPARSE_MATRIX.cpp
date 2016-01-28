// SPARSE_MATRIX.h: interface for the SPARSE_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include <map>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "SPARSE_MATRIX.h"

#include <utils/IO.h>
#include <utils/trace.h>

#include <vector>
#include <set>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX::SPARSE_MATRIX(int rows, int cols) :
    _rows(rows), _cols(cols), _dud(0.0f)
{
}

SPARSE_MATRIX::SPARSE_MATRIX() :
    _dud(0.0f)
{
}

SPARSE_MATRIX::SPARSE_MATRIX(MATRIX& matrix) :
    _dud(0.0f)
{
    _rows = matrix.rows();
    _cols = matrix.cols();

    for (int y = 0; y < _rows; y++)
        for (int x = 0; x < _cols; x++)
            (*this)(y,x) = matrix(y,x);
}

//////////////////////////////////////////////////////////////////////
// check if an entry already exists
//////////////////////////////////////////////////////////////////////
bool SPARSE_MATRIX::exists(int row, int col)
{
    pair<int,int> index(row, col);
    map<pair<int,int>, REAL>::iterator i = _matrix.find(index);
    if (i != _matrix.end())
        return true;
    return false;
}

//////////////////////////////////////////////////////////////////////
// return a reference to an entry
//////////////////////////////////////////////////////////////////////
REAL& SPARSE_MATRIX::operator()(int row, int col)
{
    /*
    // bounds check
    if (col < 0 || row < 0 || col > _cols - 1 || row > _rows - 1)
    {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " index out of bounds! " << endl;
    return _dud;
    }
     */

    TRACE_ASSERT( row < rows() && col < cols() );

    // lookup the entry
    pair<int,int> index(row, col);
    map<pair<int,int>, REAL>::iterator i = _matrix.find(index);
    if (i != _matrix.end())
        return i->second;

    // if it doesn't exist, create it
    _matrix[index] = 0.0f;
    return _matrix[index];
}

REAL SPARSE_MATRIX::operator() (int row, int col) const
{
    /*
    // bounds check
    if (col < 0 || row < 0 || col > _cols - 1 || row > _rows - 1)
    {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " index out of bounds! " << endl;
    return _dud;
    }
     */

    TRACE_ASSERT( row < rows() && col < cols() );

    // lookup the entry
    pair<int,int> index(row, col);
    map<pair<int,int>, REAL>::const_iterator i = _matrix.find(index);
    if (i != _matrix.end())
        return i->second;

    return 0.0;
}

//////////////////////////////////////////////////////////////////////
// dump all the matrix entries to a Matlab file
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::writeToMatlab(string filename, string varName)
{
    FILE* file;
    file = fopen(filename.c_str(), "w");

    // create the matrix and sparsify it
    fprintf(file, "%s = [];\n", varName.c_str());
    fprintf(file, "%s = sparse(%s);\n", varName.c_str(), varName.c_str());

    // iterate through all the entries
    map<pair<int,int>, REAL>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
    {
        pair<int,int> index = i->first;
        REAL value = i->second;
        fprintf(file, "%s(%i,%i) = %.16f;\n", varName.c_str(),
                index.first + 1, index.second + 1, value);
    }

    fclose(file);
}

void SPARSE_MATRIX::writeToBinary(string filename)
{
    ofstream f( filename.c_str(), ios::binary );
    if( f.fail() )
    {
        cerr << "Couldn't write sparse matrix to binary file: " << filename.c_str() << endl;
    }
    else
    {
        // (stevenan) 2008-11-24 13:28:02
        // We write a -1 as a header, to specify that this is a
        // struct-of-arrays-format binary sparse matrix.
        // This is what our matlab read function expects.
        // This is to maintain backwards compatibility with the old
        // array-of-structs format.
        int header = -1;
        f.write( (char*)&header, sizeof(int) );

        f.write( (char*)&_rows, sizeof(int) );
        f.write( (char*)&_cols, sizeof(int) );
        int nnz = _matrix.size();
        f.write( (char*)&nnz, sizeof(int) );

        // Write row indices
        map<pair<int,int>, REAL>::iterator i;
        for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
            pair<int,int> index = i->first;
            // This does NOT add 1 to the index. The matlab read function will do
            // that on its own.
            f.write( (char*)&index.first, sizeof(int) );
        }
        // Write columns
        for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
            pair<int,int> index = i->first;
            // This does NOT add 1 to the index. The matlab read function will do
            // that on its own.
            f.write( (char*)&index.second, sizeof(int) );
        }
        // Write values
        for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
            //pair<int,int> index = i->first;
            double valueBin = (double)i->second;
            f.write( (char*)&valueBin, sizeof(double) );
        }
    }
    f.close();
}


//////////////////////////////////////////////////////////////////////
// count the number of non-zeros per row and return them
//////////////////////////////////////////////////////////////////////
vector<int> SPARSE_MATRIX::nonZerosPerRow()
{
    vector<int> count(_rows);

    // iterate through all the entries
    map<pair<int,int>, REAL>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
        count[i->first.first]++;

    return count;
}

//////////////////////////////////////////////////////////////////////
// Print a specific row
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::printRow(int row)
{
    if (row < 0) return;
    if (row > _rows) return;

    for (int x = 0; x < _cols; x++)
        if (exists(x, row))
            cout << "A(" << row << "," << x << "): " << (*this)(x, row) << endl;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const SPARSE_MATRIX& matrix)
{
    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator i;
    const map<pair<int,int>, REAL>& data = matrix.matrix();
    for (i = data.begin(); i != data.end(); i++)
    {
        const pair<int,int> index = i->first;
        const REAL value = i->second;
        out << "K(" << index.first+1 << "," << index.second+1 << ") = " <<  value  << ";" << endl;
    }
    return out;
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//
// Note that axpy actually applies to vectors, but in this
// case we can just treat the matrix as a vector and multiply
// all its elements by alpha
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::axpy(REAL alpha, SPARSE_MATRIX& A)
{
    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator i;
    const map<pair<int,int>, REAL>& data = A.matrix();
    for (i = data.begin(); i != data.end(); i++)
    {
        const pair<int,int> index = i->first;
        const REAL value = i->second;
        (*this)(index.first, index.second) += alpha * value;
    }
}

//////////////////////////////////////////////////////////////////////
// gemv operation
//
// y = alpha * A * x
// (Don't let y == x)
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::gemv(REAL alpha, VECTOR &x, VECTOR &y)
{
    assert(this->cols() == x.size());

    if ( y.size() != this->rows() )
    {
        y.resizeAndWipe( this->rows() );
    }
    y.clear();

    // Iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator iter;
    const map<pair<int,int>, REAL>& data = this->matrix();
    for ( iter = data.begin(); iter != data.end(); iter++ )
    {
        const pair<int,int> index = iter->first;
        const REAL value = iter->second;
        int i = index.first;
        int j = index.second;
        y(i) += x(j) * value;
    }

    y *= alpha;
}

//////////////////////////////////////////////////////////////////////
// gemv operation
//
// y = A * x
// (Don't let y == x)
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::gemv(VECTOR &x, VECTOR &y)
{
    gemv(1.0, x, y);
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const SPARSE_MATRIX& A, const VECTOR& x) 
{
    assert(A.cols() == x.size());

    VECTOR y(A.rows());

    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator iter;
    const map<pair<int,int>, REAL>& data = A.matrix();
    for (iter = data.begin(); iter != data.end(); iter++)
    {
        const pair<int,int> index = iter->first;
        const REAL value = iter->second;
        int i = index.first;
        int j = index.second;
        TRACE_ASSERT( i < A.rows() && j < A.cols() );
        y(i) += x(j) * value;
    }
    return y;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const SPARSE_MATRIX::SparseColumnMatrix &A, const VECTOR &x)
{
    assert((int)A._ncol == (int)x.size());

    VECTOR y(A._nrow);

    for ( size_t col_idx = 0; col_idx < A._ncol; col_idx++ )
    {
        for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
                row_ptr++ )
        {
            y( A._i[ row_ptr ] ) += A._x[ row_ptr ] * x( col_idx );
        }
    }

    return y;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(const SPARSE_MATRIX& A, const REAL& alpha) 
{
    SPARSE_MATRIX B(A);
    B *= alpha;
    return B;
}

//////////////////////////////////////////////////////////////////////
// sparse-full matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const SPARSE_MATRIX& A, const MATRIX& B)
{
    MATRIX C(A.rows(), B.cols());

    const map<pair<int,int>, REAL>& matrix = A.matrix(); 

    map<pair<int,int>, REAL>::const_iterator iter;
    for (iter = matrix.begin(); iter != matrix.end(); iter++)
    {
        pair<int,int> index = iter->first;
        int row = index.first;
        int col = index.second;
        REAL entry = iter->second;

        for (int x = 0; x < B.cols(); x++)
            C(row, x) += B(col, x) * entry;
    }

    return C;
} 

//////////////////////////////////////////////////////////////////////
// scale matrix by a scalar
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator*=(const REAL& alpha) 
{
    // iterate through all the entries
    map<pair<int,int>, REAL>::iterator iter;
    for (iter = _matrix.begin(); iter != _matrix.end(); iter++)
        iter->second *= alpha;

    return *this;
}

//////////////////////////////////////////////////////////////////////
// add two sparse matrices together
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator+=(const SPARSE_MATRIX& A) 
{
    assert(A.rows() == _rows && A.cols() == _cols);

    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator iter;
    const map<pair<int,int>, REAL>& data = A.matrix();
    for (iter = data.begin(); iter != data.end(); iter++)
    {
        const pair<int,int> index = iter->first;
        const REAL value = iter->second;
        int i = index.first;
        int j = index.second;
        (*this)(i,j) += value;
    }

    return *this;
}

//////////////////////////////////////////////////////////////////////
// subtract one sparse matrix from another
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator-=(const SPARSE_MATRIX& A) 
{
    assert(A.rows() == _rows && A.cols() == _cols);

    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator iter;
    const map<pair<int,int>, REAL>& data = A.matrix();
    for (iter = data.begin(); iter != data.end(); iter++)
    {
        const pair<int,int> index = iter->first;
        const REAL value = iter->second;
        int i = index.first;
        int j = index.second;
        (*this)(i,j) -= value;
    }

    return *this;
}

//////////////////////////////////////////////////////////////////////
// copy A into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::copies(SPARSE_MATRIX& A)
{
    assert(A.cols() == _cols && A.rows() == _rows);

    if (A._cols == _cols && A._rows == _rows)
    {
        // first zero out what's already there, but don't do
        // a total clear, or else it will have to reallocate all the
        // entries
        map<pair<int,int>, REAL>::iterator i;
        for (i = _matrix.begin(); i != _matrix.end(); i++)
            i->second = 0.0;

        // next copy in all the entries from A. If an entry already
        // existed at the matrix entry, we will save the memory
        // allocation
        for (i = A._matrix.begin(); i != A._matrix.end(); i++)
        {
            const pair<int,int> index = i->first;
            (*this)(index.first, index.second) = i->second;
        }
        return;
    }

    // else do an expensive deep copy
    _matrix = map<pair<int,int>, REAL>(A.matrix());
}

//////////////////////////////////////////////////////////////////////
// copy A into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subcopies(MATRIX& A, int row, int col)
{
    for (int y = 0; y < A.rows(); y++)
        for (int x = 0; x < A.cols(); x++)
            (*this)(y + row, x + col) = A(y,x);
}

//////////////////////////////////////////////////////////////////////
// add A into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::add(MATRIX& A, int row, int col)
{
    for (int y = 0; y < A.rows(); y++)
        for (int x = 0; x < A.cols(); x++)
            (*this)(y + row, x + col) += A(y,x);
}

//////////////////////////////////////////////////////////////////////
// subtract A from the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subtract(MATRIX& A, int row, int col)
{
    for (int y = 0; y < A.rows(); y++)
        for (int x = 0; x < A.cols(); x++)
            (*this)(y + row, x + col) -= A(y,x);
}

//////////////////////////////////////////////////////////////////////
// All the matrix entries repackaged into vectors,
// ie _matrix(rows[x], cols[x]) = vals[x];
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::entries(vector<int>& rows, vector<int>& cols, vector<REAL>& values)
{
    // wipe previous entries
    rows.clear();
    cols.clear();
    values.clear();

    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
    {
        const pair<int,int> index = i->first;
        const REAL value = i->second;
        rows.push_back(index.first);
        cols.push_back(index.second);
        values.push_back(value);
    }
}

//////////////////////////////////////////////////////////////////////
// All the matrix entries repackaged into vectors,
// ie _matrix(rows[x], cols[x]) = vals[x];
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::entries(vector<int>& rows, vector<int>& cols, vector<const REAL*>& values)
{
    // wipe previous entries
    rows.clear();
    cols.clear();
    values.clear();

    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
    {
        const pair<int,int> index = i->first;
        //REAL value = i->second;
        rows.push_back(index.first);
        cols.push_back(index.second);
        values.push_back(&(i->second));
    }
}

REAL SPARSE_MATRIX::sum() const
{
    REAL total = 0.0;
    map<pair<int,int>, REAL>::const_iterator i;
    // Write values
    for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
        //pair<int,int> index = i->first;
        total += (REAL)i->second;
    }
    return total;
}

//////////////////////////////////////////////////////////////////////
// Set the matrix to zero. Note this will *NOT* stomp the underlying map!
// It will instead set all current entries to zero so that we are not
// forced to reallocate the sparsity structure again.
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::clear()
{
    map<pair<int,int>, REAL>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
        i->second = 0.0;
}

//////////////////////////////////////////////////////////////////////
// Clear which stomps the map structure as well
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::clearFull()
{
    _matrix.clear();
}

//////////////////////////////////////////////////////////////////////
// Convert to a full matrix
//////////////////////////////////////////////////////////////////////
MATRIX SPARSE_MATRIX::full()
{
    MATRIX matrix(_rows, _cols);

    // iterate through all the entries
    map<pair<int,int>, REAL>::const_iterator i;
    const map<pair<int,int>, REAL>& data = this->matrix();
    for (i = data.begin(); i != data.end(); i++)
    {
        const pair<int,int> index = i->first;
        const REAL value = i->second;
        matrix(index.first, index.second) = value;
    }
    return matrix;
}

//////////////////////////////////////////////////////////////////////
// Erase a row/column
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::eraseRowColumn(int rowCol)
{
    pair<int,int> entry(rowCol, rowCol);

    map<pair<int,int>, REAL>::iterator i = _matrix.find(entry);

    if (i == _matrix.end())
    {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Row/Col erasure failed!" << endl;
        return;
    }
    while (i->first.first == rowCol)
    {
        i->second = 0.0;

        // zero out the equivalent column entry
        pair<int,int> transpose(i->first.second, i->first.first);
        _matrix[transpose] = 0.0;

        i++;
    }

    map<pair<int,int>, REAL>::iterator ri = _matrix.find(entry);
    while (ri->first.first == rowCol)
    {
        ri->second = 0.0;

        // zero out the equivalent column entry
        pair<int,int> transpose(ri->first.second, ri->first.first);
        _matrix[transpose] = 0.0;
        ri--;
    }
}

//////////////////////////////////////////////////////////////////////
// 1 norm of the matrix
// - maximum absolute column sum
//////////////////////////////////////////////////////////////////////
REAL SPARSE_MATRIX::norm1()
{
    VECTOR colSums(_cols);
    map<pair<int,int>, REAL>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
    {
        pair<int,int> index = i->first;
        REAL entry = i->second;
        colSums[index.second] += fabs(entry);
    }
    return colSums.maxValue();
}

//////////////////////////////////////////////////////////////////////
// Inf norm of the matrix
// - maximum absolute row sum
//////////////////////////////////////////////////////////////////////
REAL SPARSE_MATRIX::normInf()
{
    VECTOR rowSums(_rows);
    map<pair<int,int>, REAL>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
    {
        pair<int,int> index = i->first;
        REAL entry = i->second;
        rowSums[index.first] += fabs(entry);
    }
    return rowSums.maxValue();
}

//////////////////////////////////////////////////////////////////////
// Frobenius norm of the matrix
//////////////////////////////////////////////////////////////////////
REAL SPARSE_MATRIX::normFrob()
{
    REAL total = 0.0;
    map<pair<int,int>, REAL>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
    {
        total += i->second * i->second;
    }
    return sqrt( total );
}

//////////////////////////////////////////////////////////////////////
// Sets out = this * A
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::denseMultiply( MATRIX &A, MATRIX &out )
{
    //const int n = A.rows();
    const int m = A.cols();

    TRACE_ASSERT( _cols == n, "Matrix dimensions do not match (%d, %d)", _cols, n );

    out.resizeAndWipe( _rows, m );

    // Iterate over all entries and add to the output matrix where necessary
    map<pair<int,int>, REAL>::const_iterator entryIt;
    for( entryIt = _matrix.begin(); entryIt != _matrix.end(); ++entryIt )
    {
        pair<int, int> index = entryIt->first;
        REAL value = entryIt->second;

        int row = index.first;
        int col = index.second;

        // Add to each entry in this row of the result matrix
        for ( int i = 0; i < m; i++ )
        {
            out( row, i ) += value * A( col, i );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Converts this to a dense matrix
//////////////////////////////////////////////////////////////////////
MATRIX SPARSE_MATRIX::dense()
{
    MATRIX denseVersion( _rows, _cols );

    map<pair<int,int>, REAL>::const_iterator entryIt;
    for ( entryIt = _matrix.begin(); entryIt != _matrix.end(); ++entryIt )
    {
        pair<int, int> index = entryIt->first;
        REAL value = entryIt->second;

        int row = index.first;
        int col = index.second;

        denseVersion( row, col ) = value;
    }

    return denseVersion;
}

//////////////////////////////////////////////////////////////////////
// Copies the matrix in to sparse column form
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::constructSparseColumnCopy( SparseColumnMatrix &M ) const
{
    M.clear();

    M._nrow = this->rows();
    M._ncol = this->cols();
    M._nzmax = this->size();

    M._p = (int *)malloc( (M._ncol + 1) * sizeof( int ) );
    M._i = (int *)malloc( M._nzmax * sizeof( int ) );
    M._x = (REAL *)malloc( M._nzmax * sizeof( REAL ) );

    memset( (void *)M._p, 0, (M._ncol + 1) * sizeof( int ) );

    // Start by getting counts for all columns
    for ( SPARSE_MATRIX::const_iterator i = this->begin();
            i != this->end(); i++ )
    {
        int    col_idx = i->first.second;

        M._p[ col_idx + 1 ] += 1;
    }

    // Construct the final p
    for ( size_t i = 0; i < M._ncol; i++ )
    {
        M._p[ i + 1 ] += M._p[ i ];
    }

    int *colcounts = (int *)malloc( M._ncol * sizeof( int ) );

    memset( (void *)colcounts, 0, M._ncol * sizeof( int ) );

    // Figure out where to put row pointers
    for ( const_iterator i = begin(); i != end(); i++ )
    {
        int    row_idx = i->first.first;
        int    col_idx = i->first.second;
        REAL   value = i->second;

        M._i[ M._p[ col_idx ] + colcounts[ col_idx ] ] = row_idx;
        M._x[ M._p[ col_idx ] + colcounts[ col_idx ] ] = value;

        colcounts[ col_idx ] += 1;
    }

    free( colcounts );
}

//////////////////////////////////////////////////////////////////////
// Multiplies the matrix with the given vector
//////////////////////////////////////////////////////////////////////
int SPARSE_MATRIX::matrixMultiply( const SparseColumnMatrix &A,
                                   const REAL *g, REAL *b,
                                   REAL alpha,
                                   bool clear, bool transpose,
                                   int ldaG, int ldaB )
{
    int                row;

    ldaG = max( 1, ldaG );
    ldaB = max( 1, ldaB );

    TRACE_ASSERT( ldaB == 1, "Cannot alter leading dimension of b" );

    if ( clear )
    {
        if ( transpose )
        {
            memset( (void *)b, 0, A._ncol * 1 * sizeof( REAL ) );
        }
        else
        {
            memset( (void *)b, 0, A._nrow * 1 * sizeof( REAL ) );
        }
    }

    for ( size_t col_idx = 0; col_idx < A._ncol; col_idx++ ) {
        for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
                row_idx++ )
        {
            row = A._i[ row_idx ];

            if ( transpose )
            {
                b[ col_idx ]
                    += alpha * A._x[ row_idx ] * g[ ( row ) * ldaG ];
#if 0
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( col_idx - columnStart ) * ldaB,
                        /* Align with correct row */
                        G + ( row - rowStart ) * ldaG,
                        /* Align with correct row */
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
#endif
            }
            else
            {
                b[ row ]
                    += alpha * A._x[ row_idx ] * g[ ( col_idx ) * ldaG ];
#if 0
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( row - rowStart ) * ldaB /* Align with correct row */,
                        G + ( col_idx - columnStart ) * ldaG /* "             " */,
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
#endif
            }
        }
    }

    return 0;
}

//////////////////////////////////////////////////////////////////////
// Multiplies a sub matrix of the given sparse column matrix
// with the given matrix G, placing the result in B
//////////////////////////////////////////////////////////////////////
int SPARSE_MATRIX::subMatrixMultiply( const SparseColumnMatrix &A,
                                      const REAL *G, REAL *B,
                                      int rowStart, int columnStart,
                                      int nRows, int nCols, int nColsG,
                                      REAL alpha,
                                      bool clear, bool transpose,
                                      int ldaG, int ldaB )
{
    int                row;
    int                nnz = 0;

    ldaG = max( nColsG, ldaG );
    ldaB = max( nColsG, ldaB );

    TRACE_ASSERT( ldaB == nColsG, "Cannot alter leading dimension of B" );

    TRACE_ASSERT( rowStart >= 0 && rowStart + nRows <= A._nrow,
            "Row range out of bounds" );
    TRACE_ASSERT( columnStart >= 0 && columnStart + nCols <= A._ncol,
            "Column range out of bounds" );

    if ( clear )
    {
        if ( transpose )
        {
            memset( (void *)B, 0, nCols * nColsG * sizeof( REAL ) );
        }
        else
        {
            memset( (void *)B, 0, nRows * nColsG * sizeof( REAL ) );
        }
    }

    for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ )
        for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
                row_idx++ )
        {
            row = A._i[ row_idx ];

            if ( row < rowStart || row >= rowStart + nRows )
            {
                continue;
            }

            nnz += 1;

            if ( transpose )
            {
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( col_idx - columnStart ) * ldaB,
                        /* Align with correct row */
                        G + ( row - rowStart ) * ldaG,
                        /* Align with correct row */
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
            }
            else
            {
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( row - rowStart ) * ldaB /* Align with correct row */,
                        G + ( col_idx - columnStart ) * ldaG /* "             " */,
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
            }
        }

    return nnz;
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but with vector inputs g and b
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subMatrixMultiply( const SparseColumnMatrix &A,
                                       const REAL *g, REAL *b,
                                       int rowStart, int columnStart,
                                       int nRows, int nCols,
                                       REAL alpha,
                                       bool clear, bool transpose,
                                       int ldaG, int ldaB )
{
    int                row;
    int                rowMax = rowStart + nRows;

    ldaG = max( 1, ldaG );
    ldaB = max( 1, ldaB );

    TRACE_ASSERT( ldaB == 1, "Cannot alter leading dimension of b" );

    TRACE_ASSERT( rowStart >= 0 && rowStart + nRows <= A._nrow,
            "Row range out of bounds" );
    TRACE_ASSERT( columnStart >= 0 && columnStart + nCols <= A._ncol,
            "Column range out of bounds" );

    if ( clear ) {
        if ( transpose ) {
            memset( (void *)b, 0, nCols * 1 * sizeof( REAL ) );
        } else {
            memset( (void *)b, 0, nRows * 1 * sizeof( REAL ) );
        }
    }

    for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ ) {
        for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
                row_idx++ )
        {
            row = A._i[ row_idx ];

#if 0
            if ( row < rowStart || row >= rowStart + nRows )
            {
                continue;
            }
#endif
            // Assumes rows are in sorted order
            if ( row < rowStart ) {
                continue;
            } else if ( row >= rowMax ) {
                break;
            }

            if ( transpose )
            {
                b[ col_idx - columnStart ]
                    += alpha * A._x[ row_idx ] * g[ ( row - rowStart ) * ldaG ];
#if 0
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( col_idx - columnStart ) * ldaB,
                        /* Align with correct row */
                        G + ( row - rowStart ) * ldaG,
                        /* Align with correct row */
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
#endif
            }
            else
            {
                b[ row - rowStart ]
                    += alpha * A._x[ row_idx ] * g[ ( col_idx - columnStart ) * ldaG ];
#if 0
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( row - rowStart ) * ldaB /* Align with correct row */,
                        G + ( col_idx - columnStart ) * ldaG /* "             " */,
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
#endif
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Version in which all rows after the given start rows in the
// provided column set are used
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subMatrixMultiply( const SparseColumnMatrix &A,
                                       const REAL *g, REAL *b,
                                       int rowStart, int columnStart,
                                       const IntArray &startRows,
                                       REAL alpha,
                                       bool clear, bool transpose,
                                       int ldaG )
{
    int                row;
    int                nCols = startRows.size();
    int                nRows = A._nrow - rowStart;

    ldaG = max( 1, ldaG );

    TRACE_ASSERT( rowStart >= 0 && rowStart + nRows <= A._nrow,
            "Row range out of bounds" );
    TRACE_ASSERT( columnStart >= 0 && columnStart + nCols; <= A._ncol,
            "Column range out of bounds" );

    if ( clear ) {
        if ( transpose ) {
            memset( (void *)b, 0, nCols * 1 * sizeof( REAL ) );
        } else {
            memset( (void *)b, 0, nRows * 1 * sizeof( REAL ) );
        }
    }

    if ( transpose ) {
        for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ ) {
            for ( int row_idx = A._p[ col_idx ] + startRows[ col_idx - columnStart ];
                    row_idx < A._p[ col_idx + 1 ]; row_idx++ )
            {
                row = A._i[ row_idx ];

                b[ col_idx - columnStart ]
                    += alpha * A._x[ row_idx ] * g[ ( row - rowStart ) * ldaG ];
            }
        }
    } else {
        for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ ) {
            for ( int row_idx = A._p[ col_idx ] + startRows[ col_idx - columnStart ];
                    row_idx < A._p[ col_idx + 1 ]; row_idx++ )
            {
                row = A._i[ row_idx ];

                b[ row - rowStart ]
                    += alpha * A._x[ row_idx ] * g[ ( col_idx - columnStart ) * ldaG ];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Version of the above with an extra workspace for faster (hopefully)
// vectorized operations
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subMatrixMultiply( const SparseColumnMatrix &A,
                                       const REAL *g, REAL *b,
                                       int rowStart, int columnStart,
                                       const IntArray &startRows,
                                       REAL alpha,
                                       REAL *workspace,
                                       bool clear, bool transpose,
                                       int ldaG )
{
    int                row;
    int                nCols = startRows.size();
    int                nRows = A._nrow - rowStart;

    ldaG = max( 1, ldaG );

    TRACE_ASSERT( rowStart >= 0 && rowStart + nRows <= A._nrow,
            "Row range out of bounds" );
    TRACE_ASSERT( columnStart >= 0 && columnStart + nCols; <= A._ncol,
            "Column range out of bounds" );

    if ( clear ) {
        if ( transpose ) {
            memset( (void *)b, 0, nCols * 1 * sizeof( REAL ) );
        } else {
            memset( (void *)b, 0, nRows * 1 * sizeof( REAL ) );
        }
    }

    if ( transpose ) {
        for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ ) {
            // Collect required entries for this column
            int workIdx = 0;
            for ( int row_idx = A._p[ col_idx ] + startRows[ col_idx - columnStart ];
                    row_idx < A._p[ col_idx + 1 ]; row_idx++ )
            {
                row = A._i[ row_idx ];
                workspace[ workIdx ] = g[ ( row - rowStart ) * ldaG ];
                workIdx++;
            }

            // Dot product with the matrix entries from this column
            REAL *columnData = &A._x[ A._p[ col_idx ]
                + startRows[ col_idx - columnStart ] ];
            b[ col_idx - columnStart ]
                += alpha * VECTOR::dot( workIdx, columnData, workspace );
        }
    } else {
        for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ ) {
            int columnRows = A._p[ col_idx + 1 ]
                - A._p[ col_idx ]
                - startRows[ col_idx - columnStart ];
            REAL *columnData = &A._x[ A._p[ col_idx ]
                + startRows[ col_idx - columnStart ] ];

            MATRIX::copy( workspace, columnData, columnRows, 1 );
            MATRIX::scale( workspace, columnRows, 1,
                    alpha * g[ (col_idx - columnStart) * ldaG ] );

            // Copy back to the relevant rows
            int workIdx = 0;
            for ( int row_idx = A._p[ col_idx ] + startRows[ col_idx - columnStart ];
                    row_idx < A._p[ col_idx + 1 ]; row_idx++ )
            {
                row = A._i[ row_idx ];
                b[ row - rowStart ] += workspace[ workIdx ];
                workIdx++;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but uses a specific row set
//////////////////////////////////////////////////////////////////////
int SPARSE_MATRIX::subMatrixMultiply( const SparseColumnMatrix &A,
                                      const REAL *G, REAL *B,
                                      int rowStart, int columnStart,
                                      int nRows, int nCols, int nColsG,
                                      const IntArray &inverseRowList,
                                      REAL alpha,
                                      bool clear, bool transpose,
                                      int ldaG, int ldaB )
{
    int                row;
    int                nnz = 0;
    int                nRowsFull = inverseRowList.size();

    ldaG = max( nColsG, ldaG );
    ldaB = max( nColsG, ldaB );

    TRACE_ASSERT( ldaB == nColsG, "Cannot alter leading dimension of B" );

    TRACE_ASSERT( rowStart >= 0 && rowStart + nRows <= A._nrow,
            "Row range out of bounds" );
    TRACE_ASSERT( nRows <= inverseRowList.size() );
    TRACE_ASSERT( columnStart >= 0 && columnStart + nCols <= A._ncol,
            "Column range out of bounds" );

    if ( clear ) {
        if ( transpose ) {
            memset( (void *)B, 0, nCols * nColsG * sizeof( REAL ) );
        } else {
            memset( (void *)B, 0, nRows * nColsG * sizeof( REAL ) );
        }
    }

    for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ )
        for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
                row_idx++ )
        {
            row = A._i[ row_idx ];

            if ( row < rowStart || row >= rowStart + nRowsFull ) {
                continue;
            }

            // Figure out which row this is supposed to map to in the output (if
            // not transposing) or input (if transposing)
            row -= rowStart;
            row = inverseRowList[ row ];

            // FIXME: debugging
            if ( row < 0 || row >= nRows ) {
                // FIXME: There is something about this that I don't like
#if 0
                cout << SDUMP( row ) << endl;
                TRACE_ASSERT( NULL, "Bad row index" );
#endif
                continue;
            }
            TRACE_ASSERT( row >= 0 && row < nRows,
                    "Inverse row lookup error" );

            nnz += 1;

            if ( transpose ) {
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( col_idx - columnStart ) * ldaB,
                        /* Align with correct row */
                        G + ( row ) * ldaG,
                        /* Align with correct row */
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
            } else {
                // Add row columnIndex of G to row rowIndex of B
                MATRIX::axpy( B + ( row ) * ldaB /* Align with correct row */,
                        G + ( col_idx - columnStart ) * ldaG /* "             " */,
                        1 /* one row */, nColsG,
                        alpha * A._x[ row_idx ] /* Multiply with matrix entry */ );
            }
        }

    return nnz;
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but forms B = G * A
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subMatrixLeftMultiply( const SparseColumnMatrix &A,
                                           const REAL *G, REAL *B,
                                           int rowStart, int columnStart,
                                           int nRows, int nCols, int nRowsG,
                                           bool clear, bool transpose,
                                           int ldaG, int ldaB )
{
    int                row;

    ldaG = max( transpose ? nCols : nRows, ldaG );
    ldaB = max( transpose ? nRows : nCols, ldaB );

    TRACE_ASSERT( rowStart >= 0 && rowStart + nRows <= A._nrow,
            "Row range out of bounds" );
    TRACE_ASSERT( columnStart >= 0 && columnStart + nCols <= A._ncol,
            "Column range out of bounds" );

    if ( clear )
    {
        if ( transpose )
        {
            MATRIX::clear( B, nRowsG, nRows, ldaB );
        }
        else
        {
            MATRIX::clear( B, nRowsG, nCols, ldaB );
        }
    }

    for ( int col_idx = columnStart; col_idx < columnStart + nCols; col_idx++ )
        for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
                row_idx++ )
        {
            row = A._i[ row_idx ];

            if ( row < rowStart || row >= rowStart + nRows )
            {
                continue;
            }

            if ( transpose )
            {
#if 0
                MATRIX::axpy( B + ( col_idx - columnStart ) * ldaB,
                        /* Align with correct row */
                        G + ( row - rowStart ) * ldaG,
                        /* Align with correct row */
                        1 /* one row */, nColsG,
                        A._x[ row_idx ] /* Multiply with matrix entry */ );
#endif
                // Add column colIdx of G to column rowIdx of B
                MATRIX::axpy( B + ( row - rowStart ) /* Align with column */,
                        G + ( col_idx - columnStart ) /* Align with column */,
                        nRowsG /* one column */, 1,
                        A._x[ row_idx ] /* Multiply with matrix entry */,
                        ldaB, ldaG /* Leading dimensions */ );
            }
            else
            {
                // Add column rowIndex of G to column colIndex of B
                MATRIX::axpy( B + ( col_idx - columnStart ) /* Align with column */,
                        G + ( row - rowStart ) /* Align with column */,
                        nRowsG /* one column */, 1,
                        A._x[ row_idx ] /* Multiply with matrix entry */,
                        ldaB, ldaG /* Leading dimensions */ );
#if 0
                MATRIX::axpy( B + ( row - rowStart ) * ldaB /* Align with correct row */,
                        G + ( col_idx - columnStart ) * ldaG /* "             " */,
                        1 /* one row */, nColsG,
                        A._x[ row_idx ] /* Multiply with matrix entry */ );
#endif
            }
        }
}


//////////////////////////////////////////////////////////////////////
// Converts a sparse column matrix to it's graph laplacian
// NOTE: Assumes that the matrix is symmetric
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::convertToGraphLaplacian( SparseColumnMatrix &A,
                                             bool normalize )
{
    for ( int col = 0; col < (int)A._ncol; col++ ) {
        for ( int row_ptr = A._p[ col ]; row_ptr < A._p[ col + 1 ]; row_ptr++ ) {
            if ( A._i[ row_ptr ] == col ) {
                // Set diagonal elements to node degree
                if ( normalize ) {
                    A._x[ row_ptr ] = 1.0;
                } else {
                    A._x[ row_ptr ] = A._p[ col + 1 ] - A._p[ col ];
                }
            } else {
                if ( normalize ) {
                    int row_idx = A._i[ row_ptr ];
                    REAL d1 = (REAL)( A._p[ col + 1 ] - A._p[ col ] );
                    REAL d2 = (REAL)( A._p[ row_idx + 1 ] - A._p[ row_idx ] );

                    A._x[ row_ptr ] = -1.0 / sqrt( d1 * d2 );
                } else {
                    A._x[ row_ptr ] = -1;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::projectInplace( MATRIX& basis, MATRIX& workMat, MATRIX& out )
{
    workMat.clear();
    const int r = basis.cols();
    //const int n = basis.rows();

    map<pair<int,int>, REAL>::const_iterator entryIt;
    for( entryIt = _matrix.begin(); entryIt != _matrix.end(); ++entryIt ) {
        pair<int,int> index = entryIt->first;
        REAL val = entryIt->second;
        int row = index.first;
        int col = index.second;

        for( int i = 0; i < r; i++ ) {
            workMat(i, col) += val * basis(row, i);
        }
    }

    out.clear();

    // SLOW
    out.gemm( 1.0, workMat, basis );

}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool SPARSE_MATRIX::readFromBinary( string filename )
{
    ifstream f( filename.c_str(), ios::binary );
    if( f.fail() )
    {
        cerr << "Couldn't read sparse matrix from binary file: " << filename.c_str() << endl;
        return false;
    }
    else
    {
        int header;
        f.read( (char*)&header, sizeof(header) );

        if( header != -1 )
        {
            TRACE_ERROR("Old format reading not impld");
        }

        f.read( (char*)&_rows, sizeof(_rows) );
        f.read( (char*)&_cols, sizeof(_cols) );
        int nnz = -1;
        f.read( (char*)&nnz, sizeof(nnz) );

        TRACE_ASSERT( nnz > 0 );
        TRACE_ASSERT( _rows > 0 );
        TRACE_ASSERT( _cols > 0 );

        // Write row indices
        vector<int> rows;
        for( int i = 0; i < nnz; i++ )
        {
            int r = -1;
            f.read( (char*)&r, sizeof(r) );
            rows.push_back(r);
        }
#if 0
        vector<int> rows( nnz );
        f.read( (char *)rows.data(), nnz * sizeof( int ) );
#endif

        vector<int> cols;
        for( int i = 0; i < nnz; i++ )
        {
            int c = -1;
            f.read( (char*)&c, sizeof(c) );
            cols.push_back(c);
        }
#if 0
        vector<int> cols( nnz );
        f.read( (char *)cols.data(), nnz * sizeof( int ) );
#endif

        vector<REAL> vals;
        vector< pair< pair<int, int>, REAL > > data( nnz );
        for( int i = 0; i < nnz; i++ )
        {
            REAL v = 0.0;
            f.read( (char*)&v, sizeof(v) );
            vals.push_back(v);

            data[ i ] = pair< pair<int, int>, REAL >(
                    pair<int, int>( rows[ i ], cols[ i ] ), v );
        }
#if 0
        vector<REAL> vals( nnz );
        f.read( (char *)vals.data(), nnz * sizeof( REAL ) );
#endif

#if 0
        // Write values
        for( int i = 0; i < nnz; i++ )
        {
            (*this)( rows.at(i), cols.at(i) ) = vals.at(i);
        }
#endif
        _matrix = map< pair<int, int>, REAL >( data.begin(), data.end() );
    }
    f.close();

    return true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX* SPARSE_MATRIX::makeRandom( int rows, int cols, int nnz )
{
    SPARSE_MATRIX* A = new SPARSE_MATRIX( rows, cols );

    for( int i = 0; i < nnz; i++ ) {
        int row = rand() % rows;
        int col = rand() % cols;
        (*A)(row,col) = (REAL)rand()/RAND_MAX;
    }

    return A;
}

//////////////////////////////////////////////////////////////////////
// Writes a sparse column matrix to disk
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::writeToBinary( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                   const char *filename,
                                   bool newVersion )
{
    FILE* file;
    file = fopen( filename, "wb" );

    if( file == NULL )
    {
        printf( "** WARNING ** Could not write matrix to %s\n", filename );
        return;
    }

    // Write entry count, rows and columns
    if ( newVersion ) {
#if 0
        int header = -1;
        fwrite( (void *)&header, sizeof( int ), 1, file );
#endif
        fwrite( (void *)&( A._nrow ), sizeof( int ), 1, file );
        fwrite( (void *)&( A._ncol ), sizeof( int ), 1, file );
        fwrite( (void *)&( A._nzmax ), sizeof( int ), 1, file );

        // Write row indices
        fwrite( (void *)A._i, sizeof( int ), A._nzmax, file );

        // Write column indices
        for ( size_t col_idx = 0; col_idx < A._ncol; col_idx++ ) {
            int numRows = A._p[ col_idx + 1 ] - A._p[ col_idx ];

            for ( int row_idx = 0; row_idx < numRows; row_idx++ ) {
                fwrite( (void *)&col_idx, sizeof( int ), 1, file );
            }
        }

        // Write values
        fwrite( (void *)A._x, sizeof( double ), A._nzmax, file );
    } else {
        fwrite( (void *)&( A._nzmax ), sizeof( size_t ), 1, file );
        fwrite( (void *)&( A._nrow ), sizeof( size_t ), 1, file );
        fwrite( (void *)&( A._ncol ), sizeof( size_t ), 1, file );

        // Write columns, rows and data
        fwrite( (void *)A._p, sizeof( int ), A._ncol + 1, file );
        fwrite( (void *)A._i, sizeof( int ), A._nzmax, file );
        fwrite( (void *)A._x, sizeof( REAL ), A._nzmax, file );
    }

    fclose( file );
}

//////////////////////////////////////////////////////////////////////
// Reads a sparse column matrix from disk
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::readFromBinary( SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const char *filename )
{
    FILE* file;
    file = fopen( filename, "rb" );

    if( file == NULL )
    {
        printf( "** WARNING ** Could not read matrix from %s\n", filename );
        return;
    }

    // Read entry count, rows and columns
    fread( (void *)&( A._nzmax ), sizeof( size_t ), 1, file );
    fread( (void *)&( A._nrow ), sizeof( size_t ), 1, file );
    fread( (void *)&( A._ncol ), sizeof( size_t ), 1, file );

    // Allocate
    if ( A._p )
    {
        free( A._p );
    }
    if ( A._i )
    {
        free( A._i );
    }
    if ( A._x )
    {
        free( A._x );
    }

#if 0
    A._p = new int[ A._ncol + 1 ];
    A._i = new int[ A._nzmax ];
    A._x = new REAL[ A._nzmax ];
#endif
    A._p = (int *)malloc( ( A._ncol + 1 ) * sizeof( int ) );
    A._i = (int *)malloc( A._nzmax * sizeof( int ) );
    A._x = (REAL *)malloc( A._nzmax * sizeof( REAL ) );

    // Read columns, rows and data
    fread( (void *)A._p, sizeof( int ), A._ncol + 1, file );
    fread( (void *)A._i, sizeof( int ), A._nzmax, file );
    fread( (void *)A._x, sizeof( REAL ), A._nzmax, file );

    cout << SDUMP( A._nzmax ) << endl;
    cout << SDUMP( A._nrow ) << endl;
    cout << SDUMP( A._ncol ) << endl;

    // FIXME
    cout << "First few entries from A._i:" << endl;
    for ( int i = 0; i < 10; i++ ) {
        cout << SDUMP( i ) << SDUMP( A._i[ i ] ) << SDUMP( A._x[ i ] ) << endl;
    }

    fclose( file );
}

//////////////////////////////////////////////////////////////////////
// Computes Gaussian elimination fill-in for a matrix stored in
// sparse column format.
//
// Returns the number of non-zeros in the factor matrix
//////////////////////////////////////////////////////////////////////
long int SPARSE_MATRIX::computeFillIn(
        const SPARSE_MATRIX::SparseColumnMatrix &A )
{
    TRACE_ASSERT( A._ncol == A._nrow, "Matrix must be square" );

    // Initialize row sets for each column
    vector<set<int> >          rowSets( A._ncol );

    long int                   nonZeros = 0;

    for ( size_t col_idx = 0; col_idx < A._ncol; col_idx++ )
    {
        set<int>                &columnSet = rowSets.at( col_idx );

#if 0
        if ( col_idx % 100 == 0 )
        {
#endif
            printf( "Column %lu of %d\n", col_idx + 1, (int)A._ncol );
#if 0
        }
#endif

        for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
                row_ptr++ )
        {
            int                    row_idx = A._i[ row_ptr ];

            columnSet.insert( row_idx );
        }

        for ( set<int>::const_iterator i = columnSet.begin();
                i != columnSet.end(); i++ )
        {
            set<int>              &ancestorSet = rowSets.at( *i );

            set<int>::const_iterator j( i );
            j++;

            for ( ; j != columnSet.end(); j++ )
            {
                ancestorSet.insert( *j );
            }
        }

#if 0
        // Propagate fill-in from this column
        set<int>              &ancestorSet = rowSets.at( row_idx );

        for ( int next_ptr = row_ptr + 1; next_ptr < A._p[ col_idx + 1 ];
                next_ptr++ )
        {
            int                  next_idx = A._i[ next_ptr ];

            ancestorSet.insert( next_idx );
        }
#endif

        // Count total non-zeros in this column
        nonZeros += columnSet.size();
    }
    printf( "\n" );

    return nonZeros;
}
