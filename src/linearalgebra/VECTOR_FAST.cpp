// VECTOR.h: interface for the VECTOR class.
//
//////////////////////////////////////////////////////////////////////

#include "VECTOR.h"
#ifdef USE_MKL
#include <mkl_types.h>
#include <mkl_cblas.h>
#endif
#ifdef USING_ATLAS
extern "C" {
#include <cblas.h>
}
#endif

#include <memory.h>

//////////////////////////////////////////////////////////////////////
// Constructor for the full vector
//////////////////////////////////////////////////////////////////////
VECTOR::VECTOR() :
    _size(0)
{
    _vector = NULL;
}

VECTOR::VECTOR(int size, REAL val) :
    _size(size)
{
    _vector = new REAL[_size];
    clear();

    if ( val != 0 ) {
        for ( int i = 0; i < size; i++ ) {
            _vector[ i ] = val;
        }
    }
}

VECTOR::VECTOR(const char* filename)
{
    _vector = NULL;
    read(filename);
}

VECTOR::VECTOR(const VECTOR& v) 
{
    _size = v._size;
    _vector = new REAL[_size];
    for (int x = 0; x < _size; x++)
        _vector[x] = v._vector[x];
}

VECTOR::VECTOR(FILE* file)
{
    // read dimensions
    fread((void*)&_size, sizeof(int), 1, file);

    // read data
    _vector = new REAL[_size];
    if (sizeof(REAL) == sizeof(float))
    {
        double* vecDouble = new double[_size];
        fread((void*)vecDouble, sizeof(double), _size, file);
        for (int x = 0; x < _size; x++)
            _vector[x] = vecDouble[x];
        delete[] vecDouble;
    }
    else
        fread((void*)_vector, sizeof(double), _size, file);
}

VECTOR::~VECTOR()
{
    delete[] _vector;
}

//////////////////////////////////////////////////////////////////////
// dot product with another vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::operator*(const VECTOR& vector)
{
#if BOUNDS_CHECKING_ENABLED
    if (vector._size != _size)
    {
        cout << __FILE__ << " " << __LINE__ << " VECTOR dot sizes do not match!: " << endl;
        return 123456.0f;
    }
#endif
#ifdef SINGLE_PRECISION
    return cblas_sdot(_size, _vector, 1, vector._vector, 1);
#else
    return cblas_ddot(_size, _vector, 1, vector._vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// dot product with another vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::operator*(const REAL *vector)
{
#ifdef SINGLE_PRECISION
    return cblas_sdot(_size, _vector, 1, vector, 1);
#else
    return cblas_ddot(_size, _vector, 1, vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// scale vector by a constant
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator*=(const REAL& alpha)
{
#ifdef SINGLE_PRECISION
    cblas_sscal(_size, alpha, _vector, 1);
#else
    cblas_dscal(_size, alpha, _vector, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator/=(const REAL& alpha)
{
#ifdef SINGLE_PRECISION
    cblas_sscal(_size, 1.0/alpha, _vector, 1);
#else
    cblas_dscal(_size, 1.0/alpha, _vector, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
// wipe the whole vector
//////////////////////////////////////////////////////////////////////
void VECTOR::clear()
{
    memset(_vector, 0, _size * sizeof(REAL));
}

//////////////////////////////////////////////////////////////////////
// resize and wipe the vector
//////////////////////////////////////////////////////////////////////
void VECTOR::resizeAndWipe(int size)
{
    if (_size != size)
        delete[] _vector;
    _size = size;

    _vector = new REAL[_size];
    clear();
}

//////////////////////////////////////////////////////////////////////
// write the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::write(const char* filename)
{
    FILE* file;
    file = fopen(filename, "wb");

    if( file == NULL )
    {
        printf( "** WARNING ** Could not write vector to %s\n", filename );
        return;
    }

    // write dimensions
    fwrite((void*)&_size, sizeof(int), 1, file);

    // write data
    if (sizeof(REAL) == sizeof(float))
    {
        double* vecDouble = new double[_size];
        for (int x = 0; x < _size; x++)
            vecDouble[x] = _vector[x];
        fwrite((void*)vecDouble, sizeof(double), _size, file);
        delete[] vecDouble;
    } 
    else
        fwrite((void*)_vector, sizeof(REAL), _size, file);
    fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::write(FILE* file)
{
    // write dimensions
    fwrite((void*)&_size, sizeof(int), 1, file);

    // write data
    if (sizeof(REAL) == sizeof(float))
    {
        double* vecDouble = new double[_size];
        for (int x = 0; x < _size; x++)
            vecDouble[x] = _vector[x];
        fwrite((void*)vecDouble, sizeof(double), _size, file);
        delete[] vecDouble;
    } 
    else
        fwrite((void*)_vector, sizeof(REAL), _size, file);
}

//////////////////////////////////////////////////////////////////////
// read vector from a file
//////////////////////////////////////////////////////////////////////
bool VECTOR::read(const char* filename)
{
    FILE* file;
    file = fopen(filename, "rb");
    if (file == NULL)
    {
        cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
        return false;
    }

    // read dimensions
    fread((void*)&_size, sizeof(int), 1, file);

    // read data
    delete[] _vector;
    _vector = new REAL[_size];

    if (sizeof(REAL) == sizeof(float))
    {
        double* vecDouble = new double[_size];
        fread((void*)vecDouble, sizeof(double), _size, file);
        for (int x = 0; x < _size; x++)
            _vector[x] = vecDouble[x];
        delete[] vecDouble;
    }
    else
        fread((void*)_vector, sizeof(REAL), _size, file);
    fclose(file);

    return true;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(VECTOR m)
{
    if (m.size() != _size)
    {
        delete[] _vector;
        _size= m.size();
        _vector= new REAL[_size];
    }
    memcpy (_vector, m._vector, _size * sizeof(REAL));
    return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-add operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator+=(const VECTOR& m)
{
#if BOUNDS_CHECKING_ENABLED
    if (m._size != _size)
        cout << __FILE__ << " " << __LINE__ << " : Vector sizes don't match! " << endl;
#endif
#ifdef SINGLE_PRECISION
    cblas_saxpy (_size, 1.0f, m._vector, 1, _vector, 1);
#else
    cblas_daxpy (_size, 1.0, m._vector, 1, _vector, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-subtract operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator-=(const VECTOR& m)
{
#if BOUNDS_CHECKING_ENABLED
    if (m._size != _size)
    {
        cout << __FILE__ << " " << __LINE__ << " : Vector sizes don't match! " << endl;
        cout << " this: " << _size << endl;
        cout << " input: " << m._size << endl;
    }
#endif
#ifdef SINGLE_PRECISION
    cblas_saxpy (_size, -1.0f, m._vector, 1, _vector, 1);
#else
    cblas_daxpy (_size, -1.0, m._vector, 1, _vector, 1);
#endif
    return *this;
}

//////////////////////////////////////////////////////////////////////
// compute the 2 norm
//////////////////////////////////////////////////////////////////////
REAL VECTOR::norm2()
{
#ifdef SINGLE_PRECISION
    return cblas_snrm2( _size, _vector, 1 );
#else
    return cblas_dnrm2( _size, _vector, 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// largest element in vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::max_()
{
    REAL vmax = _vector[0];

    for ( int i = 1; i < _size; i++ )
    {
        if ( _vector[i] > vmax )
            vmax = _vector[i];
    }

    return vmax;
}

//////////////////////////////////////////////////////////////////////
// smallest element in vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::min_()
{
    REAL vmin = _vector[0];

    for ( int i = 1; i < _size; i++ )
    {
        if ( _vector[i] < vmin )
            vmin = _vector[i];
    }

    return vmin;
}

//////////////////////////////////////////////////////////////////////
// largest absolute element in vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::absmax()
{
    REAL vmax = std::fabs( _vector[0] );

    for ( int i = 1; i < _size; i++ )
    {
        if ( std::fabs( _vector[i] ) > vmax )
            vmax = std::fabs( _vector[i] );
    }

    return vmax;
}

//////////////////////////////////////////////////////////////////////
// smallest absolute element in vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::absmin()
{
    REAL vmin = std::fabs( _vector[0] );

    for ( int i = 1; i < _size; i++ )
    {
        if ( std::fabs( _vector[i] ) < vmin )
            vmin = std::fabs( _vector[i] );
    }

    return vmin;
}

//////////////////////////////////////////////////////////////////////
// swap contents with another vector
//////////////////////////////////////////////////////////////////////
void VECTOR::swap(VECTOR& vec)
{
    REAL* temp = _vector;
    _vector = vec._vector;
    vec._vector = temp;
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: y += alpha * x, where y is this vector
//////////////////////////////////////////////////////////////////////
void VECTOR::axpy(REAL alpha, VECTOR& x)
{
#if BOUNDS_CHECKING_ENABLED
    if (_size != x._size)
        cout << __FILE__ << " " << __LINE__ << " : Vector size do not match! " << endl;
#endif
#ifdef SINGLE_PRECISION
    cblas_saxpy (_size, alpha, x._vector, 1, _vector, 1);
#else
    cblas_daxpy (_size, alpha, x._vector, 1, _vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// same as axpy above, but vector contents are stomped as well
//////////////////////////////////////////////////////////////////////
void VECTOR::clearingAxpy(REAL alpha, VECTOR& x)
{
#if BOUNDS_CHECKING_ENABLED
    if (_size != x._size)
        cout << __FILE__ << " " << __LINE__ << " : Vector size do not match! " << endl;
#endif
    memset(_vector, 0, _size * sizeof(REAL));
#ifdef SINGLE_PRECISION
    cblas_saxpy (_size, alpha, x._vector, 1, _vector, 1);
#else
    cblas_daxpy (_size, alpha, x._vector, 1, _vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// return the sum of the vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::sum()
{
#ifdef SINGLE_PRECISION
    return cblas_sasum(_size, _vector, 1);
#else
    return cblas_dasum(_size, _vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// in-place copy, since operator= must allocate a new VECTOR
//////////////////////////////////////////////////////////////////////
void VECTOR::copyInplace(const VECTOR& vector)
{
#if BOUNDS_CHECKING_ENABLED
    if (_size != vector._size)
        cout << __FILE__ << " " << __LINE__ << " : Vector sizes do not match! " << endl;
#endif
    memcpy (_vector, vector._vector, _size * sizeof(REAL));
}

//////////////////////////////////////////////////////////////////////
// in-place copy, starting from a given index
//////////////////////////////////////////////////////////////////////
void VECTOR::copyInplace(const VECTOR& vector, int start)
{
#if BOUNDS_CHECKING_ENABLED
    if ( (_size - start) < vector._size)
        cout << __FILE__ << " " << __LINE__ << " : Vector too small for copy! " << endl;
#endif
    memcpy (_vector + start, vector._vector, vector._size * sizeof(REAL));
}

//////////////////////////////////////////////////////////////////////
// Copies a sub vector from the input vector in to the
// given location in this vector.
//////////////////////////////////////////////////////////////////////
void VECTOR::copySubVectorInplace(const VECTOR &vector,
                                  int startLoc_in, int startLoc_out,
                                  int size)
{
    if ( size <= 0 )
        return;

#if BOUNDS_CHECKING_ENABLED
    if ( startLoc_in + size > vector._size )
    {
        cout << __FILE__ << " " << __LINE__ << " : Index out of bounds in"
            "input vector! " << endl;
    }
    if ( startLoc_out + size > _size )
    {
        cout << __FILE__ << " " << __LINE__ << " : Index out of bounds in"
            "output vector! " << endl;
    }
#endif
    memcpy(_vector + startLoc_out, vector._vector + startLoc_in,
            size * sizeof(REAL));
}

//////////////////////////////////////////////////////////////////////
// add two vectors
//////////////////////////////////////////////////////////////////////
VECTOR operator+(VECTOR& x, VECTOR& y) 
{
    VECTOR z(x.size());

    for (int i = 0; i < x.size(); i++)
        z(i) = x(i) + y(i);
    return z;
}

//////////////////////////////////////////////////////////////////////
// subtract two vectors
//////////////////////////////////////////////////////////////////////
VECTOR operator-(VECTOR& x, VECTOR& y) 
{
    VECTOR z(x.size());

    for (int i = 0; i < x.size(); i++)
        z(i) = x(i) - y(i);
    return z;
}

//////////////////////////////////////////////////////////////////////
// scale a vector
//////////////////////////////////////////////////////////////////////
VECTOR operator*(VECTOR& x, REAL& scalar) 
{
    VECTOR z(x.size());

    for (int i = 0; i < x.size(); i++)
        z(i) = x(i) * scalar;
    return z;
}

//////////////////////////////////////////////////////////////////////
// scale a vector
//////////////////////////////////////////////////////////////////////
VECTOR operator*(REAL& scalar, VECTOR& x) 
{
    VECTOR z(x.size());

    for (int i = 0; i < x.size(); i++)
        z(i) = x(i) * scalar;
    return z;
}

//////////////////////////////////////////////////////////////////////
// compute the 2 norm
//////////////////////////////////////////////////////////////////////
REAL operator^(VECTOR& x, VECTOR& y)
{
    REAL total = 0;
    for (int i = 0; i < x.size(); i++)
        total += x.data()[i] * y.data()[i];
    return total;
}

//////////////////////////////////////////////////////////////////////
// max of the vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::maxValue()
{
    if (_size == 0) return 0;

    REAL maxFound = _vector[0];
    for (int x = 0; x < _size; x++)
        if (_vector[x] > maxFound) maxFound = _vector[x];

    return maxFound;
}

//////////////////////////////////////////////////////////////////////
// min of the vector
//////////////////////////////////////////////////////////////////////
REAL VECTOR::minValue()
{
    if (_size == 0) return 0;

    REAL minFound = _vector[0];
    for (int x = 0; x < _size; x++)
        if (_vector[x] < minFound) minFound = _vector[x];

    return minFound;
}

//////////////////////////////////////////////////////////////////////
// Take the absolute value of all entries
//////////////////////////////////////////////////////////////////////
void VECTOR::fabs()
{
    for (int x = 0; x < _size; x++)
        _vector[x] = std::fabs(_vector[x]);
}

//////////////////////////////////////////////////////////////////////
// compute the root mean squared
//////////////////////////////////////////////////////////////////////
REAL VECTOR::rms()
{
    /*
       REAL sumSq = 0.0;
       for (int x = 0; x < _size; x++)
       sumSq += _vector[x] * _vector[x];
       return sqrt(sumSq / _size);
     */
    return norm2() / sqrt((REAL)_size);
}

//////////////////////////////////////////////////////////////////////
// compute the infinity norm
//////////////////////////////////////////////////////////////////////
REAL VECTOR::normInf()
{
    if (_size == 0) return 0.0;
    REAL max = std::fabs(_vector[0]);
    for (int x = 1; x < _size; x++)
        if (std::fabs(_vector[x]) > max)
            max = std::fabs(_vector[x]);
    return max;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL VECTOR::norm2( REAL *v, int size )
{
#ifdef SINGLE_PRECISION
    return cblas_snrm2( size, v, 1 );
#else
    return cblas_dnrm2( size, v, 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// Vector dot product
//////////////////////////////////////////////////////////////////////
REAL VECTOR::dot( int n, const REAL *x, const REAL *y, int incX, int incY )
{
#ifdef SINGLE_PRECISION
  return cblas_sdot( n, x, incX, y, incY );
#else
  return cblas_ddot( n, x, incX, y, incY );
#endif
}

//////////////////////////////////////////////////////////////////////
// Absolute sum of vector components
//////////////////////////////////////////////////////////////////////
REAL VECTOR::absSum( int n, const REAL *x, int incX )
{
#ifdef SINGLE_PRECISION
  return cblas_sasum( n, x, incX );
#else
  return cblas_dasum( n, x, incX );
#endif
}
