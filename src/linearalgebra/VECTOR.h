// VECTOR.h: interface for the VECTOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef VECTOR_H
#define VECTOR_H

#include <TYPES.h>

#include <map>
#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////
// An arbitrary dimension vector class
//////////////////////////////////////////////////////////////////////
class VECTOR {

    public:
        VECTOR();
        VECTOR(int size, REAL val = 0);
        VECTOR(int size, const REAL* data);
        VECTOR(const char* filename);
        VECTOR(const VECTOR& v);
        VECTOR(FILE* file);
        ~VECTOR();

        inline REAL& operator()(int index) { return _vector[index]; };
        inline REAL& operator[](int index) { return _vector[index]; };

        inline REAL operator() (int index) const { return _vector[index]; };
        inline REAL operator[] (int index) const { return _vector[index]; };

        int& size() { return _size; };
        int size() const { return _size; };

        // wipe the whole vector
        void clear();

        // write the vector to a binary file
        // everything is always written as a double
        void write(const char* filename);
        void write(FILE* file);

        // read the vector to a binary file
        // everything is always read in as a double, then
        // converted if necessary
        // Returns true if successfully read.
        bool read(const char* filename);

        // resize the vector and wipe to zero
        void resizeAndWipe(int size);

        // overloaded operators
        VECTOR& operator=(VECTOR m);
        VECTOR& operator+=(const VECTOR& v);
        VECTOR& operator-=(const VECTOR& v);
        VECTOR& operator*=(const REAL& alpha);
        VECTOR& operator/=(const REAL& alpha);

        // (stevenan) The trailing underscore is to make this windows-compatible.
        // http://polingplace.blogspot.com/2007/04/stdmin-and-stdmax-versus-visual-studio.html
        // Thanks, Microsoft.
        REAL max_();
        REAL min_();

        REAL absmax();
        REAL absmin();

        // 2 norm
        REAL norm2();

        // Infinity norm
        REAL normInf();

        REAL rms();

        // dot product
        REAL operator*(const VECTOR& vector);
        REAL operator*(const REAL *vector);

        // swap contents with another vector --
        // it's your job to ensure that they are the same size
        void swap(VECTOR& vector);

        // raw data pointer
        REAL* data() { return _vector; };
        const REAL* data() const { return _vector; }

        // BLAS axpy operation: y += alpha * x, where y is this vector
        void axpy(REAL alpha, VECTOR& x);

        // same as axpy above, but vector contents are stomped as well
        void clearingAxpy(REAL alpha, VECTOR& x);

        // sum of all the elements
        REAL sum();

        // in-place copy, since operator= must allocate a new VECTOR
        void copyInplace(const VECTOR& vector);
        void equals(const VECTOR& vector) { copyInplace(vector); }

        // In place copy from one vector to a specified location
        // in this vector.
        void copyInplace(const VECTOR &vector, int startLoc);

        // Copies a sub vector from the input vector in to the
        // given location in this vector.
        void copySubVectorInplace(const VECTOR &vector,
                int startLoc_in, int startLoc_out,
                int size);

        // mean of all the entries
        REAL mean() { return sum() / _size; };

        // max of all the elements
        REAL maxValue();

        // min of all the elements
        REAL minValue();

        // take the absolute value of all entires
        void fabs();

        static REAL norm2( REAL *v, int size );

        // Vector dot product
        static REAL dot( int n, const REAL *x, const REAL *y,
                         int incX = 1, int incY = 1 );

        // Absolute sum of vector components
        static REAL absSum( int n, const REAL *x, int incX = 1 );

    private:
        int _size;
        REAL* _vector;
};

//////////////////////////////////////////////////////////////////////
// dump vector to iostream
//////////////////////////////////////////////////////////////////////
inline ostream &operator<<(ostream &out, VECTOR& vector)
{
    out << "[" << endl; 
    for (int x = 0; x < vector.size(); x++)
        out << vector(x) << endl;
    out << "]" << endl;
    return out;
}

// overloaded operators
VECTOR operator-(const VECTOR& x, const VECTOR& y);
VECTOR operator+(const VECTOR& x, const VECTOR& y);
VECTOR operator*(const VECTOR& x, const REAL& scalar);
VECTOR operator*(const REAL& scalar, const VECTOR& x);

// x^T . y
REAL operator^(const VECTOR& x, const VECTOR& y);

#endif
