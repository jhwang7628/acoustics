//////////////////////////////////////////////////////////////////////
// STLUtil.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef STL_UTIL_H
#define STL_UTIL_H

#include <TYPES.h>

#include <vector>

//////////////////////////////////////////////////////////////////////
// Some random utility functions that I needed a place for
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Deletes all contents in a vector of pointers
//////////////////////////////////////////////////////////////////////
template <class T>
inline void clearVectorContents( std::vector<T *> &data )
{
    for ( int i = 0; i < data.size(); i++ )
    {
        delete data[ i ];
    }
}

//////////////////////////////////////////////////////////////////////
// General write function for a vector.  I assume this works.
//////////////////////////////////////////////////////////////////////
template <class T>
inline void writeVector( const char *filename, const vector<T> &data )
{
    FILE* file;
    file = fopen( filename, "wb" );

    size_t bytes_written;

    if ( file == NULL )
    {
        printf( "** WARNING ** Could not write vector to %s\n", filename );
        return;
    }

    int size = data.size();

    // Write dimensions
    bytes_written = fwrite( (void*)&size, sizeof( int ), 1, file );

    // Write data
    bytes_written = fwrite( (void*)data.data(), sizeof(T), size, file );
    fclose( file );
}

//////////////////////////////////////////////////////////////////////
// General write function for a vector.  I assume this works.
//////////////////////////////////////////////////////////////////////
template <class T>
inline void writeVector( FILE *file, const vector<T> &data )
{
    size_t bytes_written;

    if ( file == NULL )
    {
        printf( "** WARNING ** Invalid file handle\n" );
        return;
    }

    int size = data.size();

    // Write dimensions
    bytes_written = fwrite( (void*)&size, sizeof( int ), 1, file );

    // Write data
    bytes_written = fwrite( (void*)data.data(), sizeof(T), size, file );
}

#if 0
//////////////////////////////////////////////////////////////////////
// Write a list of vectors to a file
//////////////////////////////////////////////////////////////////////
inline void writeVector3Array( const char *filename, Vector3Array &data )
{
    FILE* file;
    file = fopen( filename, "wb" );

    size_t bytes_written;

    if ( file == NULL )
    {
        printf( "** WARNING ** Could not write vector to %s\n", filename );
        return;
    }

    int size = data.size();

    // Write dimensions
    bytes_written = fwrite( (void*)&size, sizeof( int ), 1, file );

    // Write data
    bytes_written = fwrite( (void*)data.data(), sizeof(VEC3F), size, file );
    fclose( file );
}
#endif

//////////////////////////////////////////////////////////////////////
// Writes vector data to a file
//////////////////////////////////////////////////////////////////////
inline void writeRealVector(const char *filename, const vector<double> &data)
{
    FILE* file;
    file = fopen(filename, "wb");

    size_t bytes_written;

    if( file == NULL )
    {
        printf( "** WARNING ** Could not write vector to %s\n", filename );
        return;
    }

    int size = data.size();

    // write dimensions
    bytes_written = fwrite((void*)&size, sizeof(int), 1, file);

    // write data
    bytes_written = fwrite((void*)data.data(), sizeof(double), size, file);
    fclose(file);
}

#if 0
//////////////////////////////////////////////////////////////////////
// Write a complex vector to a file
//////////////////////////////////////////////////////////////////////
inline void writeComplexVector(const char *filename, const vector<Complex> &data)
{
    FILE* file;
    file = fopen(filename, "wb");

    size_t bytes_written;

    if( file == NULL )
    {
        printf( "** WARNING ** Could not write vector to %s\n", filename );
        return;
    }

    // Ugly hack - write real and imaginary parts separately
    FloatArray realParts( data.size() );
    FloatArray imaginaryParts( data.size() );

    for ( int i = 0; i < data.size(); i++ )
    {
        realParts[ i ] = real( data[ i ] );
        imaginaryParts[ i ] = imag( data[ i ] );
    }

    int size = data.size();

    // write dimensions
    bytes_written = fwrite((void*)&size, sizeof(int), 1, file);

    // write data
    bytes_written = fwrite((void*)realParts.data(), sizeof(REAL), size, file);
    bytes_written = fwrite((void*)imaginaryParts.data(), sizeof(REAL), size, file);
    fclose(file);
}
#endif

//////////////////////////////////////////////////////////////////////
// General vector read
//////////////////////////////////////////////////////////////////////
template <class T>
inline void readVector( const char *filename, vector<T> &data )
{
    FILE *file;
    file = fopen(filename, "rb");

    size_t bytes_read;

    if ( file == NULL )
    {
        printf( "** WARNING ** Could not read vector %s\n", filename );
        return;
    }

    int size;

    // Read size
    bytes_read = fread((void*)&size, sizeof(int), 1, file);

    data.resize( size );

    // Read data
    bytes_read = fread((void*)data.data(), sizeof(T), size, file);
    fclose(file);
}

//////////////////////////////////////////////////////////////////////
// General vector read
//////////////////////////////////////////////////////////////////////
template <class T>
inline void readVector( FILE *file, vector<T> &data )
{
    size_t bytes_read;

    if ( file == NULL )
    {
        printf( "** WARNING ** Invalid file handle\n" );
        return;
    }

    int size;

    // Read size
    bytes_read = fread((void*)&size, sizeof(int), 1, file);

    data.resize( size );

    // Read data
    bytes_read = fread((void*)data.data(), sizeof(T), size, file);
}

#if 0
//////////////////////////////////////////////////////////////////////
// Reads a list of 3-vectors from a file
//////////////////////////////////////////////////////////////////////
inline void readVector3Array( const char *filename, Vector3Array &data )
{
    FILE *file;
    file = fopen(filename, "rb");

    size_t bytes_read;

    if ( file == NULL )
    {
        printf( "** WARNING ** Could not read vector %s\n", filename );
        return;
    }

    int size;

    // Read size
    bytes_read = fread((void*)&size, sizeof(int), 1, file);

    data.resize( size );

    // Read data
    bytes_read = fread((void*)data.data(), sizeof(VEC3F), size, file);
    fclose(file);
}
#endif

//////////////////////////////////////////////////////////////////////
// Reads vector data from a file
//////////////////////////////////////////////////////////////////////
inline void readRealVector(const char *filename, vector<double> &data)
{
    FILE *file;
    file = fopen(filename, "rb");

    size_t bytes_read;

    if ( file == NULL )
    {
        printf( "** WARNING ** Could not read vector %s\n", filename );
        return;
    }

    int size;

    // Read size
    bytes_read = fread((void*)&size, sizeof(int), 1, file);

    data.resize( size );

    // Read data
    bytes_read = fread((void*)data.data(), sizeof(double), size, file);
    fclose(file);
}

#if 0
//////////////////////////////////////////////////////////////////////
// Reads vector data from a file
//////////////////////////////////////////////////////////////////////
inline void readComplexVector(const char *filename, vector<Complex> &data)
{
    FILE *file;
    file = fopen(filename, "rb");

    size_t bytes_read;

    if ( file == NULL )
    {
        printf( "** WARNING ** Could not read vector %s\n", filename );
        return;
    }

    // Ugly hack, read real and imaginary parts separately
    FloatArray realParts( data.size() );
    FloatArray imaginaryParts( data.size() );

    int size;

    // Read size
    bytes_read = fread((void*)&size, sizeof(int), 1, file);

    data.resize( size );
    realParts.resize( size );
    imaginaryParts.resize( size );

    // Read data
    bytes_read = fread((void*)realParts.data(), sizeof(Real), size, file);
    bytes_read = fread((void*)imaginaryParts.data(), sizeof(Real), size, file);
    fclose(file);

    for ( int i = 0; i < size; i++ )
    {
        real( data[ i ] ) = realParts[ i ];
        imag( data[ i ] ) = imaginaryParts[ i ];
    }
}
#endif

//////////////////////////////////////////////////////////////////////
// Sum of entries in a vector
//////////////////////////////////////////////////////////////////////
template <class T>
inline T vectorSum( const std::vector<T> &data )
{
    T sum = 0;

    for ( int i = 0; i < data.size(); i++ )
    {
        sum += data[ i ];
    }
}

#endif
