//////////////////////////////////////////////////////////////////////
// test-superlu.cpp: Simple program for testing our SuperLU interface
//
// 
//////////////////////////////////////////////////////////////////////

#include <iostream>

#include <linearalgebra/SPARSE_MATRIX.h>

#include <superlu-interface/SuperLU_Interface.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    SPARSE_MATRIX            A( 3, 3 );

    A( 0, 0 ) = 1.0;
    A( 0, 1 ) = 1.0;

    A( 1, 0 ) = 1.0;
    A( 1, 1 ) = 1.0;
    A( 1, 2 ) = 1.0;

    A( 2, 1 ) = 1.0;
    A( 2, 2 ) = 1.0;

    SuperLU_Solver           solver( A );

    return 0;
}
