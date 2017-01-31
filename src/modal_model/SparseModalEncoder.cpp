#include "modal_model/SparseModalEncoder.h" 
#include "macros.h"
#include <iostream>

//##############################################################################
// Static class members
//##############################################################################
bool SparseModalEncoder::useEncoder = false; 
int  SparseModalEncoder::rank       = 0; 
REAL SparseModalEncoder::epsilon    = 0.; 

//##############################################################################
//##############################################################################
void SparseModalEncoder::
LeastSquareSolve(const Eigen::VectorXd &q)
{
    PRINT_FUNC_HEADER;
}

//##############################################################################
//##############################################################################
void SparseModalEncoder:: 
MinimizeSparseUpdate()
{
    PRINT_FUNC_HEADER;
}

//##############################################################################
//##############################################################################
void SparseModalEncoder:: 
Encode(const Eigen::VectorXd &q) 
{
    PRINT_FUNC_HEADER;
    LeastSquareSolve(q); 
    MinimizeSparseUpdate(); 
}

//##############################################################################
//##############################################################################
REAL SparseModalEncoder:: 
Decode(const int &vertexID) 
{
    PRINT_FUNC_HEADER;
    (void)vertexID;
    return 0.;
}

//##############################################################################
//##############################################################################
void SparseModalEncoder::
Test_PerformanceTest()
{
}
