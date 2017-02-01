#include "modal_model/SparseModalEncoder.h" 
#include "macros.h"
#include <iostream>

//##############################################################################
//##############################################################################
void BasisBuffer:: 
Initialize(const int &rows, const int &cols)
{
    assert(rows>0 && cols>0);
    A.setZero(rows, cols); 
    newestColIndex = cols-1; 
}

//##############################################################################
// Function UpdateBasis
//   Add new basis and discard the oldest
//##############################################################################
void BasisBuffer:: 
UpdateBasis(const Eigen::VectorXd &a)
{
    const int oldestColIndex = (newestColIndex+1) % A.cols(); 
    A.col(oldestColIndex) = a; 
    newestColIndex = oldestColIndex; 
}

//##############################################################################
// Static class members
//##############################################################################
bool SparseModalEncoder::useEncoder               = false; 
int  SparseModalEncoder::rank                     = 0; 
REAL SparseModalEncoder::epsilon                  = 0.; 

//##############################################################################
//##############################################################################
void SparseModalEncoder::
SetError(const double &epsilon)
{
    //// scale by 2norm of U
    //{
    //    // compute 2-norm of the modal matrix
    //    Eigen::JacobiSVD<Eigen::MatrixXd> svd(_U); 
    //    const double U_2norm = svd.singularValues()[0]; 
    //    _error_sqr_target = pow(SparseModalEncoder::epsilon/U_2norm, 2); 
    //    std::cout << "Set error for sparse encoder using ||a - \tilde{a}|| metric:\n";
    //    std::cout << "||U||_2 = " << U_2norm << "\n"; 
    //    std::cout << "error sqr target = " << _error_sqr_target << "\n"; 
    //    std::cout << std::flush; 
    //}
    //// square of epsilon
    //{
    //    _error_sqr_target = pow(SparseModalEncoder::epsilon,2);
    //    std::cout << "Set error for sparse encoder using ||a - tilde{a}||/||U|| metric:\n";
    //}
    // square of epsilon and average across modes
    {
        _error_sqr_target = pow(SparseModalEncoder::epsilon*(REAL)N_Modes(),2);
        std::cout << "Set error for sparse encoder using ||a - tilde{a}||/N||U|| metric:\n";
    }
}

//##############################################################################
//##############################################################################
void SparseModalEncoder::
LeastSquareSolve(const Eigen::VectorXd &q)
{
    PRINT_FUNC_HEADER;
    std::cout << "===== Least-Square Minimization START =====\n"; 
    _c = _Q_tilde.A.colPivHouseholderQr().solve(q); 

    bool lsq_solution_valid = true; 
    for (int ii=0; ii<_c.size(); ++ii)
        lsq_solution_valid = (lsq_solution_valid && !std::isnan(_c[ii]));
    //const bool lsq_solution_valid = (_Q_tilde.A*_c).isApprox(q, _error_sqr_target);
    if (!lsq_solution_valid)
        _c.setZero();
    _error_lsq     = (q-_Q_tilde.A*_c); 
    _error_lsq_abs = _error_lsq.array().abs(); 
    std::cout << "success = " << lsq_solution_valid << "\n"; 
    std::cout << "c       = " << _c         << "\n"; 
    std::cout << "error   = " << _error_lsq << "\n"; 
    std::cout << "===== Least-Square Minimization END =====\n"; 
}

//##############################################################################
//##############################################################################
int SparseModalEncoder:: 
MinimizeSparseUpdate()
{
    PRINT_FUNC_HEADER;
    _delta_q.setZero(); 
    double E = (_error_lsq).array().square().sum(); 
    int ii; 
    int count_sparsity=0;
    std::cout << "===== l1 Minimization START =====\n"; 
    std::cout << " current error  = " <<  E                << std::endl; 
    std::cout << " target  error  = " << _error_sqr_target << std::endl; 
    while (E >= _error_sqr_target && count_sparsity<N_Modes())
    {
        _error_lsq_abs.maxCoeff(&ii); 
        _delta_q[ii] = _error_lsq[ii]; 
        E -= pow(_error_lsq[ii],2); 
        _error_lsq[ii] = 0.0;
        _error_lsq_abs[ii] = 0.0;
        std::cout << "Iteration " << count_sparsity << ": error = " << E << "\n"; 
        //std::cout << " current q= ";
        //for (int jj=0; jj<_error_lsq.size(); ++jj)
        //    std::cout << _error_lsq[jj] << ", "; 
        //std::cout << "\n";
        ++count_sparsity;
    }
    std::cout << " N_modes = " << N_Modes() << "\n";
    std::cout << " sparsity = " << (REAL)count_sparsity/(REAL)N_Modes() << "\n";
    std::cout << "===== l1 Minimization END =====\n"; 
    return count_sparsity;
}

//##############################################################################
//##############################################################################
void SparseModalEncoder:: 
Encode(const Eigen::VectorXd &q) 
{
    PRINT_FUNC_HEADER;
    LeastSquareSolve(q); 
    MinimizeSparseUpdate(); 
    _Q_tilde.UpdateBasis(q);
    _A_tilde.UpdateBasis(_A_tilde.A*_c+_U*_delta_q);  // TODO need to use a sparsity mask on delta_q
}

//##############################################################################
//##############################################################################
REAL SparseModalEncoder:: 
Decode(const int &row) 
{
    const REAL value = _A_tilde.A(row, _A_tilde.newestColIndex); 
    return (std::isnan(value) ? 0. : value); 
}

//##############################################################################
//##############################################################################
void SparseModalEncoder::
Test_PerformanceTest()
{
}
