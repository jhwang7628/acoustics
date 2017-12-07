#include "modal_model/SparseModalEncoder.h" 
#include "wavesolver/MAC_Grid.h"
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
    // square of epsilon
    {
        _error_sqr_target = pow(SparseModalEncoder::epsilon,2);
        std::cout << "Set error for sparse encoder using ||a - tilde{a}||/||U|| metric:\n";
    }
    //// square of epsilon and average across modes
    //{
    //    _error_sqr_target = pow(SparseModalEncoder::epsilon*(REAL)N_Modes(),2);
    //    std::cout << "Set error for sparse encoder using ||a - tilde{a}||/N||U|| metric:\n";
    //}
}

//##############################################################################
//##############################################################################
void SparseModalEncoder::
ComputeWeights()
{
    const auto *interp = FreqWeighting::ISO226_Interpolator::Instance(); 
    const int N_modes = N_Modes(); 
    for (int ii=0; ii<N_modes; ++ii)
        _wi[ii] = interp->Weight(_eigenFreqs[ii]); 
    const REAL wsum = _wi.sum();
    _wi /= wsum; 
    std::cout << "Wi = " << _wi << std::endl;
    //interp->Test(); 
}

//##############################################################################
//##############################################################################
void SparseModalEncoder::
SetWeights(const Eigen::VectorXd &wi)
{
    if (wi.size() != _wi.size())
        throw std::runtime_error("**ERROR** sparse modal encoder weights wrong dim"); 
    _wi = wi; 
}

//##############################################################################
// Function LeastSquareSolve
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
    _Delta_q   = (q-_Q_tilde.A*_c); 
    _error_lsq = (_wi.array()*_Delta_q.array()).array().square(); 
    _weighted_two_norm_sqr_q = (_wi.array()*q.array()).matrix().squaredNorm(); 
    std::cout << "success = " << lsq_solution_valid << "\n"; 
    std::cout << "c       = " << _c         << "\n"; 
    std::cout << "error   = " << _error_lsq << "\n"; 
    std::cout << "===== Least-Square Minimization END =====\n"; 
}

//##############################################################################
// Function MinimizeSparseUpdate
//##############################################################################
int SparseModalEncoder:: 
MinimizeSparseUpdate()
{
    PRINT_FUNC_HEADER;
    _delta_q.resize(N_Modes());  // clear all entries
    double E = _error_lsq.sum(); 
    // using ||Wq|| to normalize
    //const REAL error_sqr_target_relative = _error_sqr_target*_weighted_two_norm_sqr_q;
    // using ||W \Delta q|| to normalize
    const REAL error_sqr_target_relative = _error_sqr_target*E; 
    int ii; 
    _q_sparsity=0;
    std::cout << "===== l1 Minimization START =====\n"; 
    std::cout << " current error  = " << E                         << std::endl; 
    std::cout << " target  error  = " << error_sqr_target_relative << std::endl; 
    while (E >= error_sqr_target_relative && _q_sparsity<N_Modes())
    {
        _error_lsq.maxCoeff(&ii); 
        _delta_q.insert(ii) = _Delta_q[ii]; 
        E -= _error_lsq[ii]; 
        _error_lsq[ii] = 0.0;
        std::cout << "Iteration " << _q_sparsity << ": error = " << E << "; turn on: " << ii << "\n"; 
        ++_q_sparsity;
    }
    std::cout << " N_modes = " << N_Modes() << "\n";
    std::cout << " sparsity = " << (REAL)_q_sparsity/(REAL)N_Modes() << "\n";
    std::cout << "===== l1 Minimization END =====\n"; 
    return _q_sparsity;
}

//##############################################################################
//##############################################################################
void SparseModalEncoder:: 
Encode(const Eigen::VectorXd &q, SimpleTimer *timers) 
{
    PRINT_FUNC_HEADER;
    if (!timers)
    {
        LeastSquareSolve(q); 
        MinimizeSparseUpdate(); 
        _q_aprox = _Q_tilde.A*_c;
        _q_aprox += _delta_q; 
        _Q_tilde.UpdateBasis(_q_aprox);
        _a_buf = _A_tilde.A*_c+_U*_delta_q; 
        _A_tilde.UpdateBasis(_a_buf);
        _q_exact = q; 
    }
    else
    {
        timers[0].Start(); 
        LeastSquareSolve(q); 
        timers[0].Pause(); 
        timers[1].Start(); 
        MinimizeSparseUpdate(); 
        timers[1].Pause(); 
        timers[2].Start(); 
        _q_aprox = _Q_tilde.A*_c;
        _q_aprox += _delta_q; 
        _Q_tilde.UpdateBasis(_q_aprox);
        timers[2].Pause(); 
        timers[3].Start(); 
        _a_buf = _A_tilde.A*_c+_U*_delta_q; 
        timers[3].Pause(); 
        timers[4].Start(); 
        _A_tilde.UpdateBasis(_a_buf);
        timers[4].Pause(); 
        _q_exact = q; 
    }
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
Debug_WriteQComparison(const std::string &filename, 
                       const std::vector<int> &qIndices) const
{
    std::ofstream stream(filename.c_str(), std::ios::app); 
    for (const auto &q_idx : qIndices) 
    {
        stream << std::setprecision(16) 
               << q_idx           << " " 
               << _q_exact(q_idx) << " " 
               << _q_aprox(q_idx) << " ";
    }
    stream << _q_sparsity << std::endl;
}
