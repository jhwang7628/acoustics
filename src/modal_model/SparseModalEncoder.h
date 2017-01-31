#ifndef SPRASE_MODAL_ENCODER_H 
#define SPRASE_MODAL_ENCODER_H 
#include "config.h"
#include <iostream>
#include <memory>
#include <Eigen/Dense> 

//##############################################################################
// Class Basis_Queue
//##############################################################################
class BasisBuffer
{
    public: 
        Eigen::MatrixXd A;  // col-major by default
        int newestColIndex; 

        void Initialize(const int &rows, const int &cols); 
        void UpdateBasis(const Eigen::VectorXd &a);
};

//##############################################################################
// Class SparseModalEncoder
//   This class uses a low-rank model to replace the dense, expensive
//   modal-displacement encoding process
//                       a = U q, 
//   where a is surface displacement, U is modal matrix, and q is modal coord. 
//   Reference and complete mathematical derivation of this encoder is given in
//   acoustic_cornell_box/Algorithms/sparse_modal_update.pdf
//##############################################################################
class SparseModalEncoder
{
    public: 
        static bool useEncoder; 
        static int  rank; 
        static REAL epsilon;

    private: 
        BasisBuffer     _Q_tilde; 
        BasisBuffer     _A_tilde; 
        Eigen::VectorXd _c; 
        Eigen::VectorXd _delta_q; 
        Eigen::VectorXd _error_lsq; 
        const Eigen::MatrixXd &_U; // modal matrix
        REAL _error_sqr_target;  // depends on modal matrix

    public: 
        SparseModalEncoder(const Eigen::MatrixXd &U)
            : _U(U)
        {
            std::cout << "_U.size() = " << _U.rows() << " " << _U.cols() << std::endl;
            std::cout << SparseModalEncoder::rank << std::endl;
            _Q_tilde.Initialize(U.cols(), SparseModalEncoder::rank); 
            _A_tilde.Initialize(U.rows(), SparseModalEncoder::rank); 
            _c.resize(SparseModalEncoder::rank); 
            _delta_q.resize(N_Modes()); 
            _error_lsq.resize(N_Modes()); 
            SetError(SparseModalEncoder::epsilon); 
        }

        inline int N_Modes()const{return _U.cols();}
        inline int N_DOFs() const{return _U.rows();}
        void SetError(const double &epsilon); 
        void LeastSquareSolve(const Eigen::VectorXd &q); 
        int  MinimizeSparseUpdate(); 
        void Encode(const Eigen::VectorXd &q); 
        REAL Decode(const int &vertexID);

        void Test_PerformanceTest(); 
};
typedef std::shared_ptr<SparseModalEncoder> SparseModalEncoderPtr; 

#endif
