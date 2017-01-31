#ifndef SPRASE_MODAL_ENCODER_H 
#define SPRASE_MODAL_ENCODER_H 
#include "config.h"
#include <memory>
#include <Eigen/Dense> 

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
        int             _N_modes; 
        Eigen::MatrixXd _Q_tilde; 
        Eigen::VectorXd _c; 
        Eigen::VectorXd _delta_q; 
        const Eigen::MatrixXd &_U; // modal matrix

    public: 
        SparseModalEncoder(const int &N_modes, const Eigen::MatrixXd &U)
            : _U(U)
        {
            _Q_tilde.setZero(N_modes, SparseModalEncoder::rank); 
            _c.resize(SparseModalEncoder::rank); 
            _delta_q.resize(N_modes); 
        }

        void LeastSquareSolve(const Eigen::VectorXd &q); 
        void MinimizeSparseUpdate(); 
        void Encode(const Eigen::VectorXd &q); 
        REAL Decode(const int &vertexID);

        void Test_PerformanceTest(); 
};
typedef std::shared_ptr<SparseModalEncoder> SparseModalEncoderPtr; 

#endif
