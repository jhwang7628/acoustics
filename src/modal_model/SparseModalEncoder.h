#ifndef SPRASE_MODAL_ENCODER_H 
#define SPRASE_MODAL_ENCODER_H 
#include "config.h"
#include "utils/FreqWeighting.hpp"
#include "utils/SimpleTimer.h"
#include <iostream>
#include <memory>
#include <Eigen/Dense> 
#include <Eigen/SparseCore>

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
        typedef Eigen::SparseVector<REAL> SparseVectord;

    private: 
        BasisBuffer     _Q_tilde; 
        BasisBuffer     _A_tilde; 
        Eigen::VectorXd _c; 
        Eigen::VectorXd _a_buf; 
        Eigen::VectorXd _q_aprox; 
        Eigen::VectorXd _q_exact; 
        SparseVectord   _delta_q;  // sparse update
        Eigen::VectorXd _Delta_q;  // lsq error q-Qc
        Eigen::VectorXd _error_lsq; // W(q-Qc) squared
        Eigen::VectorXd _wi; // modal weighting
        REAL            _error_sqr_target;  // depends on modal matrix
        REAL            _weighted_two_norm_sqr_q; 
        int             _q_sparsity; 
        const Eigen::MatrixXd &_U; // modal matrix
        const Eigen::VectorXd _eigenFreqs; 

    public: 
        SparseModalEncoder(const Eigen::MatrixXd &U, const Eigen::VectorXd &eigenFreqs)
            : _U(U), _eigenFreqs(eigenFreqs)
        {
            std::cout << "_U.size() = " << _U.rows() << " " << _U.cols() << std::endl;
            std::cout << SparseModalEncoder::rank << std::endl;
            _Q_tilde.Initialize(U.cols(), SparseModalEncoder::rank); 
            _A_tilde.Initialize(U.rows(), SparseModalEncoder::rank); 
            _c.resize(SparseModalEncoder::rank); 
            _delta_q.resize(N_Modes()); 
            _Delta_q.resize(N_Modes()); 
            _a_buf.resize(N_DOFs()); 
            _error_lsq.resize(N_Modes()); 
            _wi.setOnes(N_Modes()); 
            SetError(SparseModalEncoder::epsilon); 
            ComputeWeights(); 
        }

        inline int N_Modes()const{return _U.cols();}
        inline int N_DOFs() const{return _U.rows();}
        inline const SparseVectord &GetdQ()const{return _delta_q;}
        inline const Eigen::VectorXd &GetDQ()const{return _Delta_q;}
        inline const Eigen::VectorXd &GetQExact()const{return _q_exact;}
        inline const Eigen::VectorXd &GetQAprox()const{return _q_aprox;}
        void SetError(const double &epsilon); 
        void ComputeWeights(); 
        void SetWeights(const Eigen::VectorXd &wi);
        void LeastSquareSolve(const Eigen::VectorXd &q); 
        int  MinimizeSparseUpdate(); 
        void Encode(const Eigen::VectorXd &q, SimpleTimer *timers=nullptr); 
        REAL Decode(const int &row);

        ///// debug /////
        void Debug_WriteQComparison(const std::string &filename, 
                                    const std::vector<int> &qIndices) const; 
};
typedef std::shared_ptr<SparseModalEncoder> SparseModalEncoderPtr; 

#endif
