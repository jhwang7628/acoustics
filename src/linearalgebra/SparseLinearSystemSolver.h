#ifndef SPARSE_LINEAR_SYSTEM_SOLVER_H 
#define SPARSE_LINEAR_SYSTEM_SOLVER_H 

#include <queue> 
#include <vector>
#include <iostream>
#include <Eigen/Dense> 
#include <Eigen/Sparse> 


#include <boost/timer/timer.hpp>
#include <tools/unit_testing/testing.h>


// provide a wrapper to Eigen sparse library to solve a sparse linear system 
// Note that the sparse matrix is stored in row major format
class SparseLinearSystemSolver
{
    public: 
        typedef Eigen::SparseMatrix<double, Eigen::RowMajor>    SparseMatrix; 
        //typedef Eigen::SparseMatrix<double>    SparseMatrix; 
        typedef Eigen::SparseVector<double>                     SparseVector;
        typedef Eigen::Triplet<double>                          Triplet; 

        typedef Eigen::SparseMatrix<int, Eigen::RowMajor>       SparseMatrixPattern;  // coupled pattern
        typedef Eigen::Triplet<int>                             TripletInt;

    private: 

        // size of the system
        int                     _N;
        // number of non-zero entries
        int                     _M; 

        // linear system: Ax = b
        SparseMatrix            _A; 
        SparseVector            _b;    
        Eigen::VectorXd         _bDense; 

        std::vector<Triplet>    _stagedEntries; 

        // use these two to figure out the coupling pattern
        SparseMatrixPattern     _A_pattern;
        std::vector<TripletInt> _stagedEntriesPattern; 

        std::vector<int>        _belongComponent; 

    public: 

        inline bool HasUnstaged()   { return (_stagedEntries.size()==0 ? false : true); }
        inline int  NumberNodes()   { return _N; } 
        inline int  NumberEntries() { return _M; } 

        void StageEntry(const int &x_i, const int &x_j, const double &e_ij); 

        // Fill-in of the matrix. needs to be performed before any operation on
        // the matrix 
        void FillIn(const bool &destroyOld=true); 
        void FillInRHS(const int &index, const double &value){ _bDense(index) = value; /* _b.coeffRef(index) = value; */ } 

        void GetAllNonzeroElements(Eigen::MatrixXd &nonzeroElements); 

        // use BFS to find the connected components in the graph and return the
        // lists of nodes of each component. 
        // singleton does not count as component
        void GetAllComponents(std::vector<std::vector<int> > &allComponents); 
        void GetComponent(const int &root, std::vector<int> &component, const int &componentID=-1); 


        void Solve(Eigen::VectorXd &x); 
        void DecoupledSolve(Eigen::VectorXd &x); 
        void SparseSolve(Eigen::VectorXd &x); 
        void BiCGStabSolve(Eigen::VectorXd &x); 

        SparseLinearSystemSolver(const int &N) : 
            _N(N)
            {
                _A.resize(N,N); 
                _bDense = Eigen::VectorXd::Zero(N);
                //_b.resize(N);
                _belongComponent.resize(N); 

                _A_pattern.resize(N,N); 
            }



        //// debugging tools ////

        void PrintMatrixDense(std::ostream &os); 
        void PrintVectorDense(std::ostream &os); 
        void DenseSolve(Eigen::VectorXd &x); 
        void GetAllComponentsVerbose(std::vector<std::vector<int> > &allComponents); 


}; 



#endif
