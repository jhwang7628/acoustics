
#ifndef UNDIRECTED_GRAPH_H 
#define UNDIRECTED_GRAPH_H 


#include <queue> 
#include <vector>
#include <iostream>
#include <Eigen/Dense> 
#include <Eigen/Sparse> 


// A simple graph class that represents graph using adjacency matrix in Eigen
// Sparse matrix (col-major compressed format) 
class UndirectedGraph 
{
    public: 
        typedef Eigen::SparseMatrix<int>    SparseMatrixXi; 
        typedef Eigen::Triplet<int>         TripletInt; 

    private: 

        // number of nodes
        int                     _N;
        // number of edges 
        int                     _M; 

        SparseMatrixXi          _adjacencyMatrix; 
        std::vector<TripletInt> _stagedEdges; 

        std::vector<int>        _belongComponent; 

    public: 

        inline bool HasUnstaged() { return (_stagedEdges.size()==0 ? false : true); }
        inline int  NumberNodes() { return _N; } 
        inline int  NumberEdges() { return _M; } 

        // both (x_i, x_j), (x_j, x_i) will be added to the staged edges. 
        void AddEdgeToStage(const int &x_i, const int &x_j); 

        // put all the staged edges to the graph
        void ConstructEdges(const bool &destroyOld=true); 

        void GetAllNonzeroElements(Eigen::MatrixXd &nonzeroElements); 

        // use BFS to find the connected components in the graph and return the
        // lists of nodes of each component. 
        // singleton does not count as component
        void GetAllComponents(std::vector<std::vector<int> > &allComponents); 
        void GetComponent(const int &root, std::vector<int> &component, const int &componentID=-1); 

        UndirectedGraph(const int &N) : 
            _N(N)
            {
                _adjacencyMatrix.resize(N,N); 
                _belongComponent.resize(N); 
            }

}; 



#endif
