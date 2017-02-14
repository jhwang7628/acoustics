#ifndef TRIANGLE_MESH_GRAPH_HPP
#define TRIANGLE_MESH_GRAPH_HPP 
#include <vector>
#include "geometry/TriangleMeshKDTree.hpp"

//##############################################################################
// Class TriangleMeshGraph
//   This class uses an adjacency-list based graph to represent adjacency 
//   structure of the mesh.
//   Reference: http://www.geeksforgeeks.org/graph-and-its-representations
//##############################################################################
template <typename T> 
class TriangleMeshGraph : public TriangleMeshKDTree<T> 
{
    protected: 
        struct AdjListNode
        {
            int dest; 
            AdjListNode *next; 
            AdjListNode(const int &d) 
                : dest(d), next(nullptr)
            {}
        };
        struct AdjList
        {
            AdjListNode *head; 
            AdjList()
                : head(nullptr)
            {}
            ~AdjList();
        };
        struct Graph
        {
            int numNodes; 
            std::vector<AdjList> array; 
            void Initialize(const int &V); 
            void AddEdge(const int &src, const int &dest); 
        };

        typedef typename TriangleMesh<T>::NeighborRec NeighborRec; 
        using TriangleMesh<T>::m_vertices;
        using TriangleMesh<T>::m_triangles; 
        Graph _graph; 
        bool  _graph_built = false;

    public: 
        void BuildGraph(); 
        void NeighboursOfTriangle(const int &t_id, const size_t &maxStride, std::set<int> &neighbours) const;
};

#endif
