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
        using TriangleMesh<T>::m_vertices;
        using TriangleMesh<T>::m_triangles; 
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

        Graph _graph; 
        bool  _graph_built = false;

    public: 
        REAL ComputeClosestPointOnMesh(const int &startTriangleIndex, const Vector3d &queryPoint, Vector3d &closestPoint, int &closestTriangle, Vector3d &projectedPoint, const REAL &errorTol=0.99, const int &N_neighbours=10, const int &maxLevel=1) const; 
        void BuildGraph(); 
        void NeighboursOfTriangleRec(const int &t_id, const size_t &maxReach, std::set<int> &neighbours, std::set<int> &memo) const;
        void NeighboursOfTriangle(const int &t_id, const size_t &maxReach, std::set<int> &neighbours) const;
};

//##############################################################################
// Class TriangleDistanceComp
//##############################################################################
template <typename T>
class TriangleDistanceComp
{
    private: 
        const TriangleMeshGraph<T> *mesh; 
        const Vector3d &referencePoint; 
    public:
        TriangleDistanceComp(const TriangleMeshGraph<T> *m, const Vector3d &ref)
            : mesh(m), referencePoint(ref)
        {
            assert(mesh); 
        }
        bool operator()(const int &tid_ii, const int &tid_jj) const
        {
            const Vector3d v_ii = (mesh->TriangleCentroid(tid_ii) - referencePoint); 
            const Vector3d v_jj = (mesh->TriangleCentroid(tid_jj) - referencePoint); 
            const REAL d_ii = v_ii.lengthSqr(); 
            const REAL d_jj = v_jj.lengthSqr(); 
            return (d_ii < d_jj); 
        }
};


#endif
