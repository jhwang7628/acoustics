#ifndef TRIANGLE_MESH_GRAPH_HPP
#define TRIANGLE_MESH_GRAPH_HPP 
#include <vector>
#include "geometry/TriangleMeshKDTree.hpp"
#include "utils/timer.hpp"

//##############################################################################
// Class TriangleDistanceComp
//##############################################################################
template <typename T>
class TriangleDistanceComp
{
    private: 
        const TriangleMeshKDTree<T> *mesh; 
        const Vector3<T> &referencePoint; 
    public:
        TriangleDistanceComp(const TriangleMeshKDTree<T> *m, const Vector3<T> &ref)
            : mesh(m), referencePoint(ref)
        {
            assert(mesh); 
        }
        bool operator()(const int &tid_ii, const int &tid_jj) const
        {
            assert(mesh);
            const Vector3<T> v_ii = (mesh->TriangleCentroid(tid_ii) - referencePoint); 
            const Vector3<T> v_jj = (mesh->TriangleCentroid(tid_jj) - referencePoint); 
            const T d_ii = v_ii.lengthSqr(); 
            const T d_jj = v_jj.lengthSqr(); 
            return (d_ii < d_jj); 
        }
};

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
        struct Graph
        {
            int numNodes; 
            std::vector<std::vector<int> > array; 
            void Initialize(const int &V); 
            void AddEdge(const int &src, const int &dest); 
        };

        Graph _graph; 
        bool  _graph_built = false;

    public: 
        static std::vector<Timer<false> > timers;
        void FindKNearestTrianglesGraph(const int &k, const Vector3<T> &point, const int &maxLevel, const int &startTriangleIndex,  std::vector<int> &triangleIndices) const;
        T ComputeClosestPointOnMesh(const int &startTriangleIndex, const Vector3<T> &queryPoint, Vector3<T> &closestPoint, int &closestTriangle, Vector3<T> &projectedPoint, const T &errorTol=0.999, const int &N_neighbours=5, const int &maxLevel=1) const; 
        void BuildGraph(const T &nnRadius=0.0); 
        void NeighboursOfTriangleRec(const int &t_id, const size_t &maxReach, std::set<int> &neighbours, std::set<int> &memo) const;
        void NeighboursOfTriangle(const int &t_id, const size_t &maxReach, std::set<int> &neighbours) const;
};

#endif
