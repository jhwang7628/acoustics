#include "geometry/TriangleMeshGraph.hpp" 
#include "boost/timer/timer.hpp"
//##############################################################################
// Static variable initialize
//##############################################################################
template <typename T> 
std::vector<Timer<false> > TriangleMeshGraph<T>::timers(20);

//##############################################################################a
// Function Initialize
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::Graph::
Initialize(const int &V)
{
    numNodes = V; 
    array.resize(V); 
}

//##############################################################################
// Function AddEdge
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::Graph::
AddEdge(const int &src, const int &dest)
{
    // forward edge
    {
        std::vector<int> &list = array.at(src); 
        if (std::find(list.begin(),list.end(),dest) == list.end())
            list.push_back(dest); 
    }

    // backward edge
    {
        std::vector<int> &list = array.at(dest); 
        if (std::find(list.begin(),list.end(),src) == list.end())
            list.push_back(src);
    }
}

//##############################################################################
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
FindKNearestTrianglesGraph(const int &k, const Vector3<T> &point, const int &maxLevel, const int &startTriangleIndex, std::vector<int> &triangleIndices) const
{
    TriangleMeshGraph<T>::timers[0].start();
    std::set<int> neighbours; 
    NeighboursOfTriangle(startTriangleIndex, maxLevel, neighbours); 
    TriangleMeshGraph<T>::timers[0].pause();
    TriangleMeshGraph<T>::timers[1].start();
    triangleIndices = std::vector<int>(neighbours.begin(), neighbours.end()); 
    TriangleDistanceComp<T> sorter(this, point); 
    std::sort(triangleIndices.begin(), triangleIndices.end(), sorter); 
    TriangleMeshGraph<T>::timers[1].pause();
    if (triangleIndices.size()<k)
        throw std::runtime_error("**ERROR** not enough neighbours in the graph"); 
    TriangleMeshGraph<T>::timers[2].start();
    triangleIndices = std::vector<int>(triangleIndices.begin(), triangleIndices.begin()+k); 
    TriangleMeshGraph<T>::timers[2].pause();
}

//##############################################################################
// Function ComputeClosestPointOnMesh
//   Wrapper function using graph local search (and KD-tree as fall-back)
//##############################################################################
template <typename T> 
T TriangleMeshGraph<T>::
ComputeClosestPointOnMesh(const int &startTriangleIndex, const Vector3<T> &queryPoint, Vector3<T> &closestPoint, int &closestTriangle, Vector3<T> &projectedPoint, const T &errorTol, const int &N_neighbours, const int &maxLevel) const
{
    // find NN using graph and compute point
    std::vector<int> triangleIndices; 
    T distance;
    FindKNearestTrianglesGraph(N_neighbours, queryPoint, maxLevel, startTriangleIndex, triangleIndices);
    distance = this->ComputeClosestPointOnMeshHelper(queryPoint, triangleIndices, closestPoint, closestTriangle, projectedPoint);

    // fall back to KNN using KD-tree
    const Vector3<T> triNormal = this->triangle_normal(closestTriangle).normalized(); 
    const Vector3<T> computedNormal = (closestPoint-queryPoint).normalized(); 
    const T error = abs(triNormal.dotProduct(computedNormal)); 
    if (error < errorTol)
    {
        triangleIndices.clear(); 
#ifdef USE_OPENMP
#pragma omp critical
#endif
        this->FindKNearestTriangles(N_neighbours, queryPoint, triangleIndices); 
        distance = this->ComputeClosestPointOnMeshHelper(queryPoint, triangleIndices, closestPoint, closestTriangle, projectedPoint);
    }
    return distance; 
}

//##############################################################################
// Function BuildGraph
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
BuildGraph(const T &nnRadius)
{
#ifdef USE_BOOST
    boost::timer::auto_cpu_timer timer("Boost timer: Building Graph for mesh takes %w seconds\n" );
#endif
    typedef std::set<int>::iterator SetIterator;
    assert(m_vertices.size()>0 && m_triangles.size()>0); 
    _graph.Initialize(m_triangles.size()); 

    // grab vertex neighbours
    std::vector<std::set<int> > vertexNeighbors; 
    this->get_vtx_tgls(vertexNeighbors); 
    const int N_faces = m_triangles.size(); 
    for (int ii=0; ii<N_faces; ++ii)
    {
        const Tuple3ui &t_id = this->triangle_ids(ii); 
        std::set<int> triangleNeighbours; 
        // add all vertex neighbours
        for (int jj=0; jj<3; ++jj)
        {
            const std::set<int> &neighbours = vertexNeighbors.at(t_id[jj]); 
            triangleNeighbours.insert(neighbours.begin(), neighbours.end()); 
        }
        // add geometric neighbours (who might not be topologic neighbours)
        if (nnRadius>0 && this->_nnForest)
        {
            std::set<int> geometricNeighbours; 
            this->FindTrianglesWithinBall(this->TriangleCentroid(ii), nnRadius, geometricNeighbours); 
            triangleNeighbours.insert(geometricNeighbours.begin(), geometricNeighbours.end()); 
        }
        // add edge to graph using the set
        for (SetIterator it=triangleNeighbours.begin(); it!=triangleNeighbours.end(); ++it)
        {
            _graph.AddEdge(ii, *it); 
        }
    }
    _graph_built = true; 
}

//##############################################################################
// Function NeighboursOfTriangleRec
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
NeighboursOfTriangleRec(const int &t_id, const size_t &maxReach, std::set<int> &neighbours, std::set<int> &memo) const
{
    // base case if run out of reach, or memoized
    if (maxReach==0 || (memo.find(t_id)!=memo.end())) 
        return; 
    TriangleMeshGraph<T>::timers[3].start();
    neighbours.insert(_graph.array.at(t_id).begin(), _graph.array.at(t_id).end()); 
    TriangleMeshGraph<T>::timers[3].pause();
    TriangleMeshGraph<T>::timers[4].start();
    memo.insert(t_id);
    TriangleMeshGraph<T>::timers[4].pause();
    // recursively call all neighbours
    std::set<int> newNeighbours; 
    for (const int &n : neighbours) 
        NeighboursOfTriangleRec(n, maxReach-1, newNeighbours, memo); 
    TriangleMeshGraph<T>::timers[5].start();
    neighbours.insert(newNeighbours.begin(), newNeighbours.end());
    TriangleMeshGraph<T>::timers[5].pause();
}

//##############################################################################
// Function NeighboursOfTriangle
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
NeighboursOfTriangle(const int &t_id, const size_t &maxReach, std::set<int> &neighbours) const
{
    assert(_graph_built); 
    neighbours.clear(); 
    std::set<int> memo; 
    NeighboursOfTriangleRec(t_id, maxReach, neighbours, memo); 
}

//##############################################################################
//##############################################################################
template class TriangleMeshGraph<float>; 
template class TriangleMeshGraph<double>; 
