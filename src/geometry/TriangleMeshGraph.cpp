#include "geometry/TriangleMeshGraph.hpp" 

//##############################################################################a
// Destructor
//##############################################################################a
template <typename T> 
TriangleMeshGraph<T>::AdjList::
~AdjList()
{
    AdjListNode *toDelete = head; 
    while (toDelete != nullptr)
    {
        AdjListNode *next = toDelete->next; 
        delete toDelete; 
        toDelete = next; 
    }
}

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
    AdjListNode *newNode = new AdjListNode(dest); 
    newNode->next = array.at(src).head; 
    array.at(src).head = newNode; 

    // backward edge
    newNode = new AdjListNode(src); 
    newNode->next = array.at(dest).head; 
    array.at(dest).head = newNode; 
}

//##############################################################################
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
FindKNearestTrianglesGraph(const int &k, const Vector3d &point, const int &maxLevel, const int &startTriangleIndex, std::vector<int> &triangleIndices) const
{
    std::set<int> neighbours; 
    NeighboursOfTriangle(startTriangleIndex, maxLevel, neighbours); 
    triangleIndices = std::vector<int>(neighbours.begin(), neighbours.end()); 
    TriangleDistanceComp<T> sorter(this, point); 
    std::sort(triangleIndices.begin(), triangleIndices.end(), sorter); 
    if (triangleIndices.size()<k)
        throw std::runtime_error("**ERROR** not enough neighbours in the graph"); 
    triangleIndices = std::vector<int>(triangleIndices.begin(), triangleIndices.begin()+k); 
}

//##############################################################################
// Function ComputeClosestPointOnMesh
//   Wrapper function using graph local search (and KD-tree as fall-back)
//##############################################################################
template <typename T> 
REAL TriangleMeshGraph<T>::
ComputeClosestPointOnMesh(const int &startTriangleIndex, const Vector3d &queryPoint, Vector3d &closestPoint, int &closestTriangle, Vector3d &projectedPoint, const REAL &errorTol, const int &N_neighbours, const int &maxLevel) const
{
    // find NN using graph and compute point
    std::vector<int> triangleIndices; 
    REAL distance;
    FindKNearestTrianglesGraph(N_neighbours, queryPoint, maxLevel, startTriangleIndex, triangleIndices); 
    distance = this->ComputeClosestPointOnMeshHelper(queryPoint, triangleIndices, closestPoint, closestTriangle, projectedPoint);

    // fall back to KNN using KD-tree
    const Vector3d triNormal = this->triangle_normal(closestTriangle).normalized(); 
    const Vector3d computedNormal = (closestPoint-queryPoint).normalized(); 
    const REAL error = abs(triNormal.dotProduct(computedNormal)); 
    if (error < errorTol) // FIXME debug
    //int graph = closestTriangle; 
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
BuildGraph()
{
#ifdef USE_BOOST
    boost::timer::auto_cpu_timer timer("Boost timer: Building Graph for mesh takes %w seconds\n" );
#endif
    typedef std::set<int>::const_iterator SetIterator;
    assert(m_vertices.size()>0 && m_triangles.size()>0); 
    _graph.Initialize(m_triangles.size()); 

    // get vertex neighbours and assign them to triangles
    std::vector<std::set<int> > vertexNeighbors; 
    this->get_vtx_tgls(vertexNeighbors); 
    const int N_faces = m_triangles.size(); 
    for (int ii=0; ii<N_faces; ++ii)
    {
        const Tuple3ui &t_id = this->triangle_ids(ii); 
        std::set<int> triangleNeighbours; 
        for (int jj=0; jj<3; ++jj)
        {
            const std::set<int> &neighbours = vertexNeighbors.at(t_id[jj]); 
            triangleNeighbours.insert(neighbours.begin(), neighbours.end()); 
        }
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
    AdjListNode *node = _graph.array.at(t_id).head; 
    while (node != nullptr)
    {
        neighbours.insert(node->dest); 
        node = node->next; 
    }
    memo.insert(t_id);  // FIXME debug
    // recursively call all neighbours
    std::set<int> newNeighbours; 
    for (const int &n : neighbours) 
        NeighboursOfTriangleRec(n, maxReach-1, newNeighbours, memo); 
    neighbours.insert(newNeighbours.begin(), newNeighbours.end());
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
