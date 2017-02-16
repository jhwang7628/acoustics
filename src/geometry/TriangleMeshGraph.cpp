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
// Function ComputeClosestPointOnMesh
//   Wrapper function using graph local search (and KD-tree as fall-back)
//##############################################################################
template <typename T> 
REAL TriangleMeshGraph<T>::
ComputeClosestPointOnMesh(const int &startTriangleIndex, const Vector3d &queryPoint, Vector3d &closestPoint, int &closestTriangle, Vector3d &projectedPoint, const REAL &errorTol, const int &N_neighbours) const
{
    const int maxLevel = 1; 
    int level = 1;
    std::set<int> neighbours; 
    std::vector<int> triangleIndices;
    REAL distance; 
    {
        boost::timer::auto_cpu_timer t("search for neighbours takes %w secs\n"); 
    while (neighbours.size() < N_neighbours && level<=maxLevel) 
    {
        NeighboursOfTriangle(startTriangleIndex, level, neighbours); 
        ++level; 
    }
    }
    {
        boost::timer::auto_cpu_timer t("sort/trim vector and helper takes %w secs\n"); 
    triangleIndices = std::vector<int>(neighbours.begin(), neighbours.end()); 
    TriangleDistanceComp<T> sorter(this, queryPoint); 
    std::sort(triangleIndices.begin(), triangleIndices.end(), sorter); 
    const int len = std::min((int)triangleIndices.size(), N_neighbours); 
    if (len<1) 
        throw std::runtime_error("**ERROR** no neighbours found"); 
    triangleIndices = std::vector<int>(triangleIndices.begin(), triangleIndices.begin()+len); 
    distance = this->ComputeClosestPointOnMeshHelper(queryPoint, triangleIndices, closestPoint, closestTriangle, projectedPoint);
    }

    // fall back to KNN
    const Vector3d triNormal = this->triangle_normal(closestTriangle).normalized(); 
    const Vector3d computedNormal = (closestPoint-queryPoint).normalized(); 
    const REAL error = abs(triNormal.dotProduct(computedNormal)); 
    if (error < errorTol)
    {
        std::cout << ">>>> fall back\n";
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
// Function NeighboursOfTriangle
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
NeighboursOfTriangle(const int &t_id, const size_t &maxStride, std::set<int> &neighbours) const
{
    assert(_graph_built); 
    if (maxStride==0) 
        return; 
    AdjListNode *node = _graph.array.at(t_id).head; 
    while (node != nullptr)
    {
        neighbours.insert(node->dest); 
        node = node->next; 
    }
    // recursively call all neighbours
    for (const int &n : neighbours) 
        NeighboursOfTriangle(n, maxStride-1, neighbours); 
}

//##############################################################################
//##############################################################################
template class TriangleMeshGraph<float>; 
template class TriangleMeshGraph<double>; 
