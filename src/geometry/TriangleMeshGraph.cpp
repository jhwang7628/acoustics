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
// Function Save
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::Graph::
Save(const std::string &filename)
{
    assert(array.size()>0);
    // test if file exists, if it does, abort
    {
        ifstream f(filename.c_str());
        if (f.good()) return; 
        f.close(); 
    }
    // save to file
    ofstream stream(filename.c_str(), std::ios::out|std::ios::binary);
    int N_verts = numNodes, N_vert_neighbours; 
    stream.write((char*) &(N_verts), sizeof(int)); 
    for (int ii=0; ii<N_verts; ++ii) 
    {
        N_vert_neighbours = array.at(ii).size(); 
        stream.write((char*) &(N_vert_neighbours), sizeof(int)); 
        stream.write((char*) &(array.at(ii)[0]), sizeof(int)*N_vert_neighbours); 
    }
    stream.close(); 
}

//##############################################################################
// Function Load
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::Graph::
Load(const std::string &filename)
{
    // test if file exists, if it does not, throw exception
    {
        ifstream f(filename.c_str());
        if (!f.good()) 
            throw std::runtime_error("**ERROR** load file does not exist");  
        f.close(); 
    }
    // save to file
    array.clear(); 
    ifstream stream(filename.c_str(), std::ios::in|std::ios::binary);
    int N_verts, N_vert_neighbours; 
    stream.read((char*) &(N_verts), sizeof(int)); 
    numNodes = N_verts; 
    array.resize(numNodes);
    for (int ii=0; ii<numNodes; ++ii) 
    {
        std::vector<int> &list = array.at(ii); 
        stream.read((char*) &(N_vert_neighbours), sizeof(int)); 
        list.resize(N_vert_neighbours); 
        stream.read((char*) &(list[0]), sizeof(int)*N_vert_neighbours); 
    }
    stream.close(); 
}

//##############################################################################
// Function TryLoad
//##############################################################################
template <typename T> 
bool TriangleMeshGraph<T>::Graph::
TryLoad(const std::string &filename)
{
    try 
    {
        Load(filename); 
    }
    catch (std::runtime_error &e)
    {
        return false; 
    }
    return true; 
}

//##############################################################################
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
FindKNearestTrianglesGraph(const int &k_max, const Vector3<T> &point, const int &maxLevel, const int &startTriangleIndex, std::vector<int> &triangleIndices) const
{
    const TriangleDistanceComp<T> sorter(this, point); 
    if (maxLevel <= 1)
    {
        NeighboursOfTriangle(startTriangleIndex, triangleIndices); 
    }
    else 
    {
        std::set<int> neighbours; 
        NeighboursOfTriangle(startTriangleIndex, maxLevel, neighbours); 
        triangleIndices = std::vector<int>(neighbours.begin(), neighbours.end()); 
    }
    const int k = std::min<int>(k_max, triangleIndices.size());
    std::partial_sort(triangleIndices.begin(), triangleIndices.begin()+k, triangleIndices.end(), sorter); 
    triangleIndices = std::vector<int>(triangleIndices.begin(), triangleIndices.begin()+k); 
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
BuildGraph(const std::string &loadName, const T &nnRadius)
{
#ifdef USE_BOOST
    boost::timer::auto_cpu_timer timer("Boost timer: Building Graph for mesh takes %w seconds\n" );
#endif
    typedef std::set<int>::iterator SetIterator;
    assert(m_vertices.size()>0 && m_triangles.size()>0); 
    _graph.Initialize(m_triangles.size()); 

    const bool loaded = _graph.TryLoad(loadName);
    if (loaded) 
    {
        std::cout << "loaded name: " << loadName << std::endl;
        _graph_built = true; 
        return; 
    }

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

    if (!loaded) 
    {
        std::cout << "saved name: " << loadName << std::endl;
        _graph.Save(loadName); 
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
    neighbours.insert(_graph.array.at(t_id).begin(), _graph.array.at(t_id).end()); 
    memo.insert(t_id);
    // recursively call all neighbours
    std::set<int> newNeighbours; 
    if (maxReach-1==0)
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
// Function NeighboursOfTriangle
//   This overloading is more efficient if maxReach==1 in the above function.
//##############################################################################
template <typename T> 
void TriangleMeshGraph<T>::
NeighboursOfTriangle(const int &t_id, std::vector<int> &neighbours) const
{
    assert(_graph_built); 
    neighbours = _graph.array.at(t_id); 
}

//##############################################################################
//##############################################################################
template class TriangleMeshGraph<float>; 
template class TriangleMeshGraph<double>; 
