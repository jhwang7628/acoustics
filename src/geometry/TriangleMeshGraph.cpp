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
    AdjListNode *newNode = new AdjListNode(dest); 
    newNode->next = array.at(src).head; 
    array.at(src).head = newNode; 
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
NeighboursOfTriangle(const int &t_id, std::vector<int> &neighbours) const
{
    assert(_graph_built); 
    neighbours.clear(); 
    AdjListNode *node = _graph.array.at(t_id).head; 
    while (node != nullptr)
    {
        neighbours.push_back(node->dest); 
        node = node->next; 
    }
}

//##############################################################################
//##############################################################################
template class TriangleMeshGraph<float>; 
template class TriangleMeshGraph<double>; 
