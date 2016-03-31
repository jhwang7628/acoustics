#include <graph/UndirectedGraph.h>
#include <utils/STL_Wrapper.h> 


void UndirectedGraph::AddEdgeToStage(const int &x_i, const int &x_j)
{
    if (x_i >= _N || x_j >= _N) throw std::runtime_error("**ERROR** requested node out-of-bounds"); 
    _stagedEdges.push_back(TripletInt(x_i, x_j, 1)); 
    _stagedEdges.push_back(TripletInt(x_j, x_i, 1)); 
}

void UndirectedGraph::ConstructEdges(const bool &destroyOld)
{
    if (!destroyOld) 
    {
        for (int col=0; col<_adjacencyMatrix.outerSize(); col++) 
            for (SparseMatrixXi::InnerIterator it(_adjacencyMatrix,col); it; ++it) 
                _stagedEdges.push_back(TripletInt(it.row(), it.col(), it.value())); 
    }

    const int addEdges = _stagedEdges.size();
    if (addEdges==0) return; 

    _adjacencyMatrix.setFromTriplets(_stagedEdges.begin(), _stagedEdges.end()); 
    _stagedEdges.clear(); 

    _M = addEdges; 

}

void UndirectedGraph::GetAllNonzeroElements(Eigen::MatrixXd &nonzeroElements)
{

    nonzeroElements.resize(_adjacencyMatrix.nonZeros(), 3); 

    std::cout << _adjacencyMatrix.nonZeros() << std::endl; 

    int count = 0;
    for (int col=0; col<_adjacencyMatrix.outerSize(); col++) 
        for (SparseMatrixXi::InnerIterator it(_adjacencyMatrix,col); it; ++it) 
        {
            nonzeroElements(count,0) = it.row(); 
            nonzeroElements(count,1) = it.col(); 
            nonzeroElements(count,2) = it.value(); 

            count ++; 
        } 
}

void UndirectedGraph::GetAllComponents(std::vector<std::vector<int> > &allComponents)
{
    _belongComponent.clear(); 
    _belongComponent.resize(_N,-1); 

    allComponents.clear(); 

    int root; 

    for (int col=0; col<_adjacencyMatrix.outerSize(); col++) 
    {
        if (_adjacencyMatrix.col(col).nonZeros() != 0)
        {
            // use the first element in the column as the root for BFS 
            SparseMatrixXi::InnerIterator it(_adjacencyMatrix,col); 
            root = it.row(); 

            if (_belongComponent[root] == -1)  // not searched yet
            {
                std::vector<int> component; 
                GetComponent(root, component, allComponents.size()); 
                allComponents.push_back(component); 
            } 
        }
    }
}


// BFS graph traversal 
// Ref: https://en.wikipedia.org/wiki/Breadth-first_search
void UndirectedGraph::GetComponent(const int &root, std::vector<int> &component, const int &componentID)
{

    if (root >= _N) throw std::runtime_error("**ERROR** given root is larger than number of nodes in the graph"); 

    component.clear();

    std::vector<int> distances(_N, std::numeric_limits<int>::max()); 
    std::vector<int>   parents(_N, -1); 

    distances[root] = 0;
    component.push_back(root); 
    _belongComponent[root] = componentID; 
    std::queue<int> searchQueue; 
    searchQueue.push(root); 

    int current; 
    int currentNeighbour; 

    while(!searchQueue.empty())
    {
        current = searchQueue.front(); 
        searchQueue.pop(); 
        for (SparseMatrixXi::InnerIterator it(_adjacencyMatrix,current); it; ++it)
        {
            currentNeighbour = it.row(); 
            if (distances[currentNeighbour] == std::numeric_limits<int>::max()) // not visited yet
            { 

                _belongComponent[currentNeighbour] = componentID; 

                component.push_back(currentNeighbour); 
                distances[currentNeighbour] = distances[current] + 1; 
                parents[currentNeighbour] = current; 
                searchQueue.push(currentNeighbour); 
            }
        }
        
    }

}

void UndirectedGraph::GetAllComponentsVerbose(std::vector<std::vector<int> > &allComponents)
{
    GetAllComponents(allComponents); 

    std::cout << "there are " << allComponents.size() << " connected components found in the graph" << std::endl; 
    for (size_t ii=0; ii<allComponents.size(); ii++) 
    {
        std::cout << "component " << ii << " : "; 
        STL_Wrapper::PrintVectorContent(std::cout, allComponents[ii]); 
    }



}

