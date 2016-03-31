#include <iostream> 

#include <utils/STL_Wrapper.h>
#include <tools/unit_testing/testing.h> 
#include <graph/UndirectedGraph.h> 


int main() {

    std::cout << "testing on undrected graph class \n"; 

    UndirectedGraph graph(40); 

    graph.AddEdgeToStage(0,1); 
    graph.AddEdgeToStage(0,2); 
    graph.AddEdgeToStage(1,2); 

    graph.AddEdgeToStage(3,4); 
    graph.AddEdgeToStage(4,5); 
    graph.AddEdgeToStage(5,6); 
    graph.AddEdgeToStage(4,6); 

    graph.AddEdgeToStage(8,9);



    //graph.AddEdgeToStage(40,5); 
    graph.AddEdgeToStage(39,1); 

    graph.ConstructEdges(); 

    Eigen::MatrixXd entries; 
    graph.GetAllNonzeroElements(entries); 


    COUT_SDUMP(entries); 

    std::vector<int> component; 
    graph.GetComponent(0, component); 
    STL_Wrapper::PrintVectorContent(std::cout, component); 

    graph.GetComponent(3, component); 
    STL_Wrapper::PrintVectorContent(std::cout, component); 


    std::vector<std::vector<int> > allComponents; 
    graph.GetAllComponents(allComponents); 

    std::cout << "there are " << allComponents.size() << " connected components found in the graph" << std::endl; 
    for (size_t ii=0; ii<allComponents.size(); ii++) 
    {
        std::cout << "component " << ii << " : "; 
        STL_Wrapper::PrintVectorContent(std::cout, allComponents[ii]); 
    }


    return 0; 
}
