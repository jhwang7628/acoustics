#include "SimplePointSet.hpp" 
#include <iostream> 

using namespace std; 

int main()
{
    Vertex<double> vert(0.0,0.0,0.0); 

    SimplePointSet<double> sps; 
    sps.loadPointSet("examples/cellCentroid_sphere_fluent5y0b.dat"); 
    for (Vertex<double> & v : sps.vlist) 
    {
        //cout << v << endl; 
    }

    cout << sps << endl;
    return 0; 
}
