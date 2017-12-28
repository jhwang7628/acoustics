#include <fstream>
#include <iostream> 
#include "igl/read_triangle_mesh.h"
#include "igl/write_triangle_mesh.h"
#include "Eigen/Dense" 

//##############################################################################
// Function Main
//##############################################################################
int main(int argc, char **argv)
{
#ifndef USE_IGL
    std::cerr << "**ERROR** Need libigl\n"; 
    return 2; 
#endif
    if (argc != 3)
    {
        std::cout << "**Usage: " << argv[0] 
                                 << " <input_mesh> <output_mesh> "
                                 << "\n"; 
        return 1; 
    }

    std::cout << "Reading input mesh ...\n"; 
    Eigen::MatrixXd V; 
    Eigen::MatrixXi F, Fnew; 
    igl::read_triangle_mesh(argv[1], V, F); 

    std::cout << "Triangulating ...\n"; 
    igl::polygon_mesh_to_triangle_mesh(F, Fnew);

    std::cout << "Writing output mesh ...\n";
    igl::write_triangle_mesh(argv[2], V, Fnew);

    return 0;
}
