#include <fstream>
#include <iostream> 
#include "igl/read_triangle_mesh.h"
#include "igl/write_triangle_mesh.h"
#include "igl/remove_duplicates.h"
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
    Eigen::MatrixXd V, Vnew; 
    Eigen::MatrixXi F, Fnew; 
    Eigen::VectorXi I; 
    igl::read_triangle_mesh(argv[1], V, F); 

    std::cout << "Removing duplicates ...\n"; 
    igl::remove_duplicates(V, F, Vnew, Fnew, I, 0.000057);

    std::cout << "Writing output mesh ...\n";
    igl::write_triangle_mesh(argv[2], Vnew, Fnew);

    return 0;
}
