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
    if (argc != 2)
    {
        std::cout << "**Usage: " << argv[0]
                                 << " <input_mesh>"
                                 << "\n";
        return 1;
    }

    std::cout << "Reading input mesh ...\n";
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(argv[1], V, F);

    std::cout << "Finding Obj centroids ...\n";
    Eigen::Vector3d c = V.colwise().sum()/(double)V.rows();
    std::cout << "  centroid = " << std::setprecision(16) << std::fixed
              << c.transpose() << std::endl;

    return 0;
}
