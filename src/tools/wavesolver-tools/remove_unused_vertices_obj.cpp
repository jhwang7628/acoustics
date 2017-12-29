#include <fstream>
#include <iostream> 
#include <map>
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
    Eigen::MatrixXd V, Vnew; 
    Eigen::MatrixXi F, Fnew; 
    igl::read_triangle_mesh(argv[1], V, F); 

    std::cout << "Remap vertices and remove unused ...\n"; 
    std::map<int,int> vmap; 
    int idx = 0;
    Fnew.resize(F.rows(), F.cols()); 
    for (int ii=0; ii<F.rows(); ++ii) 
        for (int jj=0; jj<F.cols(); ++jj)
        {
            if (vmap.find(F(ii,jj)) == vmap.end())
            {
                vmap[F(ii,jj)] = idx; 
                ++ idx; 
            }
            Fnew(ii,jj) = vmap.at(F(ii,jj)); 
        }

    Vnew.resize(idx+1, 3); 
    for (const auto &m : vmap)
        Vnew.row(m.second) = V.row(m.first); 

    std::cout << "Writing output mesh ...\n"; 
    igl::write_triangle_mesh(argv[2], Vnew, Fnew); 

    return 0; 
}
