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
    if (argc < 4)
    {
        std::cout << "**Usage: " << argv[0] 
                                 << " <input_mesh> <marker_sphere_mesh> "
                                 << " <output_vertices_file> " 
                                 << " [box_size_x] [box_size_y] [box_size_z] "
                                 << " [output_box_file] "
                                 << "\n"; 
        return 1; 
    }

    std::cout << "Reading input mesh and marker ...\n"; 
    Eigen::MatrixXd V[2]; 
    Eigen::MatrixXi F[2]; 
    igl::read_triangle_mesh(argv[1], V[0], F[0]); 
    igl::read_triangle_mesh(argv[2], V[1], F[1]); 

    Eigen::Vector3d marker_centroid = V[1].colwise().sum()/(double)V[1].rows(); 
    Eigen::Vector3d marker_box[2], diff; 
    if (argc >= 7)
    {
        diff << atof(argv[4]), atof(argv[5]), atof(argv[6]);
    }
    else 
    {
        const Eigen::Vector3d vmin = V[1].colwise().minCoeff(); 
        const Eigen::Vector3d vmax = V[1].colwise().maxCoeff(); 
        const Eigen::Vector3d diff = vmax - vmin; 
    }
    marker_box[0] = marker_centroid - diff/2.0; 
    marker_box[1] = marker_centroid + diff/2.0; 

    std::cout << "Filter vertices ...\n"; 
    const auto InsideBox = 
        [&](const Eigen::Vector3d v)->bool {
                for (int dd=0; dd<3; ++dd) 
                    if (v[dd] > marker_box[1][dd] ||
                        v[dd] < marker_box[0][dd])
                        return false;
                return true;};
    std::vector<int> keep; 
    for (int ii=0; ii<V[0].rows(); ++ii)
        if (InsideBox(V[0].row(ii)))
            keep.push_back(ii); 

    std::ofstream stream(argv[3]); 
    for (const int &k : keep)
        stream << k << "\n"; 
    stream.close();

    if (argc >= 8) 
    {
        Eigen::MatrixXd Vb(8,3); 
        Eigen::MatrixXi Fb(6,4); 
        Vb << marker_box[0][0], marker_box[0][1], marker_box[0][2],
              marker_box[0][0], marker_box[1][1], marker_box[0][2],
              marker_box[0][0], marker_box[1][1], marker_box[1][2],
              marker_box[0][0], marker_box[0][1], marker_box[1][2],
              marker_box[1][0], marker_box[0][1], marker_box[0][2],
              marker_box[1][0], marker_box[1][1], marker_box[0][2],
              marker_box[1][0], marker_box[1][1], marker_box[1][2],
              marker_box[1][0], marker_box[0][1], marker_box[1][2]; 
        Fb << 0, 3, 2, 1,
              4, 7, 3, 0,
              4, 5, 6, 7,
              1, 2, 6, 5,
              7, 6, 2, 3,
              0, 1, 5, 4;
        igl::write_triangle_mesh(argv[7], Vb, Fb); 
    }

    return 0;
}
