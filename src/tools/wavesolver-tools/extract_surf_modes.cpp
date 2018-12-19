#include <iostream>
#include <fstream>
#include <cassert>
#include "geometry/TetMeshIndexToSurfaceMesh.h"
#include "deformable/ModeData.h"
int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "**Usage: " << argv[0]
                  << " <geo_file> <tet_modes> <surf_modes>\n";
        return -1;
    }
    TetMeshIndexToSurfaceMesh table;
    table.ReadFromGeoFile(argv[1]);
    ModeData volmMode, surfMode;
    volmMode.read(argv[2]);
    // copy eigenvalues
    surfMode._omegaSquared = volmMode._omegaSquared;
    // filter and copy eigenvectors
    int NV = table.N_surfaceVertices();
    for (int m=0; m<volmMode.numModes(); ++m) {
        std::vector<double> mode(NV*3);
        for (int v_s=0; v_s<NV; ++v_s) {
            const int v_t = table.GetTetIndex(v_s);
            mode.at(v_s*3+0) = volmMode.mode(m).at(v_t*3+0);
            mode.at(v_s*3+1) = volmMode.mode(m).at(v_t*3+1);
            mode.at(v_s*3+2) = volmMode.mode(m).at(v_t*3+2);
        }
        surfMode._modes.push_back(mode);
    }
    surfMode.write(argv[3]);

    return 0;
}
