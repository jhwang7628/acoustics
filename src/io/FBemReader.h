#ifndef FBEM_READER_H 
#define FBEM_READER_H 
#include <string> 
#include <config.h>
#include <geometry/Point3.hpp> 
#include <geometry/TriangleMesh.hpp> 
#include <linearalgebra/Tuple3.hpp> 
#include <modal_model/BEMSolutionMode.h> 

//##############################################################################
// This class handles all IO operations for FBem results.
// An example folder of FBem results can be found in google drive. 
//
// List of files this can handle:
//  - tet-*_*.txt: BEM solution file for the triangle mesh.
//
// Notes:
//  1) Segments of code might be grabbed from modec repo, specifically from files:
//      modec/src/fast_fit_multipole.cpp
//##############################################################################
class FBemReader
{
    private:
        const REAL _meshCheckTolerance = 1E-4; 
        const bool _checkAllowFlipNormal = true; // FBem mesh seems to have normals flipped, do a flip back when checking
        bool _ReadFBemInputToGeometry(const char *fBemInputFile, std::vector<Point3d> &vertices, std::vector<Tuple3ui> &triangles); 

    public: 
        bool CheckFBemInputAgainstMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh, const std::string &fBemInputFile); 
        bool ReadFBemOutputToInfo(std::shared_ptr<BEMSolutionMode> &solution, const std::string &fBemOutputFile); 
};

#endif
