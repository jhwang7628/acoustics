#ifndef FBEM_READER_H 
#define FBEM_READER_H 
#include <string> 
#include <modal_model/BEMSolutionMode.h> 
#include <geometry/TriangleMesh.hpp> 

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
    public: 
        bool CheckFBemInputAgainstMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh, const std::string &solutionFile); 
        bool ReadFBemOutputToInfo(std::shared_ptr<BEMSolutionMode> &solution, const std::string &solutionFile); 
};

#endif
