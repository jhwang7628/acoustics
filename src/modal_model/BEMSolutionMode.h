#ifndef BEM_SOLUTION_Mode_H 
#define BEM_SOLUTION_Mode_H 
#include <complex> 
#include <memory> 
#include <vector> 
#include <config.h> 
#include <geometry/TriangleMesh.hpp>

//##############################################################################
// This struct stores information needed for Kirchhoff integral. Might be used
// to compute multipole expansion as well. 
//
// Specifically, it stores BEM solution of pressure and velocity (prescribed
// boundary condition) on object surface. See p(y) and p_n(y) in Eq. 19 in 
// the paper: 2010, Changxi, Rigid-Body Fracture Sound with Precomputed Soundbanks. 
//
// If BEM solutions are from FastBEM, use class io/FBemReader to parse
// input/output. 
//##############################################################################
struct BEMSolutionMode
{
    private: 
        // this class is non-copyable. C++11 is needed for this syntax.  
        BEMSolutionMode(const BEMSolutionMode &other) = delete; 
        BEMSolutionMode &operator =(const BEMSolutionMode &other) = delete; 

    public: 
        BEMSolutionMode(){}
        std::shared_ptr<TriangleMesh<REAL> >    mesh; 
        std::vector<std::complex<REAL> >        pressures; 
        std::vector<std::complex<REAL> >        velocities; 
};

#endif
