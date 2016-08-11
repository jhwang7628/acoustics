#ifndef KIRCHHOFF_INTEGRAL_SOLVER_H 
#define KIRCHHOFF_INTEGRAL_SOLVER_H 
#include <config.h> 
#include <modal_model/BEMSolutionMode.h>
#include <geometry/TriangleMesh.hpp>

//##############################################################################
// This class evaluates the Kirchhoff integral. See Eq. 19 in the paper
// 2010 Changxi Zheng, Rigid-Body Fracture Sound with Precomputed Soundbanks.
//##############################################################################
class KirchhoffIntegralSolver
{
    private: 
        std::shared_ptr<TriangleMesh<REAL> >            _mesh; 
        std::vector<std::shared_ptr<BEMSolutionMode> >  _BEMSolutions; 

    public: 
        KirchhoffIntegralSolver(){}
        KirchhoffIntegralSolver(std::shared_ptr<TriangleMesh<REAL> > &mesh)
            : _mesh(mesh)
        {}
        
        inline void SetMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh){_mesh = mesh;}
        void ReadFromFBem(const std::string &FBemConfigFile, const std::string &FBemOutputFile);
};

#endif
