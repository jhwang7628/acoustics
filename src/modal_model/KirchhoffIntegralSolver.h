#ifndef KIRCHHOFF_INTEGRAL_SOLVER_H 
#define KIRCHHOFF_INTEGRAL_SOLVER_H 
#include <complex>
#include <config.h> 
#include <modal_model/BEMSolutionMode.h>
#include <geometry/TriangleMesh.hpp>

//##############################################################################
// This class evaluates the Kirchhoff integral. See Eq. 19 in the paper
// 2010 Changxi Zheng, Rigid-Body Fracture Sound with Precomputed Soundbanks.
//##############################################################################
class KirchhoffIntegralSolver
{
    public: 
        static const std::complex<REAL> j; 
    private: 
        std::shared_ptr<TriangleMesh<REAL> >            _mesh; 
        std::vector<std::shared_ptr<BEMSolutionMode> >  _BEMSolutions; 

    public: 
        KirchhoffIntegralSolver(){}
        KirchhoffIntegralSolver(std::shared_ptr<TriangleMesh<REAL> > &mesh)
            : _mesh(mesh)
        {}
        
        inline void SetMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh){_mesh = mesh;}

        std::complex<REAL> Evaluate_G(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const REAL &rCached = -1.0); 
        std::complex<REAL> Evaluate_dG_dn(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const Vector3d &normal, const REAL &rCached = -1.0); 
        void AddFBemSolution(const std::string &fBemConfigFile, const std::string &fBemOutputFile);
};

#endif
