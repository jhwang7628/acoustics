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
        REAL                                            _distanceLowBound = 1E-3; // won't be evaluated if distance lower than this value.
        std::shared_ptr<TriangleMesh<REAL> >            _mesh; 
        std::vector<std::shared_ptr<BEMSolutionMode> >  _BEMSolutions; 

        // for debug and testing
        REAL                                            _defaultTestEvaluateRadius = 1E-4;
        REAL                                            _defaultTestErrorTolerance = 1E-6;

    public: 
        KirchhoffIntegralSolver(){}
        KirchhoffIntegralSolver(std::shared_ptr<TriangleMesh<REAL> > &mesh)
            : _mesh(mesh)
        {}

        inline int N_Modes(){return _BEMSolutions.size();}
        inline REAL GetMode_omega(const int &mode){return _BEMSolutions.at(mode)->omega;}
        inline REAL GetMode_k(const int &mode){return _BEMSolutions.at(mode)->omega / _BEMSolutions.at(mode)->soundSpeed;}
        inline void SetMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh){_mesh = mesh;}
        inline void ClearSolutions(){_BEMSolutions.clear();}

        void AddFBemSolution(const std::string &fBemConfigFile, const std::string &fBemOutputFile, const REAL &omega);
        std::complex<REAL> Evaluate_G(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const REAL &rCached = -1.0) const; 
        std::complex<REAL> Evaluate_dG_dn(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const Vector3d &normal, const REAL &rCached = -1.0) const; 
        std::complex<REAL> Solve(const int &modeIndex, const Vector3d &listeningPoint) const; 

        ///// debug methods ///// 
        bool TestSolver(const int &modeIndex, const REAL &k, const Vector3d &testingPoint, const REAL &evaluateRadius = -1); 
        bool TestSolver(const int &modeIndex, const REAL &k, const Vector3d &testingPoint, REAL &residual, const REAL &evaluateRadius = -1); 
        void PrintFBemVelocityBC(const int &mode, const std::string &outFile); 
};

#endif
