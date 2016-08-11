#include <modal_model/KirchhoffIntegralSolver.h>
#include <io/FBemReader.h>

//##############################################################################
// Static member variable definitions.
//##############################################################################
const std::complex<REAL> KirchhoffIntegralSolver::j = std::complex<REAL>(0.0, 1.0); 

//##############################################################################
// This function reads in FBEM solution file. See FBemReader class for 
// file format.
//##############################################################################
void KirchhoffIntegralSolver::
AddFBemSolution(const std::string &fBemConfigFile, const std::string &fBemOutputFile)
{
    // need mesh to execute check
    assert(_mesh); 

    FBemReader fBemReader; 
    std::shared_ptr<BEMSolutionMode> solution(new BEMSolutionMode());
    fBemReader.CheckFBemInputAgainstMesh(_mesh, fBemConfigFile); 
    fBemReader.ReadFBemOutputToInfo(_mesh, fBemOutputFile, solution); 

    _BEMSolutions.push_back(solution); 
}

//##############################################################################
// This function evalutes free-space green's function.
// @param k wavenumber = 2 pi f / c, f: frequency (Hz), c: wave speed.
//##############################################################################
std::complex<REAL> KirchhoffIntegralSolver::
Evaluate_G(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const REAL &rCached)
{
    const REAL r = (rCached > 0 ? rCached : (listeningPoint - surfacePoint).length()); 
    return std::exp(-j*k*r) / (4.0*M_PI*r); 
}

//##############################################################################
// This function evalutes normal gradient of free-space green's function.
// The formula for the normal derivatives is derived on paper and matches with,
// for example, the description in 
//
//  gbenthien.net/papers/KirchhoffIntEq.pdf
//##############################################################################
std::complex<REAL> KirchhoffIntegralSolver::
Evaluate_dG_dn(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const Vector3d &normal, const REAL &rCached)
{
    const Vector3d &x = listeningPoint; 
    const Vector3d &y = surfacePoint; 
    const Vector3d &n = normal;
    const REAL r = (rCached > 0 ? rCached : (listeningPoint - surfacePoint).length()); 
    return -std::exp(-j*k*r)/(4.0*M_PI*pow(r,3))*(1.0+j*k*r)*((x-y).dotProduct(n)); 
}
