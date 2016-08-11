#include <modal_model/KirchhoffIntegralSolver.h>
#include <io/FBemReader.h>
#ifdef USE_OPENMP
    #include <omp.h>
#endif

//##############################################################################
// Static member variable definitions.
//##############################################################################
const std::complex<REAL> KirchhoffIntegralSolver::j = std::complex<REAL>(0.0, 1.0); 

//##############################################################################
// This function reads in FBEM solution file. See FBemReader class for 
// file format.
//##############################################################################
void KirchhoffIntegralSolver::
AddFBemSolution(const std::string &fBemConfigFile, const std::string &fBemOutputFile, const REAL &omega)
{
    // need mesh to execute check
    assert(_mesh); 

    FBemReader fBemReader; 
    std::shared_ptr<BEMSolutionMode> solution(new BEMSolutionMode());
    fBemReader.CheckFBemInputAgainstMesh(_mesh, fBemConfigFile); 
    fBemReader.ReadFBemOutputToInfo(_mesh, fBemOutputFile, omega, solution); 

    _BEMSolutions.push_back(solution); 
}

//##############################################################################
// This function evalutes free-space green's function.
// @param k wavenumber = 2 pi f / c, f: frequency (Hz), c: wave speed.
//##############################################################################
std::complex<REAL> KirchhoffIntegralSolver::
Evaluate_G(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const REAL &rCached) const
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
Evaluate_dG_dn(const REAL &k, const Vector3d &listeningPoint, const Vector3d &surfacePoint, const Vector3d &normal, const REAL &rCached) const
{
    const Vector3d &x = listeningPoint; 
    const Vector3d &y = surfacePoint; 
    const Vector3d &n = normal;
    const REAL r = (rCached > 0 ? rCached : (listeningPoint - surfacePoint).length()); 
    return -std::exp(-j*k*r)/(4.0*M_PI*pow(r,3))*(1.0+j*k*r)*((x-y).dotProduct(n)); 
}

//##############################################################################
// This function solves the Kirchhoff integral and return the value as complex
// value at listening point. 
//
// @param modeIndex index in the solution array
// @param listeningPoint 
//##############################################################################
std::complex<REAL> KirchhoffIntegralSolver::
Solve(const int &modeIndex, const Vector3d &listeningPoint) const
{
    // get BEM solution and calculate wavenumber
    std::shared_ptr<BEMSolutionMode> bemSolution = _BEMSolutions.at(modeIndex); 
    const REAL k = bemSolution->omega / bemSolution->soundSpeed; 

    // initialize fields
#ifdef USE_OPENMP
    const int N_threads = omp_get_max_threads(); 
#else 
    const int N_threads = 1; 
#endif
    std::vector<std::complex<REAL> > transferAllThreads(N_threads, std::complex<REAL>(0.0, 0.0)); 

    // loop through all surface elements
    const int N_elements = _mesh->num_triangles(); 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int e_idx=0; e_idx<N_elements; ++e_idx)
    {
        const Tuple3ui &triangle = _mesh->triangle_ids(e_idx); 
        const Point3d &p0 = _mesh->vertex(triangle.x); 
        const Point3d &p1 = _mesh->vertex(triangle.y); 
        const Point3d &p2 = _mesh->vertex(triangle.z); 

        // compute geometry information about element
        Vector3d normal = ((p1 - p0).crossProduct(p2 - p0)) * 0.5;  // see TriangleMesh.hpp::generate_normals()
        const REAL area = normal.length(); 
        normal.normalize(); 
        const Vector3d surfacePoint = (p0 + p1 + p2)/3.0;

        // compute necessary fields
        const std::complex<REAL> G = Evaluate_G(k, listeningPoint, surfacePoint); 
        const std::complex<REAL> dG_dn = Evaluate_dG_dn(k, listeningPoint, surfacePoint, normal); 
        const std::complex<REAL> p = bemSolution->pressures.at(e_idx); 
        const std::complex<REAL> dp_dn = -j * bemSolution->omega * bemSolution->density * bemSolution->velocities.at(e_idx); 

        // evaluate element transfer value and add to accumulator
#ifdef USE_OPENMP
        const int threadId = omp_get_thread_num(); 
#else
        const int threadId = 0; 
#endif
        transferAllThreads.at(threadId) += (G * dp_dn - dG_dn * p) * area;
    }

    return std::accumulate(transferAllThreads.begin(), transferAllThreads.end(), std::complex<double>(0, 0));
}
