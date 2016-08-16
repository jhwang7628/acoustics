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
// value at listening point. If the listening point is too close to the
// surface, as compared to variable _distanceLowBound, then return 0.
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
    bool abortCompute = false; 

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
        // The normal is flipped because this is the normal for the Helmholtz
        // exterior problem. Normal of this domain should point "into" the mesh.
        Vector3d normal = -((p1 - p0).crossProduct(p2 - p0)) * 0.5;  // see TriangleMesh.hpp::generate_normals() and fbem_input_gen.cpp:218
        const REAL area = normal.length(); 
        normal.normalize(); 
        const Vector3d surfacePoint = (p0 + p1 + p2)/3.0;

        const REAL r = (listeningPoint - surfacePoint).length(); 
        if (r < _distanceLowBound)
            abortCompute = true; // this listening point is too close to the surface

        // compute necessary fields
        const std::complex<REAL> G = Evaluate_G(k, listeningPoint, surfacePoint, r); 
        const std::complex<REAL> dG_dn = Evaluate_dG_dn(k, listeningPoint, surfacePoint, normal, r); 
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

    if (abortCompute)
        return std::complex<REAL>(0.0, 0.0);

    return std::accumulate(transferAllThreads.begin(), transferAllThreads.end(), std::complex<double>(0, 0));
}

//##############################################################################
//##############################################################################
bool KirchhoffIntegralSolver::
TestSolver(const int &modeIndex, const REAL &k, const Vector3d &testingPoint, const REAL &evaluateRadiusArg)
{
    REAL residual; 
    return TestSolver(modeIndex, k, testingPoint, residual, evaluateRadiusArg); 
}

//##############################################################################
// This function tests the solver by solving a small neighbour hood around
// testing point and check if it satisfies the Helmholtz equation using
// finite-difference. 
//
//  Helmholtz eq: 
//   laplacian(p(x)) + k^2 p(x) = 0
//
//  @param modeIndex
//  @param k wavenumber = omega / wavespeed
//  @param testingPoint
//  @param evaluateRadiusArg Used if >0, otherwise default value is used. This is
//                           the finite-difference stencil interval
//##############################################################################
bool KirchhoffIntegralSolver::
TestSolver(const int &modeIndex, const REAL &k, const Vector3d &testingPoint, REAL &residual, const REAL &evaluateRadiusArg)
{
    const REAL evaluateRadius = (evaluateRadiusArg > 0 ? evaluateRadiusArg : _defaultTestEvaluateRadius); 

    // evaluate pressure. notation:
    //     h
    //     |
    // l---0---h
    //     | 
    //     l
    const Vector3d &position000 = testingPoint; 
    const Vector3d  positionl00 = testingPoint - Vector3d(evaluateRadius, 0.0, 0.0); 
    const Vector3d  positionh00 = testingPoint + Vector3d(evaluateRadius, 0.0, 0.0); 
    const Vector3d  position0l0 = testingPoint - Vector3d(0.0, evaluateRadius, 0.0); 
    const Vector3d  position0h0 = testingPoint + Vector3d(0.0, evaluateRadius, 0.0); 
    const Vector3d  position00l = testingPoint - Vector3d(0.0, 0.0, evaluateRadius); 
    const Vector3d  position00h = testingPoint + Vector3d(0.0, 0.0, evaluateRadius); 
    const std::complex<REAL> p000 = Solve(modeIndex, position000);
    const std::complex<REAL> pl00 = Solve(modeIndex, positionl00);
    const std::complex<REAL> ph00 = Solve(modeIndex, positionh00);
    const std::complex<REAL> p0l0 = Solve(modeIndex, position0l0);
    const std::complex<REAL> p0h0 = Solve(modeIndex, position0h0);
    const std::complex<REAL> p00l = Solve(modeIndex, position00l);
    const std::complex<REAL> p00h = Solve(modeIndex, position00h);

    // evaluate laplacian(p) using 2nd central difference
    const REAL inv_r2 = 1.0 / pow(evaluateRadius, 2); 
    const std::complex<REAL> laplacian_p = ((ph00 + pl00) + (p0h0 + p0l0) + (p00h + p00l) - 6.0*p000) * inv_r2; 
    residual = std::abs(laplacian_p + pow(k, 2) * p000); 
    const bool testPassed = (residual < _defaultTestErrorTolerance); 

    std::cout << "\n\nTEST HELMHOLTZ RESIDUAL\n";
    std::cout << "TestSolver()::\n"
              << " Evaluate point = " << testingPoint << "\n"
              << " Evaluate radius = " << evaluateRadius << "\n"
              << " Tested residual = " << residual << "\n"; 
    if (!testPassed)
        std::cout << " **Test failed\n";
    else
        std::cout << " **Test passed\n";
    std::cout << std::flush; 

    return testPassed; 
}

//##############################################################################
//##############################################################################
void KirchhoffIntegralSolver::
PrintFBemVelocityBC(const int &mode, const std::string &outFile)
{
    std::shared_ptr<BEMSolutionMode> solution = _BEMSolutions.at(mode); 
    const int N_elements = _mesh->num_triangles(); 
    std::ofstream of(outFile.c_str()); 
    for (int e_idx=0; e_idx<N_elements; ++e_idx) 
    {
        const auto &velocity = solution->velocities.at(e_idx);
        of << velocity.real() << " " << velocity.imag() << "\n";
    }
    of.close(); 
}
