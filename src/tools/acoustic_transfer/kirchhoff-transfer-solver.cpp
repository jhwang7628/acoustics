#include <iostream> 
#include <fstream> 
#include <memory>
#include <string> 
#include <config.h>
#include <io/TglMeshReader.hpp>
#include <modal_model/KirchhoffIntegralSolver.h> 
#include <geometry/TriangleMesh.hpp> 

//##############################################################################
//##############################################################################
REAL Compute_q(const REAL &omega, const REAL &t) 
{
    return cos(omega * t); 
}

//##############################################################################
//##############################################################################
void Test_PerfectHarmonic()
{
    std::cout << "test\n";

    // read mesh
    const std::string meshFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/proj.tet.obj"); 
    std::shared_ptr<TriangleMesh<REAL> > mesh = std::make_shared<TriangleMesh<REAL> >(); 
    if (MeshObjReader::read(meshFile.c_str(), *mesh, false, false, 1.0) == SUCC_RETURN)
        mesh->generate_normals();

    const std::string fBemInputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/input-0_0.txt");
    const std::string fBemOutputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/ret-0_0.txt");
    const REAL omega = 2.0 * M_PI * 1020.01;
    const Vector3d listeningPoint(0.0, 0.0, 0.0825); 

    KirchhoffIntegralSolver solver(mesh); 
    solver.AddFBemSolution(fBemInputFile, fBemOutputFile, omega);
    std::complex<REAL> transferValue = solver.Solve(0, listeningPoint);

    const REAL tStart = 0.0; 
    const REAL tStep = 1.0/176400.0; 
    const REAL transferMagnitude = std::abs(transferValue); 
    const int steps = 35280; 

    std::cout << "Transfer magnitude at point " << listeningPoint << " = " << transferValue << "; its magnitude = " << transferMagnitude << std::endl;

    std::ofstream of("Test_PerfectHarmonic.txt"); 
    of << std::setprecision(16);
    for (int s_idx=0; s_idx<steps; ++s_idx)
    {
        const REAL tNow = tStart + tStep * (REAL)s_idx; 
        const REAL value = transferMagnitude * Compute_q(omega, tNow); 
        of << tNow << " " << value << std::endl;
    }
    of.close(); 
}

//##############################################################################
//##############################################################################
int main()
{
    Test_PerfectHarmonic(); 
    return 0;
}
