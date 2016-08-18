#include <iostream> 
#include <fstream> 
#include <memory>
#include <string> 
#include <config.h>
#include <io/TglMeshReader.hpp>
#include <modal_model/KirchhoffIntegralSolver.h> 
#include <geometry/TriangleMesh.hpp> 
#include <boost/program_options.hpp>

static Vector3d listeningPoint; 

//##############################################################################
// This function parses input arguments. 
//##############################################################################
static void parse_cmd(int argc, char **argv)
{
    std::cout << "\nPARSE COMMAND LINE INPUTS: \n";
    namespace po = boost::program_options; 
    po::options_description description("Allowed options");
    description.add_options()
        ("help", "Display help message")
        ("listening_point", po::value<std::vector<REAL> >()->multitoken(), "Listening Point for transfer function"); 

    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, description, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm); 
    po::notify(vm); 

    std::vector<REAL> tmp_listeningPoint; 
    if (!vm["listening_point"].empty() && (tmp_listeningPoint = vm["listening_point"].as<std::vector<REAL> >()).size() == 3)
    {
        listeningPoint.x = tmp_listeningPoint.at(0);
        listeningPoint.y = tmp_listeningPoint.at(1);
        listeningPoint.z = tmp_listeningPoint.at(2);
        std::cout << " listening point = " << listeningPoint << "\n\n";
    }
    else 
    {
        std::cerr << "**ERROR** Specify three floating point numbers after listening flag"; 
        std::cerr << description << std::endl;
        exit(1);
    }
}

//##############################################################################
//##############################################################################
REAL Compute_q(const REAL &omega, const REAL &t) 
{
    return cos(omega * t); 
}

//##############################################################################
//##############################################################################
void SingleFrequencySolve()
{
    std::cout << "\nSINGLE FREQUENCY SOLVE: \n";

    // read mesh
    const std::string meshFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/proj.tet.obj"); 
    std::shared_ptr<TriangleMesh<REAL> > mesh = std::make_shared<TriangleMesh<REAL> >(); 
    if (MeshObjReader::read(meshFile.c_str(), *mesh, false, false, 1.0) == SUCC_RETURN)
        mesh->generate_normals();

    const std::string fBemInputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/input-0_0.txt");
    const std::string fBemOutputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/ret-0_0.txt");
    const REAL omega = 2.0 * M_PI * 1020.01;
    //const Vector3d listeningPoint(0.0, 0.0, 0.0825); 

    KirchhoffIntegralSolver solver(mesh); 
    solver.AddFBemSolution(fBemInputFile, fBemOutputFile, omega);
    std::complex<REAL> transferValue = solver.Solve(0, listeningPoint);

    // additional scaling due to fbem input scaling. In particular, the 
    // fbem input in dataset is missing a 2pi scaling for velocity BC. 
    // since the helmholtz equation is linear, I can scale the output solution 
    // by this factor. Notice another bug in fbem_input_gen.cpp, when the 
    // velocity BC is constructed, there should be an extra minus sign since 
    // its the normal for the Helmholtz solve, which points inside the mesh. 
    // Added a minus sign here to compensate. 
    const REAL extraScaling = -2.0 * M_PI; 
    transferValue *= extraScaling; 
    std::cout << "Transfer value scaled by: " << extraScaling << std::endl;

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
int main(int argc, char **argv)
{
    parse_cmd(argc, argv);
    SingleFrequencySolve(); 
    return 0;
}
