#include <string>
#include <iostream>
#include <utils/STL_Wrapper.h> 
#include <deformable/ModeData.h> 
#include <io/TglMeshReader.hpp>
#include <io/RigidObjDispReader.h> 
#include <io/ImpulseSeriesReader.h> 
#include <modal_model/ModalAnalysisObject.h>
//#include <sndgen/RigidModal.h> 
#include <modal_model/ImpulseSeriesObject.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <parser/ImpulseResponseParser.h>
#include <modal_model/ModalMaterialList.h>
#include <modal_model/ModalMaterial.h>
#include <modal_model/ModalODESolver.h>
#include <modal_model/BEMSolutionMode.h> 
#include <modal_model/KirchhoffIntegralSolver.h>
#include <io/FBemReader.h>
#include <libconfig.h++> 

void TestRigidBodySim()
{
    const std::string testDisplacement("/home/jui-hsien/code/acoustics/work/plate_drop_test/displace.bin");
    const std::string testVelocity("/home/jui-hsien/code/acoustics/work/plate_drop_test/velocity.bin");
    const std::string testAcceleration("/home/jui-hsien/code/acoustics/work/plate_drop_test/acceleration.bin");
    RigidObjDispReader reader; 
    reader.ReadDisplacement(testDisplacement); 
    reader.ReadAllKinematics(testDisplacement, testVelocity, testAcceleration); 
}

void TestModal()
{
    std::cout << "test modal\n"; 
    const std::string filename("/home/jui-hsien/code/acoustics/work/plate_drop_test/example2.cfg"); 
    //ParseConfigModal parser(filename); 
    //parser.Parse(); 
    //libconfig::Config config; 
    //RigidModal rigidModal(filename); 
    //ModalAnalysis modalAnalysis(filename); 
    //modalAnalysis.BuildModalModelsFromFile(); 

    const std::string meshName("/home/jui-hsien/code/acoustics/work/plate/plate.obj"); 
    const std::string impulseFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/modalImpulses.txt"); 
    const std::string rigidsimConfigFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/default.cfg");

    std::shared_ptr<TriangleMesh<REAL> > mesh = std::make_shared<TriangleMesh<REAL> >(); 
    if (MeshObjReader::read(meshName.c_str(), *mesh, false, false, 1.0) == SUCC_RETURN)
        mesh->generate_normals();
    else 
        throw std::runtime_error("Object read failed.");
    ImpulseSeriesObject impulseSeriesObject(mesh); 
    ImpulseSeriesObject::ImpactRecord record1, record2; 
    record1.timestamp = 0.02; 
    record1.appliedVertex = 1; 
    record1.impactVector = Vector3d(1,2,3); 

    record2.timestamp = 0.03; 
    record2.appliedVertex = 4; 
    record2.impactVector = Vector3d(3,2,1); 

    impulseSeriesObject.AddImpulse(record1); 
    impulseSeriesObject.AddImpulse(record2); 
    REAL timestamp; 
    int vertex; 
    Vector3d impulse; 
    impulseSeriesObject.GetImpulse(1, timestamp, vertex, impulse); 
    REAL firstImpulseTime, lastImpulseTime; 
    impulseSeriesObject.GetRangeOfImpulses(firstImpulseTime, lastImpulseTime); 

    std::cout << impulseSeriesObject.Size() << " " << timestamp << " " << vertex << " " << impulse << " " << firstImpulseTime << " " << lastImpulseTime << std::endl;

    typedef std::shared_ptr<ImpulseSeriesObject> ImpulseSeriesObjectPtr; 
    ImpulseSeriesReader impulseSeriesReader(impulseFile, rigidsimConfigFile); 
    std::vector<ImpulseSeriesObjectPtr> objects(2, std::make_shared<ImpulseSeriesObject>()); 
    objects.at(0)->SetMesh(mesh); 
    objects.at(1)->SetMesh(mesh); 
    impulseSeriesReader.LoadImpulses(objects); 

    ImpulseSeriesObjectPtr anotherObject = std::make_shared<ImpulseSeriesObject>(mesh);
    impulseSeriesReader.LoadImpulses(0, anotherObject);
}

void TestRigidSoundObject()
{
    std::cout << "test RigidSoundObject\n";

    const std::string workingDirectory("/home/jui-hsien/code/acoustics/work/meshes/small_ball");
    const std::string prefix("small_ball"); 
    const int sdfResolution = 100; 
    const bool buildFromTet = false; 
    FDTD_RigidSoundObject object(workingDirectory, sdfResolution, prefix, buildFromTet); 
    std::cout << object.Initialized() << std::endl; 
    //object.SetID(5); 
    //object.SetMesh(object.GetMeshPtr()); 
    Vector3d impulse(1,2,3); 
    //object.AddImpulse(0.5, 10, impulse);

    REAL a, b; 
    object.GetRangeOfImpulses(a, b); 
    std::cout << a << " " << b << std::endl;
}

void TestIO()
{
    const std::string filename("/home/jui-hsien/code/acoustics/work/ball_0p1/small_ball_r_0p10.modes"); 
    ModeData modeData; 
    modeData.read(filename.c_str()); 
    std::cout << modeData << std::endl;
    for (int mm=0; mm<modeData.numModes(); ++mm) 
    {
        std::cout << "Mode " << mm << ": " << modeData._omegaSquared.at(mm) << std::endl;
        std::cout << " shape: "; 
        STL_Wrapper::PrintVectorContent(std::cout, modeData.mode(mm), 8); 
    }
}

void TestMaterialParser() 
{
    const std::string filename("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml"); 
    ImpulseResponseParser parser(filename); 
    ModalMaterialList materials; 
    parser.GetModalMaterials(materials); 
    auto material = materials.at(0); 
    std::cout << material->id << " " << material->alpha << " " << material->beta << " " << material->density << " " << material->inverseDensity << " " << material->xi(1.0) << " " << material->omega_di(1) << std::endl;
}

void TestModalODESolver() 
{
    const std::string filename("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml"); 
    ImpulseResponseParser parser(filename); 
    ModalMaterialList materials; 
    parser.GetModalMaterials(materials); 
    auto materialPtr = materials.at(0); 
    const REAL dt = 0.01; 
    const REAL omegaSquared = 0.99;
    //const REAL Q = 0.0;
    ModalODESolver solver;
    solver.Initialize(materialPtr, omegaSquared, dt); 
    //solver.StepSystem(Q);
}

void TestBEMSolution()
{
    FBemReader reader; 

    // read mesh
    const std::string meshFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/proj.tet.obj"); 
    std::shared_ptr<TriangleMesh<REAL> > mesh = std::make_shared<TriangleMesh<REAL> >(); 
    if (MeshObjReader::read(meshFile.c_str(), *mesh, false, false, 1.0) == SUCC_RETURN)
        mesh->generate_normals();
    else 
        throw std::runtime_error("Object read failed.");

    // read FBem solution
    std::shared_ptr<BEMSolutionMode> solution = std::make_shared<BEMSolutionMode>(); 
    const std::string fBemInputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/input-0_0.txt");
    const std::string fBemOutputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/ret-0_0.txt");
    const REAL omega = 2.0 * M_PI * 1020.01;
    //reader.CheckFBemInputAgainstMesh(mesh, fBemInputFile);
    //reader.ReadFBemOutputToInfo(mesh, fBemOutputFile, omega, solution); 

    KirchhoffIntegralSolver solver(mesh);
    solver.AddFBemSolution(fBemInputFile, fBemOutputFile, omega);
    solver.Solve(0, Vector3d(0, 0, 0));

    solver.PrintFBemVelocityBC(0, "allVelocityFBem.txt");

    const int N_tests = 10; 
    const REAL pushRadius = 50.0;
    bool testFailed = false; 
    for (int t_idx=0; t_idx<N_tests; ++t_idx)
    {
        // pick a random mesh point and push it away from boundary to be
        // testing point
        const int vertexId = rand() % mesh->num_vertices(); 
        const Vector3d &randomVertex = mesh->vertex(vertexId); 
        Vector3d normal = mesh->normal(vertexId); 
        normal.normalize(); 
        const Vector3d testingPoint = randomVertex + normal * pushRadius; 

        // test solver
        const bool thisFailed = !solver.TestSolver(0, omega/343.0, testingPoint); 
        testFailed = (testFailed || thisFailed); 
    }

    std::cout << "overall test " << (testFailed ? "failed" : "passed") << std::endl; 
}

int main() 
{

    std::cout << "Unit Test: Modal Sound\n"; 
    //TestIO(); 
    //TestRigidBodySim(); 
    //TestModal(); 
    //TestRigidSoundObject(); 
    //TestMaterialParser();
    //TestModalODESolver(); 
    //TestBEMSolution();
    return 0; 
}
