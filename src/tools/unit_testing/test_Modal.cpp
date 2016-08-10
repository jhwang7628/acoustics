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
    impulseSeriesObject.AddImpulse(0.02, 1, Vector3d(1,2,3)); 
    impulseSeriesObject.AddImpulse(0.03, 4, Vector3d(3,2,1)); 
    REAL timestamp; 
    int vertex; 
    Vector3d impulse; 
    impulseSeriesObject.GetImpulse(1, timestamp, vertex, impulse); 
    REAL firstImpulseTime, lastImpulseTime; 
    impulseSeriesObject.GetImpulseRange(firstImpulseTime, lastImpulseTime); 

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
    object.GetImpulseRange(a, b); 
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
    const std::string solutionFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/ret-0_0.txt");
    reader.CheckFBemInputAgainstMesh(mesh, solutionFile);
    reader.ReadFBemOutputToInfo(solution, solutionFile); 
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
    TestBEMSolution();
    return 0; 
}
