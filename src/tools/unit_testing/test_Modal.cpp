#include <string>
#include <iostream>
#include <utils/STL_Wrapper.h> 
#include <deformable/ModeData.h> 
#include <io/TglMeshReader.hpp>
#include <io/RigidObjDispReader.h> 
#include <io/ImpulseSeriesReader.h> 
#include <modal_model/ModalAnalysis.h>
#include <sndgen/RigidModal.h> 
#include <modal_model/ImpulseSeriesObject.h>
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
    ModalAnalysis modalAnalysis(filename); 
    modalAnalysis.BuildModalModelsFromFile(); 

    const std::string meshName("tmp.obj"); 
    const std::string impulseFile("modalImpulses.txt"); 

    std::shared_ptr<TriangleMesh<REAL> > mesh = std::make_shared<TriangleMesh<REAL> >(); 
    if (MeshObjReader::read(meshName.c_str(), *mesh, false, false, 1.0) == SUCC_RETURN)
        mesh->generate_normals();
    ImpulseSeriesObject impulseSeriesObject(0, mesh); 
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
    ImpulseSeriesReader impulseSeriesReader(impulseFile); 
    std::vector<ImpulseSeriesObjectPtr> objects; 
    impulseSeriesReader.LoadImpulses(objects); 
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

int main() 
{
    std::cout << "Unit Test: Modal Sound\n"; 
    //TestIO(); 
    //TestRigidBodySim(); 
    TestModal(); 
    return 0; 
}
