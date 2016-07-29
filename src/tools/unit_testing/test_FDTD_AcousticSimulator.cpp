#include <tools/unit_testing/testing.h>
#include <wavesolver/FDTD_MovableObject.h>
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/FDTD_Objects.h> 
#include <wavesolver/FDTD_AcousticSimulator.h> 
#include <wavesolver/FDTD_RigidObject_Animator.h>
#include <parser/ImpulseResponseParser.h> 
#include <linearalgebra/Vector3.hpp>
#include <geometry/BoundingBox.h> 
#include <boost/timer/timer.hpp>

//##############################################################################
// Submodules 
//##############################################################################
void TestParseMeshList()
{
    const std::string xmlName("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml"); 

    std::shared_ptr<FDTD_Objects> objectsInTheScene = std::make_shared<FDTD_Objects>(); 
    ImpulseResponseParser parser(xmlName); 

    parser.GetObjects(objectsInTheScene); 
    std::cout << *objectsInTheScene << std::endl;
    //objectsInTheScene->TestObjectDistanceField(0);
}

//##############################################################################
void TestBoundingBox()
{
    const std::string meshFileName("/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj");
    const std::string sdfFilePrefix("/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj.1m.dist");
    const int sdfResolution = 100; 
    FDTD_RigidObject object(meshFileName, sdfResolution, sdfFilePrefix, false); 
    //object.Initialize();

    SimpleTimer timer[2]; 
    timer[0].Start();
    const int N_test = 1000000;
    int p_return=0, n_return=0; 
    srand(time(NULL));
    for (int ii=0; ii<N_test; ii++) 
    {
        //const Vector3d vec(0.2,0.3,0.5); 
        const double distance = object.DistanceToMesh((double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX);
        if (distance > 0) 
            p_return ++; 
        else 
            n_return ++; 
    }
    timer[0].Pause();

    object.PrintBoundingBox(); 
    object.UpdateBoundingBox();
    object.PrintBoundingBox(); 

    COUT_SDUMP(p_return); 
    COUT_SDUMP(n_return); 

    p_return =0; 
    n_return =0;

    timer[1].Start();
    for (int ii=0; ii<N_test; ii++) 
    {
        //const Vector3d vec(0.2,0.3,0.5); 
        const double distance = object.DistanceToMesh((double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX);
        if (distance > 0) 
            p_return ++; 
        else 
            n_return ++; 
    }
    timer[1].Pause();
    std::cout << "time elapsed : " << timer[0].Duration() << std::endl;
    std::cout << "time elapsed : " << timer[1].Duration() << std::endl;

    COUT_SDUMP(p_return); 
    COUT_SDUMP(n_return); 
}

//##############################################################################
void TestAcousticSimulatorRun(const std::string &xmlName)
{
    FDTD_AcousticSimulator simulator(xmlName);
    simulator.InitializeSolver(); 
    simulator.Run();
    //simulator.TestAllComponents();
}

//##############################################################################
void TestScalarFieldSubindices()
{
    const Vector3d boxMinBound(0,0,0); 
    const Vector3d boxMaxBound(1,1,1); 
    const BoundingBox bbox(boxMinBound, boxMaxBound); 
    ScalarField field(bbox, 0.002); 
    field.TestSubindices();
}

//##############################################################################
void TestFDTD_RigidObject_Animator()
{
    const std::string displacementFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/displace.bin"); 
    FDTD_RigidObject_Animator animator; 
    animator.ReadDisplacement(displacementFile); 
    std::ofstream of("test.txt"); 
    Vector3d displacement; 
    Quaternion<REAL> quaternion;
    const REAL dt = 0.0001; 
    REAL t = 0.0; 
    while(t<1)
    {
        animator.GetObjectDisplacement(0, t, displacement, quaternion); 
        t += dt; 

        of << displacement.x << " " << displacement.y << " " << displacement.z << std::endl;
    }
    of.close();
}

//##############################################################################
int main(int argc, char ** argv)
{
    //TestBoundingBox();
    //TestParseMeshList(); 
    //TestScalarFieldSubindices();
    
    //TestFDTD_RigidObject_Animator();
    std::string xmlName("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml");
    if (argc>1) 
        xmlName = std::string(argv[1]);
    TestAcousticSimulatorRun(xmlName); 

    return 0;
}
