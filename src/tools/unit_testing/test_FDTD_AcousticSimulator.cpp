#include <tools/unit_testing/testing.h>
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/FDTD_Objects.h> 
#include <parser/ImpulseResponseParser.h> 
#include <linearalgebra/Vector3.hpp>

//##############################################################################
// Submodules 
//##############################################################################
void TestParseMeshList()
{
    const std::string xmlName("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml"); 

    FDTD_Objects objectsInTheScene; 
    ImpulseResponseParser parser(xmlName); 

    parser.GetObjects("test", objectsInTheScene); 
    std::cout << objectsInTheScene << std::endl;
}

//##############################################################################
void TestBoundingBox()
{
    const std::string meshFileName("/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj");
    const std::string sdfFilePrefix("/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj.1m.dist");
    const int sdfResolution = 100; 
    FDTD_RigidObject object(meshFileName, sdfResolution, sdfFilePrefix); 
    object.Initialize();

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
int main()
{
    //TestBoundingBox();
    TestParseMeshList(); 

    return 0;
}
