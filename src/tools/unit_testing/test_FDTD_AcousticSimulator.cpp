#include <tools/unit_testing/testing.h>
#include <wavesolver/FDTD_MovableObject.h>
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/FDTD_Objects.h> 
#include <wavesolver/FDTD_AcousticSimulator.h> 
#include <wavesolver/FDTD_RigidObject_Animator.h>
#include <wavesolver/WaterVibrationalSource.h>
#include <geometry/TriangleMeshKDTree.hpp>
#include <parser/ImpulseResponseParser.h> 
#include <linearalgebra/Vector3.hpp>
#include <geometry/BoundingBox.h> 
#include <boost/timer/timer.hpp>
#include <io/TglMeshReader.hpp>
#include <sndgen/WavReader.hpp>

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
    //srand(time(NULL));
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
    //simulator.TestAnimateObjects(1);
    //simulator.Run();
    simulator.RunForSteps(1);
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
void Test_TriangleMeshKDTree()
{
    const std::string meshFile("/home/jui-hsien/code/acoustics/work/plate/plate.tet.obj"); 
    std::shared_ptr<TriangleMesh<REAL> > mesh = std::make_shared<TriangleMeshKDTree<REAL> >(); 
    MeshObjReader::read(meshFile.c_str(), *mesh, false, false, 1.0); 
    mesh->generate_normals(); 
    std::cout << "N_vertices = " << mesh->num_vertices() << std::endl;

    std::dynamic_pointer_cast<TriangleMeshKDTree<REAL> >(mesh)->BuildKDTree();
    Vector3d queryPoint, projectedPoint, closestPoint; 
    int closestTriangle;

    const int N_tests = 50; 
    std::ofstream of("Test_TriangleMeshKDTree.txt"); 
    const double r = mesh->boundingSphereRadius(Point3d(0, 0, 0)); 
    //srand(time(NULL)); 
    for (int t_idx=0; t_idx<N_tests; ++t_idx)
    {
        std::cout << "\r" << t_idx << "/" << N_tests << std::endl;
        const double x = (double)rand()/(double)RAND_MAX * 2.0*r - r; 
        const double y = (double)rand()/(double)RAND_MAX * 2.0*r - r; 
        const double z = (double)rand()/(double)RAND_MAX * 2.0*r - r; 
        queryPoint.set(x, y, z);
        //queryPoint.set(-0.101063680118473, 0.1055460887829576, 0.006106259099252123);


        mesh->ComputeClosestPointOnMesh(queryPoint, closestPoint, closestTriangle, projectedPoint); 
        //const Tuple3ui &triangle = mesh->triangle_ids(closestTriangle);
        const Vector3d normal = closestPoint - queryPoint; 
        of << std::setprecision(16); 
        of << queryPoint.x << " " << queryPoint.y << " " << queryPoint.z << " "
           << normal.x << " " << normal.y << " " << normal.z << "\n";
    }
    std::cout << std::endl;
    of.close(); 
}

//##############################################################################
void TestWavRead()
{
    std::vector<double> data;
    WavReader<double> reader; 
    const std::string wavFile("/home/jui-hsien/code/acoustics/work/droplet_recordings/36_nr.wav");
    reader.Open(wavFile); 
    reader.Read(data); 
    reader.Close(); 
    STL_Wrapper::PrintVectorContent(std::cout, data, 10); 
    std::cout << std::endl;

    std::vector<float> f_data;
    WavReader<float> f_reader; 
    //const std::string wavFile("/home/jui-hsien/code/acoustics/work/droplet_recordings/36_nr.wav");
    f_reader.Open(wavFile); 
    f_reader.Read(f_data); 
    f_reader.Close(); 
    STL_Wrapper::PrintVectorContent(std::cout, f_data, 10); 
    std::cout << std::endl;

    std::vector<short> s_data;
    WavReader<short> s_reader; 
    //const std::string wavFile("/home/jui-hsien/code/acoustics/work/droplet_recordings/36_nr.wav");
    s_reader.Open(wavFile); 
    s_reader.Read(s_data); 
    s_reader.Close(); 
    STL_Wrapper::PrintVectorContent(std::cout, s_data, 10); 
    std::cout << std::endl;

    std::vector<int> i_data;
    WavReader<int> i_reader; 
    //const std::string wavFile("/home/jui-hsien/code/acoustics/work/droplet_recordings/36_nr.wav");
    i_reader.Open(wavFile); 
    i_reader.Read(i_data); 
    i_reader.Close(); 
    STL_Wrapper::PrintVectorContent(std::cout, i_data, 10); 
    std::cout << std::endl;
}

//##############################################################################
int main(int argc, char ** argv)
{
    //TestBoundingBox();
    //TestParseMeshList(); 
    //TestScalarFieldSubindices();
    //TestFDTD_RigidObject_Animator();
    //Test_TriangleMeshKDTree();
    //TestWaterVibrationalSource();
    //TestWavRead(); 
        
    std::string xmlName("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml");
    if (argc>1) 
        xmlName = std::string(argv[1]);
    TestAcousticSimulatorRun(xmlName); 

    return 0;
}
