#include <iostream> 
#include <fstream>
#include "Grid.h" 
#include "TetrahedronMesh.h" 
#include "OpenfoamMesh.h" 
#include "STL_Utils/STL_Wrapper.h"
#include "signalprocessing/FilterDesign.h"

#include <boost/timer/timer.hpp> 


void TestGrid(); 
void TestOpenfoamMesh();
void TestTetrahedralMesh(); 

int main() 
{

    //TestOpenfoamMesh(); 
    TestGrid(); 
    //TestTetrahedralMesh(); 

    return 0; 
} 

void TestOpenfoamMesh() 
{
    const std::string root("/hdd1/research_data/turbsound_data/OpenFOAM-run/test_mesh"); 
    const std::string zone("FLUID_SMALL"); 

    OpenfoamMesh mesh(root,zone); 
    mesh.ReinitializeMesh(); 
    mesh.ReconstructTetrahedronMesh(); 

    Eigen::MatrixXd centroids; 
    mesh.GetCentroids(centroids); 

    std::ofstream of("data/centroids.csv"); 
    for (int ii=0; ii<centroids.rows(); ii++) 
        of << centroids(ii,0) << "," << centroids(ii,1) << "," << centroids(ii,2) << std::endl;
    of.close();



    //OpenfoamCase foamCase(root,zone); 
    //foamCase.ReinitializeCase(); 
    //const int N_timesteps = foamCase.PrepareDataRead("U", 3, 0.00001, -1); 
    //if (N_timesteps>0)
    //{
    //    DataTimestep dataTimestep;
    //    Eigen::MatrixXi meshGridProjectionTable; 
    //    Eigen::MatrixXd projectedData; 
    //    Eigen::Vector3d minBound(-0.2178337591671500,-0.2017749603171500,-0.2192579190471500); 
    //    Eigen::Vector3d maxBound(0.2204221281671500,0.2364809270171500,0.2189979682871500); 
    //    Eigen::Vector3i cellCount(250,250,250); 
    //    UniformGrid<double> grid(minBound, maxBound, cellCount);

    //    foamCase.ReadTimestep(18,dataTimestep); 
    //    foamCase.GetMeshGridProjectionTable(grid, "data/meshGridToProjectionTable.indicies", meshGridProjectionTable, true); 
    //    foamCase.ProjectDataTimestep(dataTimestep, meshGridProjectionTable, grid, projectedData); 
    //    foamCase.WriteDataTimestepVTK("data/firstTS_x.vtk", "Ux", projectedData.col(0), grid);
    //    foamCase.WriteDataTimestepVTK("data/firstTS_y.vtk", "Uy", projectedData.col(1), grid);
    //    foamCase.WriteDataTimestepVTK("data/firstTS_z.vtk", "Uz", projectedData.col(2), grid);
    //    foamCase.WriteDataTimestepVTK("data/firstTS_mag.vtk", "Umag", (projectedData.col(0).array().square()+
    //                projectedData.col(1).array().square()+
    //                projectedData.col(2).array().square()).cwiseSqrt(), grid);
    //}






    //OpenfoamMesh mesh(root,zone); 

    //mesh.ReinitializeMesh(); 
    //mesh.ReconstructTetrahedronMesh();


    //mesh.PrepareDataRead("Ug", 0.00015, 0.00073); 


    ///// test listdirectory /////
    //std::vector<std::string> testfilenames = IO::listDirectoryMatch(IO::AssembleFilePath(root,"processor0").c_str(), "[[:digit:]]+\\.?[[:digit:]]*");
    //STL_Wrapper::PrintVectorContent(std::cout, testfilenames, true);

    ///// unit testing /////
    //Eigen::MatrixXd points; 
    //mesh.ReadFoamFile("points", 3, points); 
    //std::cout << "\npoints:\n"; 
    //std::cout << points.topRows(3) << std::endl;
    //std::cout << points.bottomRows(3) << std::endl;


    //Eigen::MatrixXi cellZones; 
    //mesh.ReadFoamFile_Int("cellZones", 1, cellZones, "FLUID_SMALL"); 
    //std::cout << "\ncellZones:\n"; 
    //std::cout << cellZones.topRows(3) << std::endl;
    //std::cout << cellZones.bottomRows(3) << std::endl;
    //  
    //Eigen::MatrixXi owner; 
    //mesh.ReadFoamFile_Int("owner", 1, owner); 
    //std::cout << "\nowner:\n"; 
    //std::cout << owner.topRows(3) << std::endl;
    //std::cout << owner.bottomRows(3) << std::endl;

    //CompressedList faces; 
    //mesh.ReadFoamFileFace("faces", faces); 
    //std::cout << "\nfaces:\n"; 
    //std::cout << faces[0].transpose() << std::endl;
    //std::cout << faces[1].transpose() << std::endl;
    //std::cout << faces[2].transpose() << std::endl;
    //std::cout << faces[faces.Size()-3].transpose()<< std::endl;
    //std::cout << faces[faces.Size()-2].transpose()<< std::endl;
    //std::cout << faces[faces.Size()-1].transpose()<< std::endl;

}

void TestTetrahedralMesh() 
{
    srand((unsigned int) time(0));

    /*
       {

       Eigen::Vector3d p0(0,0,0); 
       Eigen::Vector3d p1(1,0,0); 
       Eigen::Vector3d p2(0,1,0); 
       Eigen::Vector3d p3(0,0,1); 
       Tetrahedron tet(p0, p1, p2, p3);

       Eigen::Vector3d cc = tet.Centroid(); 
       Eigen::Vector4d w; 

       std::cout << cc << std::endl;
       std::cout << "Cell Center is in tet : " << tet.Inside(cc,w) << std::endl;

       const int Nquery = 1E5; 
       for (int ii=0; ii<Nquery; ii++) 
       {
       if ( ii%1000==0 ) std::cout << ii << " " << std::flush; 
       Eigen::Vector4d w; 
       Eigen::Vector3d randomPoint = Eigen::Vector3d::Random(); 
       Eigen::Vector3d n(1,1,1); 

       const double t = 1.0 - randomPoint.dot(n); 
       bool inTetrahedronAnalytical; 

       if (randomPoint(0)<0 || randomPoint(1)<0 || randomPoint(2)<0 || t<0) 
       inTetrahedronAnalytical = false; 
       else 
       inTetrahedronAnalytical = true; 

       const bool inTetrahedron = tet.Inside(randomPoint, w); 

       if (inTetrahedron != inTetrahedronAnalytical)
       {
       std::cout << "====== " << ii << " =====\n"; 
       std::cout << "Random point : " << randomPoint.transpose() << std::endl; 
       std::cout << " is in the tet numerically  : " << inTetrahedron << std::endl; 
       std::cout << " is in the tet analytically : " << inTetrahedronAnalytical << std::endl;
       if (inTetrahedron) 
       {
       std::cout << "  weights : " << w.transpose() << std::endl; 
       }
       }
       }
       std::cout << std::endl;

       }
       */

    TetrahedronMesh mesh; 
    Eigen::Vector3d listen = Eigen::Vector3d::Random(); 

    double smallestDistance = std::numeric_limits<double>::max(); 
    int smallestIndex = std::numeric_limits<int>::max();

    {
        boost::timer::auto_cpu_timer tt("constructing 250cubed tets : %w\n"); 
        const int Ntets = 250*250*250; 
        for (int ii=0; ii<Ntets; ii++) 
        {
            Eigen::Vector3d p0 = Eigen::Vector3d::Random(); 
            Eigen::Vector3d p1 = Eigen::Vector3d::Random(); 
            Eigen::Vector3d p2 = Eigen::Vector3d::Random(); 
            Eigen::Vector3d p3 = Eigen::Vector3d::Random(); 

            Tetrahedron tet(p0, p1, p2, p3); 

            mesh.AddTetrahedron(tet); 
            double distance = (tet.Centroid()-listen).norm(); 

            if (distance < smallestDistance)
            {
                smallestIndex = ii; 
                smallestDistance = distance; 
            } 
            //std::cout << ii << " : " << distance << std::endl; 

        }
        std::cout << "number of tets in the mesh = " << mesh.N_tetrahedrons() << std::endl << std::endl;

        std::cout << "smallest should be = " << smallestIndex << " at " << smallestDistance << std::endl; 
    }

    {
        boost::timer::auto_cpu_timer tt("Building trees : %w\n"); 
        mesh.BuildTree(); 
    } 

    {
        boost::timer::auto_cpu_timer tt("Query the trees for 1E6 times: %w\n"); 

        for (int ii=0; ii<1000000; ii++) 
        {
            int index = mesh.FindNearestTetrahedron(listen); 
            Tetrahedron& nearestTet = mesh.Get(index);
            if (ii%1000==0)
                std::cout << index << " : " << nearestTet << std::endl;; 
        }

    } 




}


void TestGrid_TestSmoothing()
{

    ///// test on smoothing /////
    std::cout << "testing for smoothing\n"; 
    const int N=100; 
    Eigen::Vector3d minBound(-1,-1,-1); 
    Eigen::Vector3d maxBound( 1, 1, 1); 
    Eigen::Vector3i division( N, N, N); 

    Eigen::Vector3d dx = (maxBound - minBound).cwiseQuotient( division.cast<double>() ); 

    UniformGrid<double> grid( minBound, maxBound, division ); 

    std::shared_ptr<Eigen::MatrixXd> data_random(new Eigen::MatrixXd()); 

    srand (time(NULL));
    data_random->resize( N*N*N, 1 ); 
    *data_random = Eigen::MatrixXd::Random(N*N*N,1); 

    std::ofstream of("data/data_random.csv");
    for (int ii=0; ii<N; ii++) 
    {
        for (int jj=0; jj<N; jj++) 
        {
            for (int kk=0; kk<N; kk++) 
            {
                double x = minBound[0] + ii*(maxBound[0]-minBound[0])/division[0] + dx[0]/2.; 
                double y = minBound[1] + jj*(maxBound[1]-minBound[1])/division[1] + dx[1]/2.; 
                double z = minBound[2] + kk*(maxBound[2]-minBound[2])/division[2] + dx[2]/2.; 
                of << x << "," << y << "," << z << "," << (*data_random)(ii+jj*N+kk*N*N,0) << std::endl;
            }
        }
    }
    of.close();

    Eigen::VectorXd filter(3); 
    filter << 1,2,1; 
    filter /= 4.0; 

    grid.InsertCellCenteredData( "data_random", data_random ); 
    Eigen::MatrixXd data_random_smoothed = grid.CellCenteredSmoothing( "data_random", filter ); 

    Eigen::VectorXd filter2(5); 
    filter2 << 1,2,5,2,1; 
    filter2 /= filter2.sum(); 
    //grid.InsertCellCenteredData( "data_random_", data_random_smoothed ); 
    Eigen::MatrixXd data_random_smoothed2 = grid.CellCenteredSmoothing( "data_random", filter2 ); 

    Eigen::VectorXd filter3 = SIGNAL_PROCESSING::DiscreteGaussian1D(9,2);
    //grid.InsertCellCenteredData( "data_random_", data_random_smoothed ); 
    Eigen::MatrixXd data_random_smoothed3 = grid.CellCenteredSmoothing( "data_random", filter3 ); 

    std::ofstream ofsmoothed("data/cell_centered_smoothed.csv");
    std::ofstream ofsmoothed2("data/cell_centered_smoothed2.csv");
    std::ofstream ofsmoothed3("data/cell_centered_smoothed3.csv");
    for (int kk=0; kk<N; kk++) 
    {
        for (int jj=0; jj<N; jj++) 
        {
            for (int ii=0; ii<N; ii++) 
            {
                double x = minBound[0] + (double)ii*(maxBound[0]-minBound[0])/division[0] + dx[0]/2.; 
                double y = minBound[1] + (double)jj*(maxBound[1]-minBound[1])/division[1] + dx[1]/2.; 
                double z = minBound[2] + (double)kk*(maxBound[2]-minBound[2])/division[2] + dx[2]/2.; 
                ofsmoothed << x << "," << y << "," << z << "," << data_random_smoothed(ii+jj*N+kk*N*N,0) << std::endl;
                ofsmoothed2 << x << "," << y << "," << z << "," << data_random_smoothed2(ii+jj*N+kk*N*N,0) << std::endl;
                ofsmoothed3 << x << "," << y << "," << z << "," << data_random_smoothed3(ii+jj*N+kk*N*N,0) << std::endl;
            }
        }
    }
    ofsmoothed.close();
    ofsmoothed2.close();
    ofsmoothed3.close();


}

void TestGrid_TestCellCenteredData()
{
    Eigen::Vector3d minBound; 
    Eigen::Vector3d maxBound; 
    Eigen::Vector3i division; 
    // test on cell centered data
    {
        int N = 2; 
        minBound << -1,-2,-3; 
        maxBound << 1,2,3;
        division << N,N,N;

        Eigen::Vector3d dx = (maxBound - minBound).cwiseQuotient( division.cast<double>() ); 

        UniformGrid<double> grid( minBound, maxBound, division );

        std::shared_ptr<Eigen::MatrixXd> data0(new Eigen::MatrixXd()); 
        std::shared_ptr<Eigen::MatrixXd> data1(new Eigen::MatrixXd()); 

        data0->resize(N*N*N,1);
        *data1 = Eigen::MatrixXd::Ones(N*N*N,1);

        std::ofstream of("data/raw_cellcentered.csv");
        for (int ii=0; ii<N; ii++) 
        {
            for (int jj=0; jj<N; jj++) 
            {
                for (int kk=0; kk<N; kk++) 
                {
                    (*data0)(ii+jj*N+kk*N*N,0) = (double)(ii+jj*N+kk*N*N) ; 
                    double x = minBound[0] + ii*(maxBound[0]-minBound[0])/division[0] + dx[0]/2.; 
                    double y = minBound[1] + jj*(maxBound[1]-minBound[1])/division[1] + dx[1]/2.; 
                    double z = minBound[2] + kk*(maxBound[2]-minBound[2])/division[2] + dx[2]/2.; 
                    of << x << "," << y << "," << z << "," << (*data0)(ii+jj*N+kk*N*N,0) << std::endl;
                }
            }
        }
        of.close();

        grid.InsertCellCenteredData( "data0", data0 ); 
        grid.InsertCellCenteredData( "data1", data1 );

        int Nquery = 100;
        //dx = (maxBound - minBound)/Nquery;
        std::ofstream of2("data/interpolated_cellcentered.csv");
        for (int ii=0; ii<Nquery; ii++) 
        {
            std::cout << ii << " "; 
            std::cout.flush();
            for (int jj=0; jj<Nquery; jj++) 
            {
                for (int kk=0; kk<Nquery; kk++) 
                {
                    double x = minBound[0] + ii*(maxBound[0]-minBound[0])/Nquery; 
                    double y = minBound[1] + jj*(maxBound[1]-minBound[1])/Nquery; 
                    double z = minBound[2] + kk*(maxBound[2]-minBound[2])/Nquery; 
                    Eigen::Vector3d queryPoint; 
                    queryPoint << x, y, z;
                    double interpData = grid.InterpolateCellCenteredData( "data0", queryPoint )(0);
                    of2 << x << "," << y << "," << z << "," << interpData << std::endl;
                }
            }
        }
        std::cout << std::endl;
        of2.close();
    }

    // test on vertex data
    {
        int N = 3; 
        minBound << -1,-2,-3; 
        maxBound << 1,2,3;
        division << N,N,N;

        Eigen::Vector3d dx = (maxBound - minBound).cwiseQuotient( division.cast<double>() ); 

        UniformGrid<double> grid( minBound, maxBound, division-Eigen::Vector3i::Ones() ); // cell count is vertex count -1 in each direction

        std::shared_ptr<Eigen::MatrixXd> data0(new Eigen::MatrixXd()); 
        std::shared_ptr<Eigen::MatrixXd> data1(new Eigen::MatrixXd()); 

        data0->resize(N*N*N,1);
        *data1 = Eigen::MatrixXd::Ones(N*N*N,1);

        std::ofstream of("data/raw_vertex.csv");
        for (int ii=0; ii<N; ii++) 
        {
            for (int jj=0; jj<N; jj++) 
            {
                for (int kk=0; kk<N; kk++) 
                {
                    (*data0)(ii+jj*N+kk*N*N,0) = (double)(ii+jj*N+kk*N*N) ; 
                    double x = minBound[0] + ii*(maxBound[0]-minBound[0])/division[0] + dx[0]/2.; 
                    double y = minBound[1] + jj*(maxBound[1]-minBound[1])/division[1] + dx[1]/2.; 
                    double z = minBound[2] + kk*(maxBound[2]-minBound[2])/division[2] + dx[2]/2.; 
                    of << x << "," << y << "," << z << "," << (*data0)(ii+jj*N+kk*N*N,0) << std::endl;
                }
            }
        }
        of.close();

        grid.InsertVertexData( "data0", data0 ); 
        grid.InsertVertexData( "data1", data1 );

        int Nquery = 100;
        //dx = (maxBound - minBound)/Nquery;
        std::ofstream of2("data/vertex_interpolated.csv");
        for (int ii=0; ii<Nquery; ii++) 
        {
            std::cout << ii << " "; 
            std::cout.flush();
            for (int jj=0; jj<Nquery; jj++) 
            {
                for (int kk=0; kk<Nquery; kk++) 
                {
                    double x = minBound[0] + ii*(maxBound[0]-minBound[0])/Nquery; 
                    double y = minBound[1] + jj*(maxBound[1]-minBound[1])/Nquery; 
                    double z = minBound[2] + kk*(maxBound[2]-minBound[2])/Nquery; 
                    Eigen::Vector3d queryPoint; 
                    queryPoint << x, y, z;
                    double interpData = grid.InterpolateVertexData( "data0", queryPoint )(0);
                    of2 << x << "," << y << "," << z << "," << interpData << std::endl;
                }
            }
        }
        std::cout << std::endl;
        of2.close();
    }
}

void TestGrid_TestGradient()
{
    Eigen::Vector3d minBound; 
    Eigen::Vector3d maxBound; 
    Eigen::Vector3i division; 
    // test on finite-difference 
    std::cout << "testing for finite-difference\n"; 
    const int N=10; 
    minBound << -1,-1,-1; 
    maxBound <<  1, 1, 1; 
    division <<  N, N, N; 

    Eigen::Vector3d dx = (maxBound - minBound).cwiseQuotient( division.cast<double>() ); 

    UniformGrid<double> grid( minBound, maxBound, division ); 

    std::shared_ptr<Eigen::MatrixXd> data_sine(new Eigen::MatrixXd()); 

    data_sine->resize( N*N*N, 1 ); 

    std::ofstream of("data/raw_cell_centered.csv");
    for (int ii=0; ii<N; ii++) 
    {
        for (int jj=0; jj<N; jj++) 
        {
            for (int kk=0; kk<N; kk++) 
            {
                double x = minBound[0] + ii*(maxBound[0]-minBound[0])/division[0] + dx[0]/2.; 
                double y = minBound[1] + jj*(maxBound[1]-minBound[1])/division[1] + dx[1]/2.; 
                double z = minBound[2] + kk*(maxBound[2]-minBound[2])/division[2] + dx[2]/2.; 
                (*data_sine)(ii+jj*N+kk*N*N,0) = sin( 2.*M_PI*x/(maxBound[0]-minBound[0]) ) +
                                                 sin( 2.*M_PI*y/(maxBound[1]-minBound[1]) ) +
                                                 sin( 2.*M_PI*z/(maxBound[2]-minBound[2]) ); 
                of << x << "," << y << "," << z << "," << (*data_sine)(ii+jj*N+kk*N*N,0) << std::endl;
            }
        }
    }
    of.close();

    grid.InsertCellCenteredData( "data_sine", data_sine ); 

    Eigen::MatrixXd data_sine_gradient   = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::ALL ); 
    Eigen::MatrixXd data_sine_gradientx  = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::X   ); 
    Eigen::MatrixXd data_sine_gradienty  = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::Y   ); 
    Eigen::MatrixXd data_sine_gradientz  = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::Z   ); 
    Eigen::MatrixXd data_sine_gradientxy = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::XY  ); 
    Eigen::MatrixXd data_sine_gradientyz = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::YZ  ); 
    Eigen::MatrixXd data_sine_gradientxz = grid.CellCenteredDataGradient( "data_sine", UniformGrid<double>::XZ  ); 

    std::cout << "checking x gradient components (should be all zeros): " << std::endl;
    std::cout << (data_sine_gradient.col(0) + data_sine_gradientxy.col(0) + data_sine_gradientxz.col(0) - data_sine_gradientx*3).topRows(10) << std::endl;
    std::cout << "checking y gradient components (should be all zeros): " << std::endl;
    std::cout << (data_sine_gradient.col(1) + data_sine_gradientxy.col(1) + data_sine_gradientyz.col(0) - data_sine_gradienty*3).topRows(10) << std::endl;
    std::cout << "checking z gradient components (should be all zeros): " << std::endl;
    std::cout << (data_sine_gradient.col(2) + data_sine_gradientyz.col(1) + data_sine_gradientxz.col(1) - data_sine_gradientz*3).topRows(10) << std::endl;

    std::ofstream ofx("data/cell_centered_derivative_x.csv");
    std::ofstream ofy("data/cell_centered_derivative_y.csv");
    std::ofstream ofz("data/cell_centered_derivative_z.csv");
    for (int ii=0; ii<N; ii++) 
    {
        for (int jj=0; jj<N; jj++) 
        {
            for (int kk=0; kk<N; kk++) 
            {
                double x = minBound[0] + (double)ii*(maxBound[0]-minBound[0])/division[0] + dx[0]/2.; 
                double y = minBound[1] + (double)jj*(maxBound[1]-minBound[1])/division[1] + dx[1]/2.; 
                double z = minBound[2] + (double)kk*(maxBound[2]-minBound[2])/division[2] + dx[2]/2.; 
                ofx << x << "," << y << "," << z << "," << data_sine_gradient(ii+jj*N+kk*N*N,0) << std::endl;
                ofy << x << "," << y << "," << z << "," << data_sine_gradient(ii+jj*N+kk*N*N,1) << std::endl;
                ofz << x << "," << y << "," << z << "," << data_sine_gradient(ii+jj*N+kk*N*N,2) << std::endl;
            }
        }
    }
    ofx.close();
    ofy.close();
    ofz.close();




}

void TestGrid_RotatedGrid()
{
    Eigen::Vector3d minBound(-1,-2,-3); 
    Eigen::Vector3d maxBound( 1, 2, 3); 
    Eigen::Vector3i cellCount(20,20,20); 

    UniformGrid<double> uniform(minBound,maxBound,cellCount); 
    RotatedUniformGrid<double> rotated(minBound,maxBound,cellCount); 
    rotated.SetYawDegree(30); 

    uniform.WriteVTKCellGrid("data/uniform_grid.vtk"); 
    rotated.WriteVTKCellGrid("data/rotated_grid.vtk"); 


}

void TestGrid()
{
    std::cout << "testing for grid" << std::endl;

    TestOpenfoamMesh();
    //TestGrid_TestCellCenteredData(); 
    //TestGrid_TestGradient();
    //TestGrid_TestSmoothing(); 
    //TestGrid_RotatedGrid(); 
}
