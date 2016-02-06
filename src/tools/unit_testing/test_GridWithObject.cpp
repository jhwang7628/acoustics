#include <Eigen/Dense>
#include "grid/GridWithObject.h" 
#include "utils/IO/IO.h"
#include "vtkConverter/vtkConverter.h" 

int main()
{

    UniformGridWithObject grid; 
    grid.Reinitialize("data_test/default.xml");
    //grid.ClassifyCells();
    //grid.WriteCellTypes("grid_types",1);
    //grid.Test_Reflection(); 
    //grid.ComputeInterfaceStencils();


    std::shared_ptr<Eigen::MatrixXd> data(new Eigen::MatrixXd()); 
    IO::readMatrixX<double>(*data, "data_test/test_pressure_00030.dat", IO::BINARY); 

    grid.InsertCellCenteredData("test_data", data); 
    //grid.WriteVTKCellCentered("test_pressure_00030", "test_data", "test_data");


    const Eigen::Vector3i cellCount = grid.GetCellCount(); 
    const int N_cells = cellCount[0]*cellCount[1]*cellCount[2];
    const int N_hessian=6; 
    std::vector<std::shared_ptr<Eigen::MatrixXd>> hessian(N_hessian); 
    for (int ii=0; ii<N_hessian; ii++) 
        hessian[ii].reset(new Eigen::MatrixXd(N_cells, 1)); 
    grid.CellCenteredScalarHessian("test_data", hessian); 

    for (int ii=0; ii<N_hessian; ii++) 
    {
        std::string name("ddG_"+std::to_string(ii)); 
        grid.InsertCellCenteredData(name, hessian[ii]); 
        grid.WriteVTKCellCentered(name,name,"ddG"); 
    }







    return 0; 
}



