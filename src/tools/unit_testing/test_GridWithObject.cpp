#include "grid/GridWithObject.h" 

int main()
{

    UniformGridWithObject grid; 
    grid.Reinitialize("/home/jui-hsien/code/acoustics/work/impulse-response/config/default.xml");
    //grid.ClassifyCells();
    grid.WriteCellTypes("grid_types",1);
    //grid.Test_Reflection(); 
    //grid.ComputeInterfaceStencils();


    return 0; 
}



