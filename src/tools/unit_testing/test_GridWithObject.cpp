#include "grid/GridWithObject.h" 

int main()
{

    Eigen::Vector3i cellCount(250,250,250); 
    Eigen::Vector3d minBound(-0.2178337591671500,-0.2017749603171500,-0.2192579190471500);
    Eigen::Vector3d maxBound( 0.2204221281671500, 0.2364809270171500, 0.2189979682871500); 

    UniformGridWithObject grid(minBound,maxBound,cellCount); 

    grid.Reinitialize("/home/jui-hsien/code/acoustics/work/impulse-response/config/default.xml");
    grid.ClassifyCells();
    grid.WriteCellTypes("grid_types.csv");


    return 0; 
}



