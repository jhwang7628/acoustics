#include "grid/GridWithObject.h" 

int main()
{

    UniformGridWithObject<double> grid; 

    grid.Reinitialize("/home/jui-hsien/code/acoustics/work/impulse-response/config/wavesolver_setup.xml");

    return 0; 
}



