#include "wavesolver/BoundaryInterface.h"
//##############################################################################
// Function Load
//##############################################################################
BoundaryInterfacePtr BoundaryInterface:: 
Load(const std::string &filename) 
{
    std::cout << "load boundary interfaces\n";
    return std::make_shared<BoundaryInterface>(); 
}

//##############################################################################
// Function Evaluate
//##############################################################################
REAL BoundaryInterface::
Evaluate(const Tuple3i &ind) const
{
    std::cout << "evaluate boundary interfaces with ind=\n"; 
    std::cout << ind << std::endl;
    return 1.0; 
}
