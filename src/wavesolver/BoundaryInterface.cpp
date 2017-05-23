#include "utils/IO/IO.h" 
#include "wavesolver/BoundaryInterface.h"
//##############################################################################
// Function Load
//##############################################################################
BoundaryInterfacePtr BoundaryInterface:: 
Load(const std::string &prefix) 
{
    const std::string pos_file = prefix + "_listening_position.dat"; 
    const std::string dat_file = prefix + "_listened_data.dat"; 
    BoundaryInterfacePtr interface = std::make_shared<BoundaryInterface>(); 
    IO::readMatrixX<double>(interface->Positions(), pos_file.c_str(), IO::BINARY, 1); 
    IO::readMatrixX<double>(interface->Values(), dat_file.c_str(), IO::BINARY, 1); 
    return interface; 
}

//##############################################################################
// Function Evaluate
//##############################################################################
REAL BoundaryInterface::
Evaluate(const Vector3d &pos) const
{
    if (_values.rows()==0 || _values.cols()==0)
        return 0.0; 
    // search for nearest neighbour
    REAL mindist = std::numeric_limits<REAL>::max(); 
    int minind = -1; 
    for (int r=0; r<_positions.rows(); ++r)
    {
        const REAL dist = pow(pos[0] - _positions(r,0), 2) 
                        + pow(pos[1] - _positions(r,1), 2) 
                        + pow(pos[2] - _positions(r,2), 2); 
        if (dist<mindist && dist<POS_TOL)
        {
            mindist = dist; 
            minind = r; 
        }
    }
    return (minind >= 0 ? _values(pointer, minind) : 0.0);
}
