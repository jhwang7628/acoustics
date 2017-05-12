#ifndef BOUNDARY_INTERFACE_H 
#define BOUNDARY_INTERFACE_H 
#include <vector>
#include <Eigen/Dense> 
#include <memory>
#include "config.h"
#include "linearalgebra/Vector3.hpp" 
//##############################################################################
// Class BoundaryInterface
//##############################################################################
class BoundaryInterface
{
private: 
    Tuple3i         _divisions; // [N1, N2, N3]
    Eigen::MatrixXd _positions; // (N1*N2*N3)-by-3
    Eigen::MatrixXd _values;    // (N1*N2*N3)-by-t, t is #timesteps read
public: 
    static std::shared_ptr<BoundaryInterface> Load(const std::string &filename); 
    REAL Evaluate(const Tuple3i &ind) const; 
};
using BoundaryInterfacePtr = std::shared_ptr<BoundaryInterface>; 
#endif
