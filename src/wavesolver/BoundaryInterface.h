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

    static constexpr REAL POS_TOL = 0.01;
public: 
    static std::shared_ptr<BoundaryInterface> Load(const std::string &prefix); 
    inline Eigen::MatrixXd &Positions(){return _positions;}
    inline Eigen::MatrixXd &Values(){return _values;} 
    REAL Evaluate(const Vector3d &pos) const; 

    // TODO use a pointer for now, should fix after actual box management
    int pointer = 0; 
    void AdvancePointer()
    {if (_values.cols()!=0) pointer = (pointer+1)%(int)_values.cols();}
};
using BoundaryInterfacePtr = std::shared_ptr<BoundaryInterface>; 
#endif
