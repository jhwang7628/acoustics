#ifndef BOUNDARY_INTERFACE_H 
#define BOUNDARY_INTERFACE_H 
#include <vector>
#include <Eigen/Dense> 
#include <memory>
#include "config.h"
#include "linearalgebra/Vector3.hpp" 
#include "wavesolver/SimWorldAuxDS.h"

//##############################################################################
// Class BoundaryInterface
//##############################################################################
class BoundaryInterface
{
private: 
    std::pair<ActiveSimUnit_Ptr, ActiveSimUnit_Ptr> _simUnitPair; 
    std::vector<std::pair<int,int>>                 _cellPairs;
    static REAL                                     _blendTotalTime; 
    REAL                                            _blendStartTime;
    int                                             _direction; //0:x, 1:y, 2:z

    int _WhichUnit(const std::string &id) const; 

public: 
    bool initialized = false; 

    BoundaryInterface(ActiveSimUnit_Ptr unit_a,
                      ActiveSimUnit_Ptr unit_b, 
                      const REAL &time,
                      const int &direction)
        : _simUnitPair(std::make_pair(unit_a, unit_b)),
          _blendStartTime(time),
          _direction(direction)
    {}
    bool Equal(const BoundaryInterface &rhs) const; 
    std::string GetOtherSolver(const std::string &solver_a) const; 
    bool GetOtherCell(const std::string &solver_a, 
                      const int &cell_a,
                            int &cell_b) const; 
    bool GetOtherCellPressure(const std::string &solver_a, 
                              const int &cell_a,
                                    REAL &pressure_b) const; 
    REAL GetBlendCoeff(const REAL &time) const; 
    ActiveSimUnit_Ptr GetSimUnit(const std::string &id) const;
    inline void AddCellPair(std::pair<int,int> &pair)
    {_cellPairs.push_back(std::move(pair));}
    inline auto GetCellPairsCount()
    {return _cellPairs.size();}
    inline auto GetDirection() const
    {return _direction;} 
    inline auto GetSimUnit_a()
    {return _simUnitPair.first;}
    inline auto GetSimUnit_b()
    {return _simUnitPair.second;}
    inline const auto GetSimUnit_a() const
    {return _simUnitPair.first;}
    inline const auto GetSimUnit_b() const
    {return _simUnitPair.second;}

    friend
    std::ostream &operator <<(std::ostream &os, 
                              const BoundaryInterface &interface); 
};
#endif
