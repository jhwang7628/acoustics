#include "TYPES.h"
#include "utils/IO/IO.h" 
#include "wavesolver/BoundaryInterface.h"
//##############################################################################
// Global variables
//##############################################################################
REAL BoundaryInterface::_blendTotalTime = 0.0001;

//##############################################################################
// Function _WhichUnit
//##############################################################################
int BoundaryInterface::
_WhichUnit(const std::string &id) const
{
    const auto &unit_a = _simUnitPair.first; 
    const auto &unit_b = _simUnitPair.second; 
    if (*(unit_a->simulator->GetSimulatorID()) == id)
        return 0; 
    else if (*(unit_b->simulator->GetSimulatorID()) == id) 
        return 1; 
    return -1; 
}

//##############################################################################
// Function GetSimUnit
//##############################################################################
ActiveSimUnit_Ptr BoundaryInterface::
GetSimUnit(const std::string &id) const
{
    auto unit_a = _simUnitPair.first; 
    auto unit_b = _simUnitPair.second; 
    if (*(unit_a->simulator->GetSimulatorID()) == id)
        return unit_a; 
    else if (*(unit_b->simulator->GetSimulatorID()) == id) 
        return unit_b; 
    return nullptr; 
}
//##############################################################################
// Function Equal
//##############################################################################
bool BoundaryInterface::
Equal(const BoundaryInterface &rhs) const
{
    if (_simUnitPair.first && _simUnitPair.second)
    {
        // get ids for this
        const auto &unit_a = _simUnitPair.first; 
        const auto &unit_b = _simUnitPair.second; 
        const std::string *id_a = unit_a->simulator->GetSimulatorID(); 
        const std::string *id_b = unit_b->simulator->GetSimulatorID(); 
        if (!id_a || !id_b) 
            return false; 

        // get ids for rhs
        const auto &unit_a_rhs = rhs.GetSimUnit_a(); 
        const auto &unit_b_rhs = rhs.GetSimUnit_b(); 
        const std::string *id_a_rhs = unit_a_rhs->simulator->GetSimulatorID(); 
        const std::string *id_b_rhs = unit_b_rhs->simulator->GetSimulatorID(); 
        if (!id_a_rhs || !id_b_rhs) 
            return false; 

        if ((((*id_a) == (*id_a_rhs) && (*id_b) == (*id_b_rhs)) ||
             ((*id_a) == (*id_b_rhs) && (*id_b) == (*id_a_rhs))) && 
            (_direction == rhs.GetDirection()))
            return true; 
    }
    return false; 
}

//##############################################################################
// Function GetOtherSolver
//##############################################################################
std::string BoundaryInterface:: 
GetOtherSolver(const std::string &solver_a) const
{
    assert(_simUnitPair.first && _simUnitPair.second); 
    const int u = _WhichUnit(solver_a); 
    if (u == 0)      return *(_simUnitPair.second->simulator->GetSimulatorID());
    else if (u == 1) return *(_simUnitPair.first->simulator->GetSimulatorID());
    return ""; 
}

//##############################################################################
// Function GetOtherCell
//##############################################################################
bool BoundaryInterface::
GetOtherCell(const std::string &solver_a, 
             const int &cell_a,
                   int &cell_b) const
{
    const int u = _WhichUnit(solver_a); 
    bool keyExists = false; 
    for (const auto &pair : _cellPairs)
    {
        if      (u==0 && pair.first==cell_a)
        {
            keyExists = true; 
            cell_b = pair.second; 
            break;
        }
        else if (u==1 && pair.second==cell_a)
        {
            keyExists = true; 
            cell_b = pair.first; 
            break; 
        }
    }
    return keyExists; 
}

//##############################################################################
// Function GetOtherCellPressure
//##############################################################################
bool BoundaryInterface::
GetOtherCellPressure(const std::string &solver_a, 
                     const int &cell_a,
                           REAL &pressure_b) const
{
    int cell_b; 
    const bool keyExists = GetOtherCell(solver_a, cell_a, cell_b); 
    if (keyExists)
    {
        auto &unit_a = (_WhichUnit(solver_a) == 0 ? 
                _simUnitPair.first ->simulator 
              : _simUnitPair.second->simulator); 
        auto &unit_b = (_WhichUnit(solver_a) == 0 ? 
                _simUnitPair.second->simulator
              : _simUnitPair.first ->simulator); 
        auto &grid_b = unit_b->GetGrid(); 
        const Tuple3i cellIndices = 
            grid_b.pressureField().cellIndex(cell_b); 
        // TODO START
        // check if its bulk or ghost cells .. 
        // also need to modify classify cells routine to take the neighbours into account
        //
        //
        // TODO END
        if (grid_b.IsPressureCellBulk(cell_b))
        {
            VECTOR vPressure; 
            unit_b->GetSolver()->vertexPressure(cellIndices, vPressure);
            pressure_b = vPressure[0];
        }
        else if (grid_b.IsPressureCellGhost(cell_b))
        {
            grid_b.BoundaryGhostCellPressure(solver_a, cell_a, cell_b, pressure_b); 
        }
        else 
        {
            pressure_b = 0.0;
        }
    }
    return keyExists;
}
//##############################################################################
// Function GetBlendCoeff
//##############################################################################
REAL BoundaryInterface::GetBlendCoeff(const REAL &time) const
{
    if (time < _blendStartTime)
        return 0.0; 
    else if (time > (_blendStartTime + _blendTotalTime))
        return 1.0; 
    return (time - _blendStartTime)/_blendTotalTime;
}


//##############################################################################
//##############################################################################
std::ostream &operator <<(std::ostream &os, 
                          const BoundaryInterface &interface)
{
    os << "------------------------------------------------------------------\n" 
       << "Class BoundaryInterface\n" 
       << "------------------------------------------------------------------\n"
       << " num cell pairs   : " << interface._cellPairs.size() << "\n"
       << " direction        : " << interface._direction << "\n"
       << " blend total time : " << BoundaryInterface::_blendTotalTime << "\n"
       << " blend start time : " << interface._blendStartTime << "\n"
       << "------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
