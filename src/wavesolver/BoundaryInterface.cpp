#include "TYPES.h"
#include "utils/IO/IO.h" 
#include "wavesolver/BoundaryInterface.h"
//##############################################################################
// Global variables
//##############################################################################
REAL BoundaryInterface::_blendTotalTime = 0.01;

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
        VECTOR vPressure; 
        auto &unit_a = (_WhichUnit(solver_a) == 0 ? 
                _simUnitPair.first ->simulator 
              : _simUnitPair.second->simulator); 
        auto &unit_b = (_WhichUnit(solver_a) == 0 ? 
                _simUnitPair.second->simulator
              : _simUnitPair.first ->simulator); 
        const Tuple3i cellIndices = 
            unit_b->GetGrid().pressureField().cellIndex(cell_b); 
        unit_b->GetSolver()->vertexPressure(cellIndices, vPressure);
        pressure_b = vPressure[0];
    }
    return keyExists;
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
