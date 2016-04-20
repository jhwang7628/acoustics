#include <wavesolver/FDTD_Objects.h> 

//##############################################################################
// Handles a list of objects embedded in the wave solver
// note this is not thread safe
//##############################################################################
void FDTD_Objects::
AddObject(const std::string &objectName, RigidObjectPtr &object)
{
    if (_meshIDMap.find(objectName) == _meshIDMap.end())
    {
        _rigidObjects.push_back(object); 
        _meshIDMap[objectName] = _rigidObjects.size()-1; 
    }
}

//##############################################################################
//##############################################################################
std::ostream &operator<<(std::ostream &os, const FDTD_Objects &objects) 
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Class FDTD_Objects\n" 
       << "--------------------------------------------------------------------------------\n";
    for (int ii=0; ii<objects.N(); ++ii) 
    {
        const std::string &meshName = objects.GetMeshName(ii); 
        os << " Object " << ii << ": " << meshName
           << " <ID:" << objects.GetMeshID(meshName) << ">\n"; 
    }
    os << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
