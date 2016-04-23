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
int FDTD_Objects::
OccupyByObject(const Vector3d &positionWorld) 
{
    const int N_objects = N(); 
    for (int ii=0; ii<N_objects; ii++) 
    {
        const double distance = _rigidObjects[ii]->DistanceToMesh(positionWorld.x, positionWorld.y, positionWorld.z); 
        if (distance < DISTANCE_TOLERANCE) 
            return ii; 
    }
    return -1;
}

//##############################################################################
//##############################################################################
REAL FDTD_Objects::
ObjectDistance(const int &objectIndex, const Vector3d &positionWorld) 
{
    return _rigidObjects[objectIndex]->DistanceToMesh(positionWorld.x,positionWorld.y,positionWorld.z);
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
ObjectNormal(const int &objectIndex, const Vector3d &positionWorld, Vector3d &queriedNormal) 
{
    _rigidObjects[objectIndex]->NormalToMesh(positionWorld.x,positionWorld.y,positionWorld.z, queriedNormal);
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
AddVibrationalSourceToObject(const int &objectIndex, VibrationalSourcePtr &sourcePtr) 
{
    _rigidObjects.at(objectIndex)->AddVibrationalSource(sourcePtr);
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
TestObjectDistanceField(const size_t &ind)
{
    _rigidObjects.at(ind)->TestQueryDistance();
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

