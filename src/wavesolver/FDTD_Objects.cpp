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
REAL FDTD_Objects::
LowestObjectDistance(const Vector3d &positionWorld) 
{
    REAL distance = std::numeric_limits<REAL>::max(); 
    for (int ii=0; ii<N(); ++ii) 
        distance = std::min<REAL>(distance, _rigidObjects[ii]->DistanceToMesh(positionWorld.x, positionWorld.y, positionWorld.z)); 
    return distance; 
}

//##############################################################################
// Find the closest object with respect to positionWorld
//##############################################################################
void FDTD_Objects::
LowestObjectDistance(const Vector3d &positionWorld, REAL &distance, int &objectID) 
{
    REAL queriedDistance; 
    distance = std::numeric_limits<REAL>::max(); 
    objectID = std::numeric_limits<int>::max(); 
    const int N_objects = N(); 
    for (int ii=0; ii<N_objects; ++ii) 
    {
        queriedDistance = _rigidObjects[ii]->DistanceToMesh(positionWorld.x, positionWorld.y, positionWorld.z); 
        if (queriedDistance < distance)
        {
            distance = queriedDistance; 
            objectID = ii;
        }
    } 
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
AddPressureSource(PressureSourcePtr &sourcePtr) 
{
    _pressureSources.push_back(std::move(sourcePtr)); 
}

//##############################################################################
//##############################################################################
REAL FDTD_Objects::
EvaluateNearestVibrationalSources(const Vector3d &position, const Vector3d &normal, const REAL &time) 
{
    REAL distance; 
    int objectID; 
    LowestObjectDistance(position, distance, objectID); 
    return _rigidObjects.at(objectID)->EvaluateBoundaryAcceleration(position, normal, time);
}

//##############################################################################
//##############################################################################
REAL FDTD_Objects::
EvaluatePressureSources(const Vector3d &position, const Vector3d &normal, const REAL &time) 
{
    const int N_sources = _pressureSources.size(); 
    REAL sourceValue = 0; 
    for (int ii=0; ii<N_sources; ++ii) 
    {
        sourceValue += _pressureSources.at(ii)->Evaluate(position, normal, time); 
    }
    return sourceValue; 
}

//##############################################################################
// Compute the reflected stencils until outside all objects. 
// 
// Suppose k reflections are needed to get the stencil, whose original position
// is at originalPoint, outside all objects in the scene, then from the
// discretization, we have 
//
//      1st reflection: P_IP^(1) - P_GC     = h_1 * dp/dn|_BI^(1) 
//      2nd reflection: P_IP^(2) - P_IP^(1) = h_2 * dp/dn|_BI^(2) 
//      3rd reflection: P_IP^(3) - P_IP^(2) = h_3 * dp/dn|_BI^(3) 
//                              .
//                              .
//                              .
//      kth reflection: P_IP^(k) - P_IP^(k-1) = h_k * dp/dn|_BI^(k) 
//
//          GC: original ghost cell
//          IP: image point
//          BI: boundary point
//          h_i: distance travelled at i-th reflection
//          dp/dn|_x: Neumann boundary condition evaluated at position x
//
// Summing all these up, we get
//      P_IP^(k) - P_GC = sum_{i=1}^{k} h_i*dp/dn|_BI^(i)
//
// Return values: 
//  bool: successfully find a point outside all objects
//  reflectedPoint: reflected point found at last reflection 
//  boundaryPoint: boundary point found at last reflection
//  erectedNormal: erected normal found at last reflection
//  accumulatedBoundaryConditionValue: sum of boundary condition multiplied 
//                                     by travelled distance for all steps
//##############################################################################
bool FDTD_Objects::
ReflectAgainstAllBoundaries(const Vector3d &originalPoint, const REAL &time, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &accumulatedBoundaryConditionValue, const REAL &density, const int &maxIteration)
{
    int objectID = OccupyByObject(originalPoint);
    accumulatedBoundaryConditionValue = 0.0;
    if (objectID == -1)  // the input point is not inside any objects in the scene
    {
        reflectedPoint = originalPoint; 
        return true; 
    }

    Vector3d intermediatePoint; 
    REAL distanceTravelled; 
    intermediatePoint = originalPoint; 
    for (int ii=0; ii<maxIteration; ++ii)
    {
        _rigidObjects[objectID]->ReflectAgainstBoundary(intermediatePoint, reflectedPoint, boundaryPoint, erectedNormal, distanceTravelled); 
        accumulatedBoundaryConditionValue += -distanceTravelled*_rigidObjects[objectID]->EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, time) *(-density);
        objectID = OccupyByObject(reflectedPoint); 
        if (objectID == -1)
            return true; 
        intermediatePoint = reflectedPoint; 
    }
    return false; 
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
void FDTD_Objects::
PrintAllSources(std::ostream &os)
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Function FDTD_Objects::PrintAllSources \n" 
       << "--------------------------------------------------------------------------------\n"
       << " number of sources = " << N_sources() << "\n"
       << " source list BEGIN\n"; 
    for (const auto &p : _pressureSources) 
        p->PrintSourceInfo(os);
    os << " source list END\n"
       << "--------------------------------------------------------------------------------\n"; 
    os << std::flush;
}

//##############################################################################
//##############################################################################
std::ostream &operator <<(std::ostream &os, const FDTD_Objects &objects) 
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

