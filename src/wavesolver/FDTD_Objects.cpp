#include <wavesolver/FDTD_Objects.h> 

//##############################################################################
// Handles a list of objects embedded in the wave solver
// note this is not thread safe
//##############################################################################
void FDTD_Objects::
AddObject(const int &objectName, RigidObjectPtr &object)
{
    _rigidObjects[objectName] = object; 
}

//##############################################################################
//##############################################################################
int FDTD_Objects::
AddConstraint(const std::string &id, FDTD_PlaneConstraint_Ptr &constraint)
{
    _constraints[id] = constraint; 
    return _constraints.size();
}

//##############################################################################
//##############################################################################
int FDTD_Objects::
AddConstraints(FDTD_Objects_Ptr &rhs)
{
    auto &c = rhs->GetConstraints(); 
    _constraints.insert(c.begin(), c.end()); 
    return _constraints.size();
}

//##############################################################################
//##############################################################################
int FDTD_Objects::
OccupyByObject(const Vector3d &positionWorld) 
{
    const int N_objects = N(); 
    for (const auto &m : _rigidObjects)
    {
        const double distance = m.second->DistanceToMesh(positionWorld.x, positionWorld.y, positionWorld.z); 
        if (distance < DISTANCE_TOLERANCE) 
            return m.first; 
    }
    return -1;
}

//##############################################################################
//##############################################################################
bool FDTD_Objects::
OccupyByConstraint(const Vector3d &pos)
{
    FDTD_PlaneConstraint_Ptr ptr;
    return OccupyByConstraint(pos, ptr);
}

//##############################################################################
//##############################################################################
bool FDTD_Objects::
OccupyByConstraint(const Vector3d &pos, FDTD_PlaneConstraint_Ptr &constraint)
{
    for (const auto &p : _constraints)
        if (p.second->Inside(pos))
        {
            constraint = p.second; 
            return true;
        }
    return false; 
}

//##############################################################################
//##############################################################################
REAL FDTD_Objects::
ObjectDistance(const int &objectIndex, const Vector3d &positionWorld) 
{
    return _rigidObjects.at(objectIndex)->DistanceToMesh(positionWorld.x,positionWorld.y,positionWorld.z);
}

//##############################################################################
//##############################################################################
REAL FDTD_Objects::
LowestObjectDistance(const Vector3d &positionWorld) 
{
    REAL distance = std::numeric_limits<REAL>::max(); 
    for (const auto &m : _rigidObjects) 
        distance = std::min<REAL>(distance, m.second->DistanceToMesh(positionWorld.x, positionWorld.y, positionWorld.z)); 
    return distance; 
}

//##############################################################################
// Find the closest object with respect to positionWorld
//
// Note: this might fail if positionWorld is too far from the region where sdf
// is defined. 
//##############################################################################
bool FDTD_Objects::
LowestObjectDistance(const Vector3d &positionWorld, REAL &distance, int &objectID) 
{
    REAL queriedDistance; 
    bool found = false; 
    distance = std::numeric_limits<REAL>::max(); 
    objectID = std::numeric_limits<int>::max(); 
    for (const auto &m : _rigidObjects) 
    {
        queriedDistance = m.second->DistanceToMesh(positionWorld.x, positionWorld.y, positionWorld.z); 
        if (queriedDistance < distance)
        {
            distance = queriedDistance; 
            objectID = m.first;
            found = true; 
        }
    } 
    return found; 
}

//##############################################################################
//##############################################################################
bool FDTD_Objects::
LowestConstraintDistance(const Vector3d &position, REAL &unsignedDistance,
                         FDTD_PlaneConstraint_Ptr &constraint)
{
    REAL d; 
    bool found = false; 
    unsignedDistance = std::numeric_limits<REAL>::max(); 
    constraint.reset(); 
    for (const auto &m : _constraints) 
    {
        d = m.second->UnsignedDistanceToPlane(position); 
        if (d < unsignedDistance)
        {
            unsignedDistance = d; 
            constraint = m.second;
            found = true; 
        }
    } 
    return found; 
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
ObjectNormal(const int &objectIndex, const Vector3d &positionWorld, Vector3d &queriedNormal) 
{
    _rigidObjects.at(objectIndex)->NormalToMesh(positionWorld.x,positionWorld.y,positionWorld.z, queriedNormal);
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

////##############################################################################
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
//
//
//
// THIS FUNCTION IS NO LONGER USED
//##############################################################################
bool FDTD_Objects::
ReflectAgainstAllBoundaries(const int &startObjectID, const Vector3d &originalPoint, const REAL &time, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &accumulatedBoundaryConditionValue, const REAL &density, const int &maxIteration)
{
    assert(startObjectID >= 0);

    // reflection for fixed amount of times, can break early if any
    // intermediate point is outside of all objects (not first iteration
    // though)
    REAL distance; 
    int objectID = startObjectID;
    accumulatedBoundaryConditionValue = 0.0;
    Vector3d intermediatePoint = originalPoint; 
    for (int ii=0; ii<maxIteration; ++ii)
    {
        _rigidObjects.at(objectID)->ReflectAgainstBoundary(intermediatePoint, reflectedPoint, boundaryPoint, erectedNormal, distance); 
        if (distance < 0) // originalPoint was inside the object
            accumulatedBoundaryConditionValue += -2.0*distance*_rigidObjects.at(objectID)->EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, time) *density;
        else if (distance >= 0) // originalPoint was outside the object
            accumulatedBoundaryConditionValue += distance*_rigidObjects.at(objectID)->EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, time) *density;

        // if new reflectedPoint is outside object, break the loop
        objectID = OccupyByObject(reflectedPoint); 
        if (objectID == -1)
            return true; 
        else
            intermediatePoint = reflectedPoint; 
    }
    return false; 
}

//##############################################################################
// @return ODE time (should be the same across all ODEs registered)
//##############################################################################
REAL FDTD_Objects::
AdvanceAllModalODESolvers(const int &N_steps)
{
    REAL t = -1.0;
    for (auto &object : _rigidObjects)
        if (object.second->Type() == RIGID_SOUND_OBJ)
        {
            t = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object.second)
                ->AdvanceModalODESolvers(N_steps); 
        }
    return t; 
}

//##############################################################################
//##############################################################################
REAL FDTD_Objects::
GetEarliestImpactEvent()
{
    REAL earliestTime = std::numeric_limits<REAL>::max(); 
    for (const auto &object : _rigidObjects)
    {
        if (object.second->Type() == RIGID_SOUND_OBJ)
        {
            auto o = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object.second); 
            if (o->N_Impulses() > 0)
                earliestTime = min(earliestTime, o->GetFirstImpulseTime());
        }
    }
    return earliestTime; 
}

//##############################################################################
// Function SetObjectStates
//##############################################################################
void FDTD_Objects::
SetObjectStates(const REAL time)
{
    // update modal vectors for the next time step
    for (auto &m : _rigidObjects) 
        if (m.second->Type() == RIGID_SOUND_OBJ) 
        {
            std::dynamic_pointer_cast<FDTD_RigidSoundObject>(m.second)
                ->UpdateQPointers(); 
        }
    AnimateObjects(time);
}

//##############################################################################
// Function InitializeAnimator
//##############################################################################
void FDTD_Objects::
InitializeAnimator(const std::string &fileDisplacement,
                   const std::string &fileVelocity,
                   const std::string &fileAcceleration)
{
    _objectAnimator = std::make_shared<FDTD_RigidObject_Animator>(); 
    _objectAnimator->ReadAllKinematics(fileDisplacement, 
                                       fileVelocity, 
                                       fileAcceleration); 
    AnimateObjects(0.0); // apply the transformation right away
}

//##############################################################################
// Function AnimateObjects
//##############################################################################
void FDTD_Objects::
AnimateObjects(const REAL toTime)
{
    if (!_objectAnimator)
        return; 
    Point3d newCOM; 
    Quaternion<REAL> quaternion; 
    for (auto &m : _rigidObjects)
    {
        const auto &object = m.second;
        if (object->Animated())
        {
            _objectAnimator->GetRigidObjectTransform(m.first, toTime, newCOM, quaternion);
            quaternion.normalize();
            object->SetRigidBodyTransform(newCOM, quaternion);
        }
    }
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
Join(FDTD_Objects_Ptr &toBeJoined) 
{
    auto &to_be_joined_obj_list = toBeJoined->GetRigidObjects(); 
    auto &to_be_joined_src_list = toBeJoined->GetPressureSources(); 
    for (auto &m : to_be_joined_obj_list)
    {
        _rigidObjects[m.first] = m.second; 
    }
    for (auto &p : to_be_joined_src_list)
    {
        _pressureSources.push_back(p); 
    }
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
void FDTD_Objects::
WriteFailedReflections(const std::string &filename)
{
    for (auto &object : _rigidObjects) 
    {
        const std::string objFilename = filename + object.second->GetMeshName();
        object.second->WriteDebugArrow(objFilename); 
    }
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
ClearFailedReflections()
{
    for (auto &object : _rigidObjects) 
        object.second->ClearDebugArrow(); 
}

//##############################################################################
//##############################################################################
void FDTD_Objects::
DebugWriteModalQ(const int steps, const std::string &filename)
{
    std::ofstream tmp("tmp.dat", std::ostream::out | std::ostream::binary);
    if (!of_q) 
    {
        of_q = new std::ofstream(filename.c_str(), 
                                 std::ostream::out | std::ostream::binary); 
    }
    for (auto &object : _rigidObjects)
    {
        if (object.second->Type() == RIGID_SOUND_OBJ &&
            std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object.second)->IsModalObject())
            std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object.second)->AdvanceModalODESolvers(steps, 0, tmp, *of_q); 
    }
    tmp.close();
    of_q->close();
}

//##############################################################################
//##############################################################################
std::ostream &operator <<(std::ostream &os, const FDTD_Objects &objects) 
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Class FDTD_Objects\n" 
       << "--------------------------------------------------------------------------------\n";
    for (const auto &m : objects._rigidObjects)
    {
        os << " Object <" << m.first << ", " << m.second->GetMeshName() << "\n"; 
    }
    os << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}

