#include <wavesolver/FDTD_RigidObject_Animator.h>
#include <utils/STL_Wrapper.h>

//##############################################################################
//##############################################################################
void FDTD_RigidObject_Animator::
ReadDisplacement(const std::string &filename)
{
    _rigidsimResults = std::make_shared<RigidObjDispReader>(); 
    _rigidsimResults->ReadDisplacement(filename); 
    const auto &timeStamps = _rigidsimResults->Timesteps(); 
    _timeStart = timeStamps.at(0); 
    _timeStop = timeStamps.at(timeStamps.size()-1); 
    _timestep = timeStamps.at(1) - timeStamps.at(0);
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject_Animator::
ReadAllKinematics(const std::string &fileDisplacement, const std::string &fileVelocity, const std::string &fileAcceleration)
{
    _rigidsimResults = std::make_shared<RigidObjDispReader>(); 
    _rigidsimResults->ReadAllKinematics(fileDisplacement, fileVelocity, fileAcceleration); 
    const auto &timeStamps = _rigidsimResults->Timesteps(); 
    _timeStart = timeStamps.at(0); 
    _timeStop = timeStamps.at(timeStamps.size()-1); 
    _timestep = timeStamps.at(1) - timeStamps.at(0);
}

//##############################################################################
// Because of how it is stored, the current implementation might not have good 
// data locality. Might want to consider refactor some parts for a much larger
// scene.
//##############################################################################
void FDTD_RigidObject_Animator::
GetObjectDisplacement(const int &objectID, const REAL &time, Vector3d &displacement, Quaternion<REAL> &quaternion)
{
    const int N_frames = _rigidsimResults->N_Frames(); 
    const auto &timeStamps = _rigidsimResults->Timesteps(); 
    const std::vector<std::vector<Point3<REAL> > > &positions  = _rigidsimResults->Positions();
    const std::vector<std::vector<Quaternion<REAL> > > &quaternions  = _rigidsimResults->Rotations();
    if (time <= _timeStart)
    {
        displacement = positions.at(0).at(objectID); 
        quaternion = quaternions.at(0).at(objectID); 
    }
    else if (time >= _timeStop)
    {
        displacement = positions.at(N_frames-1).at(objectID); 
        quaternion = quaternions.at(N_frames-1).at(objectID); 
    }
    else 
    {
        // compute frames
        const int lFrame = std::min<int>((int)floor((time - _timeStart) / _timestep), N_frames-1); 
        const int hFrame = std::min<int>(lFrame+1, N_frames-1);

        // interpolate points
        const Point3<REAL> &lPoint = positions.at(lFrame).at(objectID);
        const Point3<REAL> &hPoint = positions.at(hFrame).at(objectID);
        const REAL ratio = (time - timeStamps.at(lFrame)) / _timestep; 
        displacement.x = lPoint.x*(1.0-ratio) + hPoint.x*ratio; 
        displacement.y = lPoint.y*(1.0-ratio) + hPoint.y*ratio; 
        displacement.z = lPoint.z*(1.0-ratio) + hPoint.z*ratio; 

        // interpolate quaternions
        const Quaternion<REAL> &lQuat = quaternions.at(lFrame).at(objectID); 
        const Quaternion<REAL> &hQuat = quaternions.at(hFrame).at(objectID); 
        quaternion = lQuat.slerp(ratio, hQuat); 
    }
}
