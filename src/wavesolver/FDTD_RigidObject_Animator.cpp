#include <fstream>
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
// Function GetRigidObjectTransform
//##############################################################################
void FDTD_RigidObject_Animator::
GetRigidObjectTransform(const int &objectID, const REAL &time, 
                        Point3d &newCOM, Quaternion<REAL> &quaternion)
{
    const int N_frames = _rigidsimResults->N_Frames(); 
    const auto &timeStamps = _rigidsimResults->Timesteps(); 
    const std::vector<std::vector<Point3<REAL> > > &positions  = _rigidsimResults->Positions();
    const std::vector<std::vector<Quaternion<REAL> > > &quaternions  = _rigidsimResults->Rotations();
    if (time <= _timeStart)
    {
        newCOM = positions.at(0).at(objectID); 
        quaternion = quaternions.at(0).at(objectID); 
    }
    else if (time >= _timeStop)
    {
        newCOM = positions.at(N_frames-1).at(objectID); 
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
        newCOM.x = lPoint.x*(1.0-ratio) + hPoint.x*ratio; 
        newCOM.y = lPoint.y*(1.0-ratio) + hPoint.y*ratio; 
        newCOM.z = lPoint.z*(1.0-ratio) + hPoint.z*ratio; 

        // interpolate quaternions
        const Quaternion<REAL> &lQuat = quaternions.at(lFrame).at(objectID); 
        const Quaternion<REAL> &hQuat = quaternions.at(hFrame).at(objectID); 
        quaternion = lQuat.slerp(ratio, hQuat); 
    }
}

//##############################################################################
// Function ReadObjectKinFile
//   file format: 
//     frame_id translation_x translation_y translation_z quaternion_w  ...
//     quaternion_y quaternion_z 
//##############################################################################
void FDTD_RigidObject_Blender_Animator::
ReadObjectKinData(const std::string &objId, 
                  const KinematicsMetadata &meta,
                  const bool blenderRotateYUp)
{
    const std::string kinFile = meta.fileKinematics; 
    std::cout << "Reading object kinematic file: " << kinFile << std::endl; 
    std::ifstream stream(kinFile.c_str()); 
    if (!stream) 
        throw std::runtime_error("**ERROR** Cannot read object kinematics file: " + kinFile); 

    auto &data = _objectsData[objId].data; 
    _objectsData.at(objId).meta = meta; 
    std::string line; 
    double buf_w;
    Vector3d buf_v; 
    while(std::getline(stream, line))
    {
        StepData sd;
        std::istringstream iss(line); 
        iss >> sd.frame; 
        iss >> sd.translation.x;
        iss >> sd.translation.y;
        iss >> sd.translation.z;
        iss >> buf_w; 
        iss >> buf_v[0]; 
        iss >> buf_v[1]; 
        iss >> buf_v[2]; 
        // to prevent numerical error, clamp w
        buf_w = std::max<REAL>(-1.0, std::min<REAL>(1.0, buf_w));  
        sd.rotation = Quaternion<REAL>(buf_w, buf_v); 
        // blender uses z-up as default so we apply -90 rot about x-axis
        if (blenderRotateYUp)
        {
            Quaternion<REAL> corr = 
                Quaternion<REAL>::fromAxisRotD(Vector3<REAL>(1.0, 0.0, 0.0),
                                               -90.);
            sd.rotation = sd.rotation*corr; 
        }
        std::cout << "read step data: "  // FIXME debug
                  << sd.frame << " " 
                  << sd.translation << " " 
                  << sd.rotation << std::endl; 
        data.push_back(sd); 
    }
}

//##############################################################################
// Function GetRigidObjectTransform
//##############################################################################
void FDTD_RigidObject_Blender_Animator::
GetRigidObjectTransform(const int &objectID, const REAL &time, 
                        Point3d &newCOM, Quaternion<REAL> &quaternion)
{
    const std::string objName = std::to_string(objectID); 
    const auto &it = _objectsData.find(objName); 
    const auto &objData = it->second; 
    if (it == _objectsData.end() || objData.data.size() == 0)
        return; 

    const int idx_0 = std::min<int>(
            std::max<int>(0, (time-_timeStart)/objData.meta.stepSize), 
            objData.data.size()-1); 
    const int idx_1 = std::min<int>(idx_0+1, objData.data.size()-1); 
    const REAL t_0 = idx_0*objData.meta.stepSize; 
    const REAL t_1 = idx_1*objData.meta.stepSize; 

    const REAL alpha = (idx_0 == idx_1 ? 1.0 : (time - t_0)/(t_1 - t_0)); 
    newCOM = objData.data.at(idx_0).translation*(1.0-alpha)
           + objData.data.at(idx_1).translation*(    alpha); 
    quaternion = objData.data.at(idx_0).rotation.slerp(alpha, 
                 objData.data.at(idx_1).rotation); 
}
