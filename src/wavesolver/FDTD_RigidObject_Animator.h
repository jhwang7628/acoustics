#ifndef FDTD_RIGIDOBJECT_ANIMATOR_H
#define FDTD_RIGIDOBJECT_ANIMATOR_H 
#include <vector> 
#include <string>
#include <memory>
#include <config.h> 
#include <geometry/Point3.hpp> 
#include <linearalgebra/Vector3.hpp> 
#include <linearalgebra/Quaternion.hpp> 
#include <io/RigidObjDispReader.h>

//##############################################################################
// This class reads the rigidsim results and can be queried for interpolated
// kinematics of the objects in the scene. Should be own by managing class such
// as FDTD_AcousticSimulator. 
//
// Linear interpolation will be used.
//##############################################################################
class FDTD_RigidObject_Animator
{
    private: 
        std::shared_ptr<RigidObjDispReader> _rigidsimResults; 
        REAL _timeStart; 
        REAL _timeStop; 
        REAL _timestep; 

    public: 
        FDTD_RigidObject_Animator(){}

        void ReadDisplacement(const std::string &filename); 
        void ReadAllKinematics(const std::string &fileDisplacement, const std::string &fileVelocity, const std::string &fileAcceleration);
        void GetRigidObjectTransform(const int &objectID, const REAL &time, Point3d &displacement, Quaternion<REAL> &quaterion);
};

#endif
