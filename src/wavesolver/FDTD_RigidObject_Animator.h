#ifndef FDTD_RIGIDOBJECT_ANIMATOR_H
#define FDTD_RIGIDOBJECT_ANIMATOR_H 
#include <map>
#include <vector> 
#include <string>
#include <memory>
#include <config.h> 
#include <geometry/Point3.hpp> 
#include <linearalgebra/Vector3.hpp> 
#include <linearalgebra/Quaternion.hpp> 
#include <io/RigidObjDispReader.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

//##############################################################################
// Enum AnimatorType
//##############################################################################
enum AnimatorType
{
    RBD_SOLVER_ANIMATOR = 0, 
    BLENDER_ANIMATOR
}; 

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

    protected: 
        AnimatorType _type; 
        REAL _timeStart; 
        REAL _timeStop; 
        REAL _timestep; 

    public: 
        FDTD_RigidObject_Animator(const AnimatorType &type = RBD_SOLVER_ANIMATOR) 
            : _type(type)
        {}

        void ReadDisplacement(const std::string &filename); 
        void ReadAllKinematics(const std::string &fileDisplacement, const std::string &fileVelocity, const std::string &fileAcceleration);
        AnimatorType Type() const {return _type;} 
        virtual void GetRigidObjectTransform(const int &objectID, const REAL &time, Point3d &displacement, Quaternion<REAL> &quaterion);
};
using FDTD_RigidObject_Animator_Ptr = std::shared_ptr<FDTD_RigidObject_Animator>; 

//##############################################################################
// Class FDTD_RigidObject_Blender_Animator
//   Manages rigid transformation produced by Blender
//##############################################################################
class FDTD_RigidObject_Blender_Animator : public FDTD_RigidObject_Animator
{
    private: 
        struct StepData
        {
            int frame; 
            Vector3<REAL> translation;
            Quaternion<REAL> rotation; 
        };
        struct ObjectData
        {
            KinematicsMetadata meta; 
            std::vector<StepData> data;
        };
        std::map<std::string, ObjectData> _objectsData; 

    public: 
        FDTD_RigidObject_Blender_Animator()
            : FDTD_RigidObject_Animator(BLENDER_ANIMATOR)
        {}

        void ReadObjectKinData(const std::string &objId, 
                               const KinematicsMetadata &meta,
                               const bool blenderRotateYUp); 

        virtual void GetRigidObjectTransform(const int &objectID, 
                const REAL &time, 
                Point3d &displacement, 
                Quaternion<REAL> &quaterion);
}; 
using FDTD_RigidObject_Blender_Animator_Ptr = 
    std::shared_ptr<FDTD_RigidObject_Blender_Animator>; 

#endif
