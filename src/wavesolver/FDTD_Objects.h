#ifndef FDTD_OBJECTS_H 
#define FDTD_OBJECTS_H 
#include <unordered_map> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/PressureSource.h>
#include <wavesolver/FDTD_RigidObject_Animator.h>

//##############################################################################
// Forward declaration
//##############################################################################
class SimWorld;

//##############################################################################
// Handles a list of objects embedded in the wave solver
//##############################################################################
class FDTD_Objects;
using FDTD_Objects_Ptr = std::shared_ptr<FDTD_Objects>; 
class FDTD_Objects
{
    private: 
        std::unordered_map<int, RigidSoundObjectPtr> _rigidObjects; 
        std::vector<PressureSourcePtr>               _pressureSources; 
        FDTD_RigidObject_Animator_Ptr                _objectAnimator; 

    public: 
        inline int N() const {return _rigidObjects.size();} 
        inline int N_sources() const {return _pressureSources.size();}
        inline FDTD_RigidSoundObject &Get(const int &ind){return *(_rigidObjects.at(ind));} 
        inline RigidSoundObjectPtr GetPtr(const int &ind){return _rigidObjects.at(ind);} 
        inline PressureSourcePtr &GetPressureSourcePtr(const int &ind){return _pressureSources.at(ind);} 
        inline std::string GetMeshName(const int &ind) const {return _rigidObjects.at(ind)->GetMeshName();}
        //inline int GetMeshID(const string &meshName) const {return _meshIDMap.at(meshName);}
        inline std::vector<PressureSourcePtr> &GetPressureSources(){return _pressureSources;}
        inline const auto &GetRigidSoundObjects() const {return _rigidObjects;}
        inline auto &GetRigidSoundObjects(){return _rigidObjects;}
        inline bool HasModalObject() const 
        {
            bool has=false; 
            for (const auto &m : _rigidObjects) 
                has = (has || m.second->IsModalObject()); 
            return has;
        }
        inline bool HasExternalPressureSources(){return _pressureSources.size()>0;}
        // add object if objectName is not in the map
        void AddObject(const int &objectName, RigidSoundObjectPtr &object); 
        // return index of the object that occupies the position, -1 if none. 
        // in the case where multiple objects are occupying that position (due
        // to numerical errors), return the first in the vector
        int OccupyByObject(const Vector3d &positionWorld); 
        REAL ObjectDistance(const int &objectIndex, const Vector3d &positionWorld); 
        REAL LowestObjectDistance(const Vector3d &positionWorld); 
        void LowestObjectDistance(const Vector3d &positionWorld, REAL &distance, int &objectID); 
        void ObjectNormal(const int &objectIndex, const Vector3d &positionWorld, Vector3d &queriedNormal); 
        void AddVibrationalSourceToObject(const int &objectIndex, VibrationalSourcePtr &sourcePtr);
        void AddPressureSource(PressureSourcePtr &sourcePtr); 
        REAL EvaluateNearestVibrationalSources(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        REAL EvaluatePressureSources(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        bool ReflectAgainstAllBoundaries(const int &startObjectID, const Vector3d &originalPoint, const REAL &time, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &accumulatedBoundaryConditionValue, const REAL &density, const int &maxIteration=1); 
        REAL AdvanceAllModalODESolvers(const int &N_steps); 
        REAL GetEarliestImpactEvent();
        void SetObjectStates(const REAL time); 
        void InitializeAnimator(const std::string &fileDisplacement,
                                const std::string &fileVelocity,
                                const std::string &fileAcceleration); 
        void AnimateObjects(const REAL toTime); 
        void Join(FDTD_Objects_Ptr &toBeJoined); 

        //// debug methods ////
        void TestObjectDistanceField(const size_t &ind); 
        void PrintAllSources(std::ostream &os); 
        void WriteFailedReflections(const std::string &file);
        void ClearFailedReflections();

    friend std::ostream &operator <<(std::ostream &os, const FDTD_Objects &objects);
    friend SimWorld; 
};

#endif 
