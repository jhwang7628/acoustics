#ifndef FDTD_OBJECTS_H 
#define FDTD_OBJECTS_H 

#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/PressureSource.h>

//##############################################################################
// Handles a list of objects embedded in the wave solver
//##############################################################################
class FDTD_Objects
{
    private: 
        std::vector<RigidSoundObjectPtr>    _rigidObjects; 
        std::vector<PressureSourcePtr>      _pressureSources; 
        std::map<std::string, int>          _meshIDMap; 

    public: 
        inline int N() const {return _rigidObjects.size();} 
        inline int N_sources() const {return _pressureSources.size();}
        inline FDTD_RigidSoundObject &Get(const int &ind){return *(_rigidObjects[ind]);} 
        inline RigidSoundObjectPtr GetPtr(const int &ind){return _rigidObjects.at(ind);} 
        inline PressureSourcePtr &GetPressureSourcePtr(const int &ind){return _pressureSources.at(ind);} 
        inline std::string GetMeshName(const int &ind) const {return std::string(_rigidObjects.at(ind)->GetMeshName());}
        inline int GetMeshID(const string &meshName) const {return _meshIDMap.at(meshName);}
        inline std::vector<PressureSourcePtr> &GetPressureSources(){return _pressureSources;}
        inline const std::vector<RigidSoundObjectPtr> &GetRigidSoundObjects() const {return _rigidObjects;}
        inline std::vector<RigidSoundObjectPtr> &GetRigidSoundObjects(){return _rigidObjects;}
        inline bool HasExternalPressureSources(){return _pressureSources.size()>0;}
        // add object if objectName is not in the map
        void AddObject(const std::string &objectName, RigidSoundObjectPtr &object); 
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
        bool ReflectAgainstAllBoundaries(const Vector3d &originalPoint, const REAL &time, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &accumulatedBoundaryConditionValue, const REAL &density, const int &maxIteration=1); 
        REAL AdvanceAllModalODESolvers(const int &N_steps); 
        REAL GetEarliestImpactEvent();

        //// debug methods ////
        void TestObjectDistanceField(const size_t &ind); 
        void PrintAllSources(std::ostream &os); 

    friend std::ostream &operator <<(std::ostream &os, const FDTD_Objects &objects);
};

#endif 
