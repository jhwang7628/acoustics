#ifndef FDTD_OBJECTS_H 
#define FDTD_OBJECTS_H 

#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/FDTD_RigidObject.h> 

//##############################################################################
// Handles a list of objects embedded in the wave solver
//##############################################################################
class FDTD_Objects
{
    private: 
        std::vector<RigidObjectPtr> _rigidObjects; 
        std::map<std::string, int>  _meshIDMap; 

    public: 
        inline int N() const {return _rigidObjects.size();} 
        inline FDTD_RigidObject& Get(const int &ind){return *(_rigidObjects[ind]);} 
        inline RigidObjectPtr GetPtr(const int &ind){return _rigidObjects[ind];} 
        inline std::string GetMeshName(const int &ind) const {return std::string(_rigidObjects.at(ind)->GetMeshName());}
        inline int GetMeshID(const string &meshName) const {return _meshIDMap.at(meshName);}
        // add object if objectName is not in the map
        void AddObject(const std::string &objectName, RigidObjectPtr &object); 
        // return index of the object that occupies the position, -1 if none. 
        // in the case where multiple objects are occupying that position (due
        // to numerical errors), return the first in the vector
        int OccupyByObject(const Vector3d &positionWorld); 
        REAL ObjectDistance(const int &objectIndex, const Vector3d &positionWorld); 
        void ObjectNormal(const int &objectIndex, const Vector3d &positionWorld, Vector3d &queriedNormal); 
        void AddVibrationalSourceToObject(const int &objectIndex, VibrationalSourcePtr &sourcePtr);

        //// debug methods ////
        void TestObjectDistanceField(const size_t &ind); 

    friend std::ostream &operator <<(std::ostream &os, const FDTD_Objects &objects);
};

#endif 
