#ifndef FDTD_OBJECTS_H
#define FDTD_OBJECTS_H
#include <unordered_map>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/FDTD_ShellObject.h>
#include <wavesolver/PressureSource.h>
#include <wavesolver/FDTD_RigidObject_Animator.h>
#include <wavesolver/FDTD_PlaneConstraint.h>

//##############################################################################
// Forward declaration
//##############################################################################
class SimWorld;

//##############################################################################
// Struct TriangleIdentifier
//##############################################################################
struct TriangleIdentifier
{
    int objectID;
    int triangleID;
    bool active;
    Vector3d centroid;
    Vector3d normal;
    TriangleIdentifier()
        : objectID(-1), triangleID(-1)
    {}
    TriangleIdentifier(const int &o_id, const int &t_id)
        : objectID(o_id), triangleID(t_id)
    {}
    TriangleIdentifier(const int &o_id, const int &t_id,
                       const Vector3d &c, const Vector3d &n)
        : objectID(o_id), triangleID(t_id), centroid(c), normal(n)
    {}
};
using TriangleIdentifier_UPtr = std::unique_ptr<TriangleIdentifier>;

//##############################################################################
// Struct TIComp
//   Comparator for TriangleIdentifier class
//##############################################################################
struct TIComp
{
    bool operator() (const TriangleIdentifier &lhs,
                     const TriangleIdentifier &rhs)
    {
        // Szudzik's function
        int lhsid, rhsid;
        {
            const unsigned long long _a = std::min(lhs.objectID, lhs.triangleID);
            const unsigned long long _b = std::max(lhs.objectID, lhs.triangleID);
            lhsid = (_a>=_b ? _a*_a+_a+_b : _a+_b*_b);
        }
        {
            const unsigned long long _a = std::min(rhs.objectID, rhs.triangleID);
            const unsigned long long _b = std::max(rhs.objectID, rhs.triangleID);
            rhsid = (_a>=_b ? _a*_a+_a+_b : _a+_b*_b);
        }
        return lhsid < rhsid;
    }
};

//##############################################################################
// Handles a list of objects embedded in the wave solver
//##############################################################################
class FDTD_Objects;
using FDTD_Objects_Ptr = std::shared_ptr<FDTD_Objects>;
class FDTD_Objects
{
    private:
        std::unordered_map<std::string, FDTD_PlaneConstraint_Ptr> _constraints;
        std::unordered_map<int, RigidObjectPtr>                   _rigidObjects;
        std::vector<PressureSourcePtr>                            _pressureSources;
        FDTD_RigidObject_Animator_Ptr                             _objectAnimator;
        bool                                                      _disableAllEvals = false;

    public:
        inline int N() const {return _rigidObjects.size();}
        inline int N_sources() const {return _pressureSources.size();}
        inline int N_constraints() const {return _constraints.size();}
        inline FDTD_RigidObject &Get(const int &ind){return *(_rigidObjects.at(ind));}
        inline RigidObjectPtr GetPtr(const int &ind){return _rigidObjects.at(ind);}
        inline PressureSourcePtr &GetPressureSourcePtr(const int &ind){return _pressureSources.at(ind);}
        inline std::string GetMeshName(const int &ind) const {return _rigidObjects.at(ind)->GetMeshName();}
        //inline int GetMeshID(const string &meshName) const {return _meshIDMap.at(meshName);}
        inline std::vector<PressureSourcePtr> &GetPressureSources(){return _pressureSources;}
        inline const auto &GetRigidObjects() const {return _rigidObjects;}
        inline auto &GetRigidObjects(){return _rigidObjects;}
        inline auto &GetConstraint(const std::string &id){return _constraints;}
        inline auto &GetConstraints(){return _constraints;}
        inline bool HasModalObject() const
        {
            bool has=false;
            for (const auto &m : _rigidObjects)
                if (m.second->Type() == RIGID_SOUND_OBJ)
                    has = (has || std::dynamic_pointer_cast<FDTD_RigidSoundObject>(m.second) ->IsModalObject());
            return has;
        }
        inline bool HasExternalPressureSources(){return _pressureSources.size()>0;}
        // add object if objectName is not in the map
        void AddObject(const int &objectName, RigidObjectPtr &object);
        int AddConstraint(const std::string &id, FDTD_PlaneConstraint_Ptr &constraint);
        int AddConstraints(FDTD_Objects_Ptr &rhs);
        // return index of the object that occupies the position, -1 if none.
        // in the case where multiple objects are occupying that position (due
        // to numerical errors), return the first in the vector
        int OccupyByObject(const Vector3d &positionWorld);
        bool OccupyByConstraint(const Vector3d &pos);
        bool OccupyByConstraint(const Vector3d &pos, FDTD_PlaneConstraint_Ptr &constraint);
        bool TriangleCubeIntersection(const Vector3d &cubeCenter,
                                      const Vector3d &cubeHalfSize,
                                      std::set<TriangleIdentifier, TIComp> &outTris,
                                      const BoundingBox *bbox = nullptr);
        bool TriangleSetCubeIntersection(const Vector3d &maxCubeBound,
                                      const Vector3d &minCubeBound,
                                      const std::set<TriangleIdentifier, TIComp> &candidateTris,
                                      std::set<TriangleIdentifier, TIComp> &outTris);
        REAL ObjectDistance(const int &objectIndex, const Vector3d &positionWorld);
        REAL LowestObjectDistance(const Vector3d &positionWorld);
        bool LowestObjectDistance(const Vector3d &positionWorld, REAL &distance, int &objectID);
        bool LowestThinObjectDistance(const Vector3d &positionWorld, REAL &distance, int &objectID, const std::set<TriangleIdentifier, TIComp> &candidates);
        bool LowestConstraintDistance(const Vector3d &positionWorld, REAL &unsignedDistance, FDTD_PlaneConstraint_Ptr &constraint);
        void ObjectNormal(const int &objectIndex, const Vector3d &positionWorld, Vector3d &queriedNormal);
        void AddVibrationalSourceToObject(const int &objectIndex, VibrationalSourcePtr &sourcePtr);
        void AddPressureSource(PressureSourcePtr &sourcePtr);
        REAL EvaluateNearestVibrationalSources(const Vector3d &position, const Vector3d &normal, const REAL &time);
        REAL EvaluatePressureSources(const Vector3d &position, const Vector3d &normal, const REAL &time);
        bool ReflectAgainstAllBoundaries(const int &startObjectID, const Vector3d &originalPoint, const REAL &time, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &accumulatedBoundaryConditionValue, const REAL &density, const int &maxIteration=1);
        REAL AdvanceAllModalODESolvers(const int &N_steps);
        REAL GetEarliestImpactEvent();
        REAL GetEarliestEventTime(const REAL &startTime) const;
        void SetObjectStates(const REAL time);
        bool UpdateSourceTimes(const REAL time);
        void InitializeAnimator(const std::map<std::string, KinematicsMetadata> &objKinFiles);
        void InitializeAnimator(const std::string &fileDisplacement,
                                const std::string &fileVelocity,
                                const std::string &fileAcceleration);
        void AnimateObjects(const REAL toTime);
        void Join(FDTD_Objects_Ptr &toBeJoined);
        // for time parallelization - disable evaluations if we are in the overlap region
        void DisableAllEvals();
        void EnableAllEvals();
        bool AreEvalsDisabled() { return _disableAllEvals; }
        bool ShaderIsZero(const REAL t) const;

        //// debug methods ////
        std::ofstream *of_q;
        void TestObjectDistanceField(const size_t &ind);
        void PrintAllSources(std::ostream &os);
        void WriteFailedReflections(const std::string &file);
        void ClearFailedReflections();
        void DebugWriteModalQ(const int steps, const std::string &filename);

    friend std::ostream &operator <<(std::ostream &os, const FDTD_Objects &objects);
    friend SimWorld;
};

#endif
