#ifndef FDTD_SHELLOBJECT_H
#define FDTD_SHELLOBJECT_H
#include <queue>
#include "config.h"
#include "boost/filesystem.hpp"
#include "geometry/BoundingBox.h"
#include "wavesolver/FDTD_RigidObject.h"
//##############################################################################
// Class FDTD_ShellObject
//##############################################################################
class FDTD_ShellObject : public FDTD_RigidObject
{
private:
    struct _DataStep
    {
        static const REAL stepSize; 
        int               frameID; 
        REAL              frameTime; 
        std::vector<REAL> acceleration; // length: N_total_vertices*3 
        std::vector<REAL> displacement;
        std::vector<REAL> constraint_displacement; 
        std::vector<REAL> constraint_acceleration; 
        std::vector<int>  constrainedVertices; 
    };
    using _DataStep_UPtr = std::unique_ptr<_DataStep>; 
    struct _DataBuffer
    {
        static std::queue<boost::filesystem::path> qPrefixes; // to-be-read file prefixes
        const std::string                          accSuffix = "wsacc"; 
        const std::string                          disSuffix = "displacement"; 
        const std::string                          cdisSuffix = "constraint_displacement"; 
        const std::string                          caccSuffix = "constraint_acceleration"; 
        const int                                  bufLen = 10; 
        unsigned int                               bufValidLen;
        REAL                                       bufStartTime = -1.0;
        FDTD_ShellObject                          *owner; 
        std::vector<_DataStep>                     buf; 
        REAL                                       lastQueriedTime; 
        _DataBuffer()
            : owner(nullptr)
        {} 
        _DataBuffer(FDTD_ShellObject *o)
            : owner(o)
        {
            buf.resize(bufLen);
            bufValidLen = 0;
        }
        bool Empty() const {return bufValidLen != 0;}
        const _DataStep &Data(const int i) const {return buf.at(i);}
        bool FindData(const REAL time, _DataStep **data); 
        void ReadNextBuffer(); 
        void ReadMetaData(); 
    };

    bool                   _init = false; 
    std::vector<int>       _o2iMap; // vertex map original->internal
    std::vector<int>       _i2oMap; // vertex map internal->original
    std::vector<int>       _constrainedVertices; 
    std::string            _dataDir;
    std::string            _vertexMapFile = "vertex_map.txt";
    mutable _DataStep     *_currentData = nullptr; // allocated on stack, mere pointer (don't delete)
    mutable _DataBuffer    _dataBuffer; 
public:
    FDTD_ShellObject()
        : FDTD_RigidObject(SHELL_OBJ),
          _dataBuffer(this)
    {}
    FDTD_ShellObject(const std::string &workingDirectory,
                     const std::string &shellDataDirectory,
                     const int &resolution, 
                     const std::string &objectPrefix,
                     const std::shared_ptr<PML_WaveSolver_Settings> &sset,
                     const std::string &meshName)
        : FDTD_RigidObject(SHELL_OBJ,
                           workingDirectory,
                           resolution, 
                           objectPrefix,
                           false,
                           sset,
                           meshName),
          _dataDir(shellDataDirectory), 
          _dataBuffer(this)
    {}

    virtual REAL DistanceToMesh(const double &x, const double &y, const double &z)
    {const Vector3d pos(x,y,z); return DistanceToMesh(pos);}
    virtual REAL DistanceToMesh(const Vector3d &position); 
    virtual bool NormalToMesh(const double &x, const double &y, const double &z, 
                              Vector3d &queriedNormal)
    {return false;}
    virtual bool NormalToMesh(const Vector3d &position, Vector3d &queriedNormal)
    {return NormalToMesh(position.x, position.y, position.z, queriedNormal);}
    virtual int ReflectAgainstBoundary(const Vector3d &originalPoint, 
                                       Vector3d &reflectedPoint, 
                                       Vector3d &boundaryPoint, 
                                       Vector3d &erectedNormal, 
                                       REAL &distanceTravelled, 
                                       const int &startFromTriangle=-1);
    virtual void UpdateBoundingBox(); 

    void Initialize(); 
    void UpdatePosAcc(const REAL time);
    void GetAllVertexPos(std::vector<Point3d> &allp)const; 
    bool GetDisplacements(std::vector<REAL> &alld)const; 
    Vector3d GetVertexPos(const int v_idx) const; 
    Vector3d GetVertexAcc(const int v_idx) const; 
    int ReflectAgainstBoundary(const Vector3d &originalPoint, 
                               const std::set<int> &triangles,
                               Vector3d &reflectedPoint, 
                               Vector3d &boundaryPoint, 
                               Vector3d &erectedNormal, 
                               REAL &distance);
};
using FDTD_ShellObject_Ptr = std::shared_ptr<FDTD_ShellObject>; 
#endif
