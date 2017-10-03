#ifndef FDTD_SHELLOBJECT_H
#define FDTD_SHELLOBJECT_H
#include "config.h"
#include "geometry/BoundingBox.h"
#include "wavesolver/FDTD_RigidObject.h"
//##############################################################################
// Class FDTD_ShellObject
//##############################################################################
class FDTD_ShellObject : public FDTD_RigidObject
{
private:
    bool _init = false; 
    BoundingBox _bbox; 
public:
    FDTD_ShellObject()
        : FDTD_RigidObject(SHELL_OBJ)
    {}
    FDTD_ShellObject(const std::string &workingDirectory,
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
                           meshName)
    {}

    virtual REAL DistanceToMesh(const double &x, const double &y, const double &z)
    {const Vector3d pos(x,y,z); return DistanceToMesh(pos);}
    virtual REAL DistanceToMesh(const Vector3d &position); 
    virtual bool NormalToMesh(const double &x, const double &y, const double &z, 
                              Vector3d &queriedNormal)
    {return false;}
    virtual bool NormalToMesh(const Vector3d &position, Vector3d &queriedNormal)
    {return NormalToMesh(position.x, position.y, position.z, queriedNormal);}

    void Initialize(); 
};
#endif
