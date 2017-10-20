#ifndef SHELL_VIBRATION_SOURCE_H
#define SHELL_VIBRATION_SOURCE_H
#include <TYPES.h> 
#include <wavesolver/FDTD_ShellObject.h>
#include <wavesolver/VibrationalSource.h> 
//##############################################################################
// Class ShellVibrationalSource
//##############################################################################
class ShellVibrationalSource : public VibrationalSource
{
    private: 
        FDTD_ShellObject_Ptr _owner; 

    public: 
        ShellVibrationalSource(FDTD_ShellObject_Ptr owner)
            : _owner(owner)
        {}
        virtual REAL Evaluate(const Vector3d &position, 
                              const Vector3d &normal, 
                              const REAL &time, const int &hintTriangle=-1); 
        virtual REAL Evaluate(const int &vertexID, 
                              const Vector3d &vertexNormal, 
                              const REAL &time); 
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time);
        virtual REAL EvaluateVelocity(const Vector3d &position, 
                                      const Vector3d &normal, 
                                      const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, 
                                          const Vector3d &normal, 
                                          const REAL &time); 
};

#endif
