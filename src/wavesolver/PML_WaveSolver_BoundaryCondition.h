#ifndef PML_WAVESOLVER_BOUNDARY_CONDITION_H 
#define PML_WAVESOLVER_BOUNDARY_CONDITION_H 

#include <TYPES.h> 
#include <linearalgebra/Vector3d.hpp> 
//##############################################################################
// Represents the surface acceleration that causes sound. 
//##############################################################################
class VibrationalSource
{
    private:
        std::shared_ptr<FDTD_RigidObject> _owner; 

    public:
        VibrationalSource(){}
        VibrationalSource(std::shared_ptr<FDTD_RigidObject> owner)
            : _owner(owner)
        {}

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)=0; 
};

#endif
