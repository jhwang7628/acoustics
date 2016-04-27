#ifndef VIBRATIONAL_SOURCE_H 
#define VIBRATIONAL_SOURCE_H 

#include <TYPES.h> 
#include <linearalgebra/Vector3.hpp> 
#include <wavesolver/Source.h> 
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 

//##############################################################################
// Represents the surface acceleration that causes sound. 
//##############################################################################
class FDTD_RigidObject; 
class VibrationalSource : public Source
{
    protected:
        std::shared_ptr<FDTD_RigidObject> _owner; 

    public:
        VibrationalSource(){}
        VibrationalSource(RigidObjectPtr owner)
            : _owner(owner)
        {}

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)=0; 
};

#endif
