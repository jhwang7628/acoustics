#ifndef MODAL_VIBRATIONAL_SOURCE_H 
#define MODAL_VIBRATIONAL_SOURCE_H 

#include <TYPES.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/FDTD_RigidSoundObject.h> 

//##############################################################################
// This class handles source evaluation for modal objects.
//##############################################################################
class ModalVibrationalSource : public VibrationalSource
{
    private: 
        RigidSoundObjectPtr _modalObjectOwner; 

    public:
        ModalVibrationalSource(RigidObjectPtr owner);

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time); 
};

#endif
