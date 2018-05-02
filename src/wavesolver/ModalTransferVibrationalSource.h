#ifndef MODAL_TRANSFER_VIBRATIONAL_SOURCE
#define MODAL_TRANSFER_VIBRATIONAL_SOURCE
#include "wavesolver/VibrationalSource.h"
//##############################################################################
// Class ModalTransferVibrationalSource
//##############################################################################
class ModalTransferVibrationalSource : public VibrationalSource
{
    private:
        RigidSoundObjectPtr _owner;
        bool _useAllModes = true;
        std::vector<int> _modes;

    public:
        ModalTransferVibrationalSource(RigidObjectPtr owner);
        ModalTransferVibrationalSource(RigidObjectPtr owner,
                                       const std::vector<int> &modes);
        inline void SetOwner(const RigidSoundObjectPtr &ptr){_owner = ptr;}
        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal,
                              const REAL &time, const int &hintTriangle=-1);
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal,
                              const REAL &time);
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time);
};
#endif
