#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include "wavesolver/ModalTransferVibrationalSource.h"

//##############################################################################
//##############################################################################
ModalTransferVibrationalSource::
ModalTransferVibrationalSource(RigidObjectPtr owner)
    : VibrationalSource(owner)
{
    _owner = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(owner);
}

//##############################################################################
//##############################################################################
ModalTransferVibrationalSource::
ModalTransferVibrationalSource(RigidObjectPtr owner,
                               const std::vector<int> &modes)
    : VibrationalSource(owner),
      _useAllModes(false),
      _modes(modes)
{
    _owner = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(owner);
    for (const int &m : modes)
        if (m<0)
            _useAllModes = true;
}

//##############################################################################
//##############################################################################
REAL ModalTransferVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal,
         const REAL &time, const int &hintTriangle)
{
    throw std::runtime_error("**ERROR** not implemented");
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL ModalTransferVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &normal,
         const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
    return 0.0;
}

//##############################################################################
//##############################################################################
Vector3d ModalTransferVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    Vector3d r;
    if (_useAllModes)
        for (int m=0; m<_owner->N_Modes(); ++m)
            r += _owner->PerfectHarmonics_SampleModalAcceleration(m, vertexID, time);
    else
        for (const int m : _modes)
            r += _owner->PerfectHarmonics_SampleModalAcceleration(m, vertexID, time);
    return r;
}
