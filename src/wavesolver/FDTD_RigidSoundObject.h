#ifndef FDTD_RIGID_SOUND_OBJECT_H 
#define FDTD_RIGID_SOUND_OBJECT_H 
#include <modal_model/ModalAnalysis.h> 
#include <modal_model/ImpulseSeriesObject.h> 
#include <wavesolver/FDTD_RigidObject.h>

//##############################################################################
// Class that manages object that is characterized by surface mesh, level-set
// functions (FDTD_RigidObject) and has impulse (ImpulseSeriesObject) 
// prescribed by the rigidsim tool
//##############################################################################
class FDTD_RigidSoundObject : public FDTD_RigidObject, public ImpulseSeriesObject
{
    public: 
        FDTD_RigidSoundObject()
            : FDTD_RigidObject(), ImpulseSeriesObject()
        {
        }
        FDTD_RigidSoundObject(const std::string &meshFileName, const int &resolution, const std::string &sdfFilePrefix, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : FDTD_RigidObject(meshFileName, resolution, sdfFilePrefix, meshName, scale), ImpulseSeriesObject(GetMeshPtr())
        {
        }
};

#endif
