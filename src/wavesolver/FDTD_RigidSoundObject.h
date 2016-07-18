#ifndef FDTD_RIGID_SOUND_OBJECT_H 
#define FDTD_RIGID_SOUND_OBJECT_H 
#include <Eigen/Dense>
#include <modal_model/ModalAnalysisObject.h> 
#include <modal_model/ImpulseSeriesObject.h> 
#include <wavesolver/FDTD_RigidObject.h>

//##############################################################################
// Class that manages object that is characterized by surface mesh, level-set
// functions (FDTD_RigidObject) and has impulse (ImpulseSeriesObject) 
// prescribed by the rigidsim tool. 
//##############################################################################
class FDTD_RigidSoundObject : public FDTD_RigidObject, public ImpulseSeriesObject, public ModalAnalysisObject
{
    public: 
        struct ModeAttribute{ REAL modeMax; REAL modeMin; REAL absMax; };

    private: 
        Eigen::VectorXd _activeModeValues; 
        ModeAttribute   _activeModeAttributes; 

    public: 
        // build object
        FDTD_RigidSoundObject()
            : FDTD_RigidObject(), 
              ImpulseSeriesObject(), 
              ModalAnalysisObject()
        {
        }

        // build object with mesh, sdf
        FDTD_RigidSoundObject(const std::string &meshFileName, const int &resolution, const std::string &sdfFilePrefix, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : FDTD_RigidObject(meshFileName, resolution, sdfFilePrefix, meshName, scale), 
              ImpulseSeriesObject(GetMeshPtr()), 
              ModalAnalysisObject()
        {
        }

        // build object with mesh, sdf, modes
        FDTD_RigidSoundObject(const std::string &meshFileName, const int &resolution, const std::string &sdfFilePrefix, const std::string &modeFile, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : FDTD_RigidObject(meshFileName, resolution, sdfFilePrefix, meshName, scale), 
              ImpulseSeriesObject(GetMeshPtr()), 
              ModalAnalysisObject(modeFile)
        {
        }

        void SetVertexModeValues(const int &modeIndex); 
        void GetVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues); 
        void GetVertexModeValuesNormalized(const int &modeIndex, Eigen::VectorXd &modeValues); 
};

#endif
