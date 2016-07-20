#ifndef FDTD_RIGID_SOUND_OBJECT_H 
#define FDTD_RIGID_SOUND_OBJECT_H 
#include <iostream>
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

        // modal displacement
        Eigen::VectorXd _qOld; 
        Eigen::VectorXd _qNew; 

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

        inline void InitializeModeVectors(){_qOld.setZero(N_Modes()); _qNew.setZero(N_Modes());}
        void SetVertexModeValues(const int &modeIndex); 
        void GetVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues); 
        void GetVertexModeValuesNormalized(const int &modeIndex, Eigen::VectorXd &modeValues); 
        void GetForceInModalSpace(const ImpactRecord &record, Eigen::VectorXd &forceInModalSpace); 
        void GetModalDisplacementAux(const int &mode, Eigen::VectorXd &displacement);
        void GetModalDisplacement(const int &mode, Eigen::VectorXd &displacement);  // only transform the mode quried
        void GetModalDisplacement(Eigen::VectorXd &displacement); // transform all the mode displacements
        void AdvanceModalODESolvers(const int &N_steps);
        void AdvanceModalODESolvers(const int &N_steps, std::ofstream &os);
};

#endif
