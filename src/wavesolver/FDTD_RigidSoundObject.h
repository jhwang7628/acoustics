#ifndef FDTD_RIGID_SOUND_OBJECT_H 
#define FDTD_RIGID_SOUND_OBJECT_H 
#include <iostream>
#include <Eigen/Dense>
#include <modal_model/ModalAnalysisObject.h> 
#include <modal_model/ImpulseSeriesObject.h> 
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/AccelerationNoiseVibrationalSource.h>
#include <utils/SimpleTimer.h>

//##############################################################################
// Class that manages object that is characterized by surface mesh, level-set
// functions (FDTD_RigidObject) and has impulse (ImpulseSeriesObject) 
// prescribed by the rigidsim tool. 
//
// Note: be careful when fetching velocity and acceleration data, see comments
// below.
//
// The stored fields follows the timeline
// tn-1         tn         tn+1        tn+2
//  n-1   n-1    n     n    n+1   n+1   n+2
// --p-----v-----p-----v-----p-----v-----p--
//  q_p         q_c         q_n         q_nn
//    qDot_c_minus
//            qDDot_c
//               qDDot_c_plus
//
// When AdvanceModalODESolvers() is called, it will compute q_nn and update all
// derivatives needed. Therefore it is expected the impact forces are fetched
// in the interval [tn+1, tn+2] is used. This should be called before acoustic
// time stepping so that derivatives are updated properly.
//
// After the acoustic time stepping, call UpdateQPointers() so that the
// internal timestamps for ODE solvers are updated and pointers for q arrays
// are shifted properly. 
//
// For safety, I limit the timestamp access when sampling velocity and
// acceleration to make sure its called correctly. See individual functions. 
//##############################################################################
class FDTD_RigidSoundObject : public FDTD_RigidObject, public ModalAnalysisObject, public ImpulseSeriesObject
{
    public: 
        struct ModeAttribute{ REAL modeMax; REAL modeMin; REAL absMax; };

    private: 
        Eigen::VectorXd _activeModeValues; 
        ModeAttribute   _activeModeAttributes; 

        Eigen::VectorXd _q_p;  // previous displacement
        Eigen::VectorXd _q_c;  // current displacement
        Eigen::VectorXd _q_n;  // next displacement
        Eigen::VectorXd _q_nn; // next next displacement
        Eigen::VectorXd _qDot_c_minus; // previous half step velocity
        Eigen::VectorXd _qDDot_c; // current acceleration
        Eigen::VectorXd _qDDot_c_plus; // next half step acceleration

        // for vectorized IIR
        Eigen::VectorXd _coeff_qNew;  // 2 epsilon cos(theta)
        Eigen::VectorXd _coeff_qOld;  // -epsilon^2
        Eigen::VectorXd _coeff_Q;     // 2(epsilon cos(theta+gamma) - epsilon^2 cos(2theta + gamma))/(3 omega omega_d m)

        // Timers
        SimpleTimer _timer_mainstep[3];  // advance, q->u, IO
        SimpleTimer _timer_substep_advanceODE[3]; // interpolate force, transform force, step system
        SimpleTimer _timer_substep_q2u[1]; 

        // cached properties
        bool _invInertiaTensorCached=false; 
        Matrix3<REAL> _invInertiaTensor; 

        bool _animated = false; // if an object is being animated by rbd sim, then this should be set to true

    public: 
        // build object
        FDTD_RigidSoundObject()
            : FDTD_RigidObject(), 
              ModalAnalysisObject(),
              ImpulseSeriesObject()
        {
        }

        // build object with mesh, sdf
        FDTD_RigidSoundObject(const std::string &workingDirecotry, const int &resolution, const std::string &objectPrefix, const bool &buildFromTetMesh, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : FDTD_RigidObject(workingDirecotry, resolution, objectPrefix, buildFromTetMesh, meshName, scale), 
              ModalAnalysisObject(),
              ImpulseSeriesObject(GetMeshPtr())
        {
        }

        // build object with mesh, sdf, modes
        FDTD_RigidSoundObject(const std::string &workingDirecotry, const int &resolution, const std::string &objectPrefix, const std::string &modeFile, const bool &buildFromTetMesh, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : FDTD_RigidObject(workingDirecotry, resolution, objectPrefix, buildFromTetMesh, meshName, scale), 
              ModalAnalysisObject(modeFile),
              ImpulseSeriesObject(GetMeshPtr())
        {
        }

        inline bool Animated(){return _animated;}
        inline bool IsModalObject(){return N_Modes()>0;}
        inline void SetAnimated(const bool &is){_animated = is;}
        inline Vector3d PremultiplyInvInertiaTensor(const Vector3d &x){return _invInertiaTensor * x;}

        void Initialize(); 
        void InitializeModeVectors(); 
        void GetVertexModeValuesNormalized(const int &modeIndex, Eigen::VectorXd &modeValues); 
        void GetForceInModalSpace(const ImpactRecord &record, Eigen::VectorXd &forceInModalSpace); 
        void GetModalDisplacementAux(const int &mode, Eigen::VectorXd &displacement);
        void GetModalDisplacement(const int &mode, Eigen::VectorXd &displacement);  // only transform the mode quried
        void GetModalDisplacement(Eigen::VectorXd &displacement); // transform all the mode displacements
        void AdvanceModalODESolvers(const int &N_steps);
        void AdvanceModalODESolvers(const int &N_steps, const int &mode, std::ofstream &of_displacement, std::ofstream &of_q);
        void UpdateQPointers(); 
        const Point3d &CenterOfMass() const {return _volumeCenter;} // volume center = mass center
        REAL Mass() const; 
        void InvInertiaTensor(Matrix3<REAL> &I_inv, const bool &compute=false); 
        REAL EffectiveMass(const Vector3d &x, const Vector3d &n); 
        // Since velocity and acceleration are estimated using central difference, their values correspond to qOld 
        // and thus when fetching, we should be getting solution values at time t=_time - _ODEStepSize. 
        REAL SampleModalDisplacement(const Vector3d &samplePoint, const Vector3d &samplePointNormal, const REAL &sampleTime); 
        REAL SampleModalVelocity(const Vector3d &samplePoint, const Vector3d &samplePointNormal, const REAL &sampleTime); 
        REAL SampleModalAcceleration(const Vector3d &samplePoint, const Vector3d &samplePointNormal, const REAL &sampleTime); 
        REAL SampleModalAcceleration(const int &vertexID, const Vector3d &vertexNormal, const REAL &sampleTime); 

        REAL EstimateContactTimeScale(const int &vertex_a, const REAL &contactSpeed, const Vector3d &impulse_a); 
        REAL EstimateContactTimeScale(const std::shared_ptr<FDTD_RigidSoundObject> &object_b, const int &vertex_a, const int &vertex_b, const REAL &contactSpeed, const Vector3d &impulse_a); 

        ///// debug methods /////

        // Perfect harmonics: q(t) = cos(omega * t), omega is angular frequency
        // of mode
        REAL PerfectHarmonics_SampleModalVelocity(const int &mode, const Vector3d &samplePoint, const Vector3d &samplePointNormal, const REAL &sampleTime); 
        REAL PerfectHarmonics_SampleModalAcceleration(const int &mode, const Vector3d &samplePoint, const Vector3d &samplePointNormal, const REAL &sampleTime); 
        void PrintAllVelocity(const std::string &filename, const int &mode) const;

    friend class AccelerationNoiseVibrationalSource; 
};

#endif
