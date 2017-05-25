#ifndef FDTD_ACOUSTIC_SIMULATOR_H
#define FDTD_ACOUTSIC_SIMULATOR_H 

#include <string> 
#include <parser/ImpulseResponseParser.h> 
#include <geometry/BoundingBox.h> 
#include <wavesolver/PML_WaveSolver.h> 
#include <wavesolver/PML_WaveSolver_Settings.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/HarmonicVibrationalSource.h>
#include <wavesolver/ModalVibrationalSource.h>
#include <wavesolver/AccelerationNoiseVibrationalSource.h>
#include <wavesolver/FDTD_Objects.h> 
#include <wavesolver/FDTD_RigidObject_Animator.h>

//#############################################################################
// Sturct SimBox
//#############################################################################
struct SimBox
{
    Tuple3i  rasterizedCenter; 
    Vector3d continuousCenter;
};

//##############################################################################
// A managing class to coordinate wavesolver, rigid body simulator, parser, 
// GUI, IO
//
// Reference for the wavesolver: 
//  [1] Chadwick 2012, Precomputed Acceleration Noise for Improved Rigid-Body Sound
//  [2] Liu 1997, The perfectly matched layer for acoustic waves in absorptive media
//##############################################################################
class FDTD_AcousticSimulator
{
    protected: 
        // meta objects
        std::shared_ptr<ImpulseResponseParser>      _parser; 
        std::shared_ptr<PML_WaveSolver>             _acousticSolver; 
        std::shared_ptr<PML_WaveSolver_Settings>    _acousticSolverSettings; 
        std::shared_ptr<FDTD_Objects>               _sceneObjects; 
        std::shared_ptr<FDTD_RigidObject_Animator>  _sceneObjectsAnimator; 
        SimBox                                      _simBox;

        // state representation
        bool                    _canInitializeSolver; 

        // necessary fields 
        int                     _stepIndex;
        int                     _snapshotIndex;
        double                  _simulationTime; 
        std::string             _configFile; 
        std::string            *_simulatorID = nullptr;  // unique id that is optional

    private: 
        void _ParseSolverSettings();
        void _SetBoundaryConditions(); 
        void _SetPressureSources(); 
        void _SetListeningPoints(); 
        std::string _CompositeFilename(const std::string filename); // composite output path and prefix
        void _SaveSolverSettings(const std::string &filename); 
        void _SavePressureCellPositions(const std::string &filename); 
        void _SaveVelocityCellPositions(const std::string &filename, const int &dim); 
        void _SaveListeningPositions(const std::string &filename); 
        void _SavePressureTimestep(const std::string &filename); 
        void _SaveVelocityTimestep(const std::string &filename, const int &dim); 
        void _SaveListeningData(const std::string &filename);
        void _SaveModalFrequencies(const std::string &filename); 
        void _SaveSimulationSnapshot(const std::string &numbering);

    public: 
        FDTD_AcousticSimulator(std::string *simulatorID=nullptr)
            : FDTD_AcousticSimulator("", simulatorID)
        {}
        FDTD_AcousticSimulator(const std::string &configFile, std::string *simulatorID=nullptr)
            : _stepIndex(0), _snapshotIndex(0), _simulationTime(0.0), 
              _configFile(configFile),
              _simulatorID(simulatorID)
        {} 
        ~FDTD_AcousticSimulator()
        {
            if (_simulatorID) delete _simulatorID; 
        }

        inline bool CanInitializeSolver() const {return _parser && _sceneObjects && _acousticSolverSettings;}
        inline bool SceneHasModalObject() const {return _sceneObjects->HasModalObject();}
        inline bool ShouldContinue() const {return _simulationTime < _acousticSolverSettings->timeEnd;}
        inline void SetParser(const ImpulseResponseParser_Ptr &parser){_parser = parser;}
        inline void SetSolverSettings(PML_WaveSolver_Settings_Ptr rhs){_acousticSolverSettings = rhs;} 
        inline void SetSceneObjects(FDTD_Objects_Ptr objects){_sceneObjects = objects;}
        inline const std::shared_ptr<PML_WaveSolver> &GetSolver() const {return _acousticSolver;}
        inline const std::shared_ptr<PML_WaveSolver_Settings> &GetSolverSettings() const {return _acousticSolverSettings;}
        inline const std::shared_ptr<FDTD_Objects> &GetSceneObjects() const {return _sceneObjects;} 
        inline std::shared_ptr<FDTD_Objects> &GetSceneObjects(){return _sceneObjects;} 
        inline MAC_Grid &GetGrid(){return _acousticSolver->GetGrid();}
        inline REAL GetSimulationTime(){return _simulationTime;}
        inline std::string *GetSimulatorID(){return _simulatorID;}

        // parse, instance grid and solver, read mesh 
        void InitializeSolver(); // wrapper
        void InitializeSolver(const BoundingBox &solverBox, const PML_WaveSolver_Settings_Ptr &settings); 
        void ResetStartTime(const REAL &startTime); 
        bool RunForSteps(const int &N_steps); 
        void Run(); 
        bool RunHalfStep(const int &flag); 
        void PostStepping(const REAL &odeTime); 
        void EndStepping(); 
        void PreviewStepping(const uint &speed); 
        void Pause(); 
        void SaveSolverConfig(); 
        void LoadSolverResult(const std::string &dataDirectory); 

        // scene and simbox kinematics
        void AnimateObjects(const REAL newTime=-1); 
        bool MoveSimBox(const Vector3d &amount);  

        //// debug method //// 
        void TestAllComponents(); 
        void TestMoveObjects(); 
        void TestAnimateObjects(const int &N_steps); 
};
using FDTD_AcousticSimulator_Ptr = std::shared_ptr<FDTD_AcousticSimulator>; 

#endif
