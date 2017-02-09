#ifndef FDTD_ACOUSTIC_SIMULATOR_H
#define FDTD_ACOUTSIC_SIMULATOR_H 

#include <string> 
#include <parser/ImpulseResponseParser.h> 
#include <wavesolver/PML_WaveSolver.h> 
#include <wavesolver/PML_WaveSolver_Settings.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/HarmonicVibrationalSource.h>
#include <wavesolver/ModalVibrationalSource.h>
#include <wavesolver/AccelerationNoiseVibrationalSource.h>
#include <wavesolver/FDTD_Objects.h> 
#include <wavesolver/FDTD_RigidObject_Animator.h>

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

        // state representation
        bool                    _canInitializeSolver; 

        // necessary fields 
        int                     _stepIndex;
        int                     _snapshotIndex;
        double                  _simulationTime; 
        std::string             _configFile; 

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
        FDTD_AcousticSimulator()
            : _canInitializeSolver(false), _stepIndex(0), _snapshotIndex(0), _simulationTime(0.0)
        {}
        FDTD_AcousticSimulator(const std::string &configFile)
            : _canInitializeSolver(false), _stepIndex(0), _snapshotIndex(0),  _simulationTime(0.0), _configFile(configFile)
        {} 

        inline bool SceneHasModalObject() const {return _sceneObjects->HasModalObject();}
        inline bool ShouldContinue() const {return _simulationTime < _acousticSolverSettings->timeEnd;}
        inline const std::shared_ptr<PML_WaveSolver> &GetSolver() const {return _acousticSolver;}
        inline const std::shared_ptr<PML_WaveSolver_Settings> &GetSolverSettings() const {return _acousticSolverSettings;}
        inline const std::shared_ptr<FDTD_Objects> &GetSceneObjects() const {return _sceneObjects;} 
        inline std::shared_ptr<FDTD_Objects> &GetSceneObjects(){return _sceneObjects;} 
        inline MAC_Grid &GetGrid(){return _acousticSolver->GetGrid();}
        inline REAL GetSimulationTime(){return _simulationTime;}

        // parse, instance grid and solver, read mesh 
        void InitializeSolver(); 
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
        void AnimateObjects(const REAL newTime=-1); 

        //// debug method //// 
        void TestAllComponents(); 
        void TestMoveObjects(); 
        void TestAnimateObjects(const int &N_steps); 
};

#endif
