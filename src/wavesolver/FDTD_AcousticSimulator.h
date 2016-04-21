#ifndef FDTD_ACOUSTIC_SIMULATOR_H
#define FDTD_ACOUTSIC_SIMULATOR_H 

#include <string> 

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
    private: 
        // meta objects
        ImpulseResponseParser   _parser; 
        PML_WaveSolver          _acousticSolver; 
        PML_WaveSolver_Settings _acousticSolverSettings; 
        FDTD_Objects            _sceneObjects; 

        // state representation
        bool                    _canInitializeSolver; 

        // necessary fields 
        std::string             _configFile; 

    private: 
        void GetBasicSolverSettings();
        void GetOptionalSolverSettings();
        void ReadMesh(); 

    public: 
        //FDTD_AcousticSimulator(){}
        FDTD_AcousticSimulator(const std::string &configFile)
            : _configFile(configFile)
        {} 

        // parse, instance grid and solver, read mesh 
        void InitializeSolver(); 
        void Run(); 
        void Pause(); 
        void SaveSolverResult(); 
        void LoadSolverResult(const std::string &dataDirectory); 

        // GUI, OpenGL handling 
        void Draw(); 

        //// debug method //// 
};

#endif
