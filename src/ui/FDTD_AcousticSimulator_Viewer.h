#ifndef FDTD_ACOUSTIC_SIMULATOR_VIEWER_H 
#define FDTD_ACOUSTIC_SIMULATOR_VIEWER_H 
#include <memory>
#include <QGLViewer/qglviewer.h> 
#include <wavesolver/FDTD_AcousticSimulator.h>
#include <wavesolver/FDTD_RigidObject_Animator.h> 
#include <modal_model/KirchhoffIntegralSolver.h>
#include <linearalgebra/Vector3.hpp>
#include <linearalgebra/Vector2.hpp>
#include <colormap/ColorMap.h>
#include <config.h>

//##############################################################################
// Class that renders Acoustic Simulator results and helps debugging.
//##############################################################################
class FDTD_AcousticSimulator_Viewer : public QGLViewer
{
    public: 
        typedef std::shared_ptr<FDTD_AcousticSimulator> SimulatorPtr; 
        typedef std::shared_ptr<KirchhoffIntegralSolver> BEMSolverPtr; 
        struct Arrow{Vector3f start; Vector3f normal;}; 
        struct Sphere{Vector3f origin; REAL scale;};
        struct Slice{int dim; Vector3d origin; Vector3Array samples; Eigen::MatrixXd data; Vector3Array gridLines; int N_sample_per_dim; Vector3d minBound; Vector3d maxBound; std::vector<MAC_Grid::Cell> cells; bool dataReady = false;};

    private: 
        SimulatorPtr            _simulator; 
        int                     _halfStepFlag = 0;

        bool                    _remoteConnection = false; 
        int                     _currentFrame = 0;
        QString                 _message; 
        QString                 _messageSelection;
        QString                 _messageColormap; 
        int                     _wireframe;
        int                     _sliceWireframe;
        bool                    _drawBox; 
        std::vector<Vector3f>   _objectColors; 
        std::vector<Sphere>     _sphereCin; 
        std::vector<Arrow>      _arrowCin; 
        std::vector<Slice>      _sliceCin; 
        REAL                    _drawAbsMax; 
        MAC_Grid::Cell          _listenedCell; 

        // slice related fields
        int                         _sliceDataPointer; // 0: pressure; 1: cell id; 2: vx; 3: vy; 4:vz; 5: p_x; 6: p_y; 7: p_z; 8: frequency transfer; 9: frequency transfer residual
        std::shared_ptr<ColorMap>   _sliceColorMap; 
        Vector2d                    _sliceColorMapRange; 
        int                         _sliceDivision = 80; 
        bool                        _fixedSliceColorMapRange = false; 

        int                         _meshDataPointer = 0; // 0: nothing; 1: curvature

        // frequency transfer solver
        BEMSolverPtr    _bemSolver;
        int             _bemModePointer = 0; 

        void SetAllKeyDescriptions(); 
        void DrawMesh(); 
        void DrawImpulses(); 
        void DrawBox(); 
        void DrawListeningPoints(); 
        void DrawSelection(); 
        void DrawLights(); 
        void DrawSlices(const int &dataPointer);
        void DrawDebugCin();
        inline void ResetSliceColormap(){_sliceColorMap->set_interpolation_range(_sliceColorMapRange.x, _sliceColorMapRange.y);}

    protected: 
        virtual void draw(); 
        virtual void drawWithNames(); 
        virtual void init();
        virtual void init_gl();
        virtual void animate();
        virtual void keyPressEvent(QKeyEvent *e);
        virtual QString helpString() const;
        virtual void postSelection(const QPoint &point); 

    public: 

        FDTD_AcousticSimulator_Viewer(const std::string &simulationXMLFile)
            : _simulator(new FDTD_AcousticSimulator(simulationXMLFile))
        {
            RestoreDefaultDrawOptions();
            _simulator->InitializeSolver(); 
            if (_simulator->SceneHasModalObject() && _simulator->GetSolverSettings()->validateUsingFBem)
                InitializeBEMSolver();
        }

        inline void SetAllSliceDataReady(const bool &isReady){for (auto &slice : _sliceCin) slice.dataReady = isReady;}
        void InitializeBEMSolver(); 
        void ConstructSliceSamples(Slice &slice); 
        void ComputeAndCacheSliceData(const int &dataPointer, Slice &slice); 
        void DrawOneFrameForward(); 
        void DrawHalfFrameForward(); 
        void RestoreDefaultDrawOptions(); 
        void PrintFrameInfo(); 
        void PrintDrawOptions(); 
        void Push_Back_ReflectionArrows(const std::string &filename); 
};

#endif
