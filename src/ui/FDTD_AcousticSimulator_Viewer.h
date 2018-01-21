#ifndef FDTD_ACOUSTIC_SIMULATOR_VIEWER_H
#define FDTD_ACOUSTIC_SIMULATOR_VIEWER_H
#include <memory>
#include <QWidget>
#include <QGLViewer/qglviewer.h>
#include <wavesolver/FDTD_AcousticSimulator.h>
#include <wavesolver/FDTD_RigidObject_Animator.h>
#include <wavesolver/PML_WaveSolver_Settings.h>
#include <wavesolver/SimWorld.h>
#include <modal_model/KirchhoffIntegralSolver.h>
#include <linearalgebra/Vector3.hpp>
#include <linearalgebra/Vector2.hpp>
#include <colormap/ColorMap.h>
#include <config.h>

//##############################################################################
// Forward declaration
//##############################################################################
class FDTD_AcousticSimulator_Widget;

//##############################################################################
// Class that renders Acoustic Simulator results and helps debugging.
//##############################################################################
class FDTD_AcousticSimulator_Viewer : public QGLViewer
{
    public:
        typedef std::shared_ptr<KirchhoffIntegralSolver> BEMSolverPtr;
        struct Arrow{Vector3f start; Vector3f normal;};
        struct Sphere{Vector3f origin; REAL scale;};
        struct Slice
        {
            int dim;
            Vector3d origin;
            Vector3Array samples;
            Eigen::MatrixXd data;
            Vector3Array gridLines;
            Tuple3i N_sample_per_dim;
            Vector3d minBound;
            Vector3d maxBound;
            std::vector<MAC_Grid::Cell> cells;
            bool dataReady = false;
            ActiveSimUnit_Ptr intersectingUnit;
        };
        struct TextureMeta
        {
            QString name
                = "/home/jui-hsien/code/acoustics/src/ui/assets/checkerboard.png";
            float blockSizeUV = 0.2;
            bool groundLoaded = false;
        } textureMeta;

    private:
        SimWorld_UPtr                            _simWorld;
        std::shared_ptr<PML_WaveSolver_Settings> _solverSettings;

        uint                    _previewSpeed = 0;
        bool                    _remoteConnection = false;
        int                     _halfStepFlag = 0;
        int                     _currentFrame = 0;
        QString                 _message;
        QString                 _messageSelection;
        QString                 _messageColormap;
        int                     _wireframe;
        std::bitset<2>          _sliceWireframe;
        bool                    _drawBoxLis;
        int                     _drawGround;
        bool                    _drawHashedCells;
        bool                    _takeSnapshots = false;
        std::vector<Vector3f>   _objectColors;
        std::vector<Sphere>     _sphereCin;
        std::vector<Arrow>      _arrowCin;
        std::vector<Slice>      _sliceCin;
        REAL                    _drawAbsMax;
        REAL                    _drawImpulseScaling = 1.;
        ActiveSimUnit_Ptr       _listenedUnit;
        MAC_Grid::Cell          _listenedCell;

        // slice related fields
        int                         _sliceDataPointer; // 0: pressure; 1: cell id; 2: vx; 3: vy; 4:vz; 5: p_x; 6: p_y; 7: p_z; 8: frequency transfer; 9: frequency transfer residual
        std::shared_ptr<ColorMap>   _sliceColorMap;
        Vector2d                    _sliceColorMapRange;
        int                         _sliceDivision = -1;
        bool                        _fixedSliceColorMapRange = false;
        int                         _meshDataPointer = 0; // 0: nothing; 1: curvature; 2: surface acceleration

        // frequency transfer solver
        BEMSolverPtr    _bemSolver;
        int             _bemModePointer = 0;

        // scene parameters
        BoundingBox *_sceneBox = nullptr;

        void SetAllKeyDescriptions();
        void UpdateVertexShader();
        void DrawMesh();
        void DrawImpulses();
        void DrawBox();
        void LoadGround();
        void DrawGround();
        void DrawListeningPoints();
        void DrawSelection();
        void DrawLights();
        void DrawSlices(const int &dataPointer);
        void DrawDebugCin();
        void DrawHashedCells();
        inline void ResetSliceColormap(){_sliceColorMap->set_interpolation_range(_sliceColorMapRange.x, _sliceColorMapRange.y);}

        void updateGL() {update();}

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
        FDTD_AcousticSimulator_Viewer(SimWorld_UPtr world)
            : _simWorld(std::move(world))
        {
            RestoreDefaultDrawOptions();
            _solverSettings = _simWorld->GetSolverSettings();
            //if (_simulator->SceneHasModalObject() && _solverSettings->validateUsingFBem)
            //    InitializeBEMSolver();
        }
        ~FDTD_AcousticSimulator_Viewer()
        {
            if (_sceneBox) delete _sceneBox;
        }

        inline void SetAllSliceDataReady(const bool &isReady)
        {for (auto &slice : _sliceCin) slice.dataReady = isReady;}
        inline void SetPreviewSpeed(const uint &speed)
        {
            _previewSpeed = speed;
            _solverSettings->isPreviewMode = true;
        }
        void InitializeBEMSolver();
        void AddSlice(const int &dim, const REAL &offset);
        void ConstructSliceSamples(Slice &slice);
        void ComputeAndCacheSliceData(const int &dataPointer);
        void DrawOneFrameForward();
        void RestoreDefaultDrawOptions();
        void PrintFrameInfo();
        void PrintDrawOptions();
        void Push_Back_ReflectionArrows(const std::string &filename);
        void MoveSceneCenter(const int &dim, const double &displacement);

    friend FDTD_AcousticSimulator_Widget;
};

#endif
