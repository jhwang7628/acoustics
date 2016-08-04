#ifndef FDTD_ACOUSTIC_SIMULATOR_VIEWER_H 
#define FDTD_ACOUSTIC_SIMULATOR_VIEWER_H 
#include <memory>
#include <QGLViewer/qglviewer.h> 
#include <wavesolver/FDTD_AcousticSimulator.h>
#include <wavesolver/FDTD_RigidObject_Animator.h> 
#include <linearalgebra/Vector3.hpp>
#include <config.h>

//##############################################################################
// Class that renders Acoustic Simulator results and helps debugging.
//##############################################################################
class FDTD_AcousticSimulator_Viewer : public QGLViewer
{
    public: 
        typedef std::shared_ptr<FDTD_AcousticSimulator> SimulatorPtr; 
        struct Arrow{Vector3f start; Vector3f normal;}; 
        struct Slice{int dim; Vector3d origin; Vector3Array samples; Vector3Array gridLines; int N_sample_per_dim; REAL minBound; REAL maxBound;};

    private: 
        SimulatorPtr    _simulator; 

        int             _currentFrame = 0;
        QString         _message; 
        QString         _messageSelection;
        int             _wireframe;
        int             _sliceWireframe;
        bool            _drawBox; 
        std::vector<Vector3f> _objectColors; 
        std::vector<Vector3f> _sphereCin; 
        std::vector<Arrow> _arrowCin; 
        std::vector<Slice> _sliceCin; 
        REAL            _drawAbsMax; 
        int             _sliceDataPointer; // 0: pressure; 1: cell id

        void SetAllKeyDescriptions(); 
        void DrawMesh(); 
        void DrawBox(); 
        void DrawListeningPoints(); 
        void DrawLights(); 
        void DrawSlices(const int &dataPointer);
        void DrawDebugCin();

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
        }

        void ConstructSliceSamples(Slice &slice); 
        void DrawOneFrameForward(); 
        void RestoreDefaultDrawOptions(); 
        void PrintFrameInfo(); 
        void PrintDrawOptions(); 
};

#endif
