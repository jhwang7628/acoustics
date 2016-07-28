#ifndef FDTD_ACOUSTIC_SIMULATOR_VIEWER_H 
#define FDTD_ACOUSTIC_SIMULATOR_VIEWER_H 
#include <memory>
#include <QGLViewer/qglviewer.h> 
#include <wavesolver/FDTD_AcousticSimulator.h>
#include <config.h>

//##############################################################################
// Class that renders Acoustic Simulator results and helps debugging.
//##############################################################################
class FDTD_AcousticSimulator_Viewer : public QGLViewer
{
    public: 
        typedef std::shared_ptr<FDTD_AcousticSimulator> SimulatorPtr; 

    private: 
        SimulatorPtr    _simulator; 
        QString         _message; 
        QString         _messageSelection;

        void SetAllKeyDescriptions(); 
        void DrawMesh(); 

    protected: 
        virtual void draw(); 
        virtual void drawWithNames(); 
        virtual void init();
        virtual void animate();
        virtual void keyPressEvent(QKeyEvent *e);
        virtual QString helpString() const;
        virtual void postSelection(const QPoint &point); 

    public: 

        FDTD_AcousticSimulator_Viewer(const std::string &simulationXMLFile)
            : _simulator(new FDTD_AcousticSimulator(simulationXMLFile))
        {
            _simulator->InitializeSolver(); 
        }

        void RestoreDefaultDrawOptions(); 
        void PrintDrawOptions(); 
};

#endif
