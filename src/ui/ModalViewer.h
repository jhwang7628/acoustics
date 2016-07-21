#ifndef MODAL_VIEWER_H 
#define MODAL_VIEWER_H 
#include <Eigen/Dense> 
#include <QGLViewer/qglviewer.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <io/ImpulseSeriesReader.h>
#include <modal_model/ImpulseSeriesObject.h>
#include <config.h>

//##############################################################################
// Class that helps debugging FDTD_RigidSoundObject and associated classes
//##############################################################################
class ModalViewer : public QGLViewer
{
    public: 
        typedef std::shared_ptr<FDTD_RigidSoundObject> RigidSoundObjectPtr; 

    private: 
        struct TimeInterval{ REAL start; REAL stop; }; 

        RigidSoundObjectPtr _rigidSoundObject; 
        TimeInterval        _impulseRange; 
        bool                _drawImpulse; 
        bool                _displayMessage; 
        int                 _wireframe;  // 0: both wire and face; 1: wire; 2: face
        int                 _drawModes; 
        int                 _currentFrame; 
        int                 _currentImpulseFrame; 
        REAL                _timeStepSize; 
        REAL                _ODEStepSize; 
        REAL                _startTime; 
        QString             _message; 
        QString             _messageSelection; 
        Eigen::VectorXd     _vertexValues; 

        void SetAllKeyDescriptions(); 
        void DrawMesh(); 
        void DrawImpulses();
        void UpdateVertexValues(); 

    protected: 
        virtual void draw(); 
        virtual void drawWithNames(); 
        virtual void init();
        virtual void animate();
        virtual void keyPressEvent(QKeyEvent *e);
        virtual QString helpString() const;
        virtual void postSelection(const QPoint &point); 

        virtual void DrawOneFrameForward(); 
        virtual void DrawOneFrameBackward(); 

    public: 

        ModalViewer(RigidSoundObjectPtr rigidSoundObject)
            : _rigidSoundObject(rigidSoundObject), _currentFrame(0), _currentImpulseFrame(0), 
              _timeStepSize(0.001), _startTime(0.0)
        {}

        inline REAL CurrentTime(){return _startTime + _currentFrame * _timeStepSize;} 
        inline void SetRigidSoundObject(RigidSoundObjectPtr &rigidSoundObject){_rigidSoundObject = rigidSoundObject;} 
        void PrepareImpulses(); 
        void PrepareModes(); 
        void RestoreDefaultDrawOptions();
        void StepODEAndStoreResults(); 
        void PrintDrawOptions(); 
        void PrintFrameInfo(); 
};

#endif 
