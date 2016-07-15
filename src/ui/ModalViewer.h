#ifndef MODAL_VIEWER_H 
#define MODAL_VIEWER_H 
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
        bool                _wireframe;
        int                 _currentFrame; 
        int                 _currentImpulseFrame; 
        QString             _message; 

        void SetAllKeyDescriptions(); 
        void DrawMesh(); 
        void DrawImpulses();

    protected: 
        virtual void draw(); 
        virtual void init();
        virtual void animate();
        virtual void keyPressEvent(QKeyEvent *e);
        virtual QString helpString() const;

        virtual void DrawOneFrameForward(); 
        virtual void DrawOneFrameBackward(); 

    public: 

        ModalViewer(RigidSoundObjectPtr rigidSoundObject)
            : _rigidSoundObject(rigidSoundObject), _currentFrame(0), _currentImpulseFrame(0)
        {}

        inline void SetRigidSoundObject(RigidSoundObjectPtr &rigidSoundObject){_rigidSoundObject = rigidSoundObject;} 
        void PrepareImpulses(); 
        void PrepareModes(); 
        void RestoreDefaultDrawOptions();
        void PrintDrawOptions(); 
        void PrintFrameInfo(); 
};

#endif 
