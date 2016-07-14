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
        RigidSoundObjectPtr _rigidSoundObject; 
        bool                _drawImpulse; 

        void SetAllKeyDescriptions(); 

    protected: 
        virtual void draw(); 
        virtual void init();
        virtual void animate();
        virtual void keyPressEvent(QKeyEvent *e);
        virtual QString helpString() const;

    public: 

        ModalViewer(RigidSoundObjectPtr rigidSoundObject)
            : _rigidSoundObject(rigidSoundObject)
        {}

        inline void SetRigidSoundObject(RigidSoundObjectPtr &rigidSoundObject){_rigidSoundObject = rigidSoundObject;} 
        void PrepareImpulses(); 
        void RestoreDefaultDrawOptions();
        void PrintDrawOptions(); 
};

#endif 
