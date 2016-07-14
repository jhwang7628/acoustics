#ifndef MODAL_VIEWER_H 
#define MODAL_VIEWER_H 
#include <QGLViewer/qglviewer.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
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

    protected: 
        virtual void draw(); 
        virtual void init();
        virtual void animate();
        virtual QString helpString() const;

    public: 

        ModalViewer(RigidSoundObjectPtr rigidSoundObject)
            : _rigidSoundObject(rigidSoundObject)
        {}

        inline void SetRigidSoundObject(RigidSoundObjectPtr &rigidSoundObject){_rigidSoundObject = rigidSoundObject;} 
};

#endif 
