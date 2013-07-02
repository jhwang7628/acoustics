#ifndef DEMO_DROP_OBJS_H
#   define DEMO_DROP_OBJS_H

#include "demo.h"
#include "rigid/Simulator.h"
#include "visual/RigidBodyRenderer.hpp"

class QGLViewer;

class DemoDropObjs : public XDemo
{
    public:
        typedef FixVtxTetMesh<REAL>                 TTetMesh;
        typedef RigidBodySimulator::TRigidBody      TRigidBody;

        DemoDropObjs(const char* file, QGLViewer* canvas); 

        int start();
        int step();
        void draw();

    private:
        QGLViewer*          canvas_;

        double              stepSize_;
        double              timeLen_;

        RigidBodySimulator  rsim_;
        RigidBodyRenderer   rbRender_;
};

#endif
