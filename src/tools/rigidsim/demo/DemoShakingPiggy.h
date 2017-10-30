#ifndef DEMO_SHAKING_PIGGY_H
#   define DEMO_SHAKING_PIGGY_H

#include "demo.h"
#include "rigid/Simulator.h"
#include "visual/RigidBodyRenderer.hpp"

class QGLViewer;

class DemoShakingPiggy : public XDemo
{
    public:
        typedef FixVtxTetMesh<REAL>                 TTetMesh;
        typedef RigidBodySimulator::TRigidBody      TRigidBody;

        DemoShakingPiggy(const char* file, QGLViewer* canvas);

        int start();
        int step();
        void draw();

    private:
        QGLViewer*          canvas_;

        double              stepSize_;
        double              timeLen_;

        RigidBodySimulator  rsim_;
        RigidBodyRenderer   rbRender_;
        TRigidBody*         piggy_;
        Point3d             piggyInitPos_;
        Quat4d              piggyInitRot_;
};

#endif
