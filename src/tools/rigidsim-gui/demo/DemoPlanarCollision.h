#ifndef DEMO_PLANAR_COLLISION_H
#define DEMO_PLANAR_COLLISION_H

#include "demo.h"
#include "rigid/Simulator.h"
#include "visual/RigidBodyRenderer.hpp"

//##############################################################################
// Forward declaration
//##############################################################################
class QGLViewer;

//##############################################################################
// This class abstracts a demo scene where multiple objects can collide on the
// 2D plane of the scene. 
//##############################################################################
class DemoPlanarCollision : public XDemo
{
    public:
        typedef FixVtxTetMesh<REAL>                 TTetMesh;
        typedef RigidBodySimulator::TRigidBody      TRigidBody;

        DemoPlanarCollision(const char* file, QGLViewer* canvas); 

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
