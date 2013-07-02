/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  DemoDropObjsWithFixed.h
 *
 *        Version:  1.0
 *        Created:  01/06/12 01:18:55
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef DEMO_DROP_OBJS_WITH_FIXED_INC
#   define DEMO_DROP_OBJS_WITH_FIXED_INC

#include "demo.h"
#include "rigid/Simulator_QS.h"
#include "visual/RigidBodyRenderer.hpp"

class QGLViewer;

class DemoDropObjsWithFixed : public XDemo
{
    public:
        typedef FixVtxTetMesh<REAL>             TTetMesh;
        typedef RigidSimulator_QS::TRigidBody   TRigidBody;

        DemoDropObjsWithFixed(const char* file, QGLViewer* canvas); 

        int start();
        int step();
        void draw();

    private:
        QGLViewer*          canvas_;

        double              stepSize_;
        double              timeLen_;

        RigidSimulator_QS   rsim_;
        RigidBodyRenderer   rbRender_;
};

#endif
