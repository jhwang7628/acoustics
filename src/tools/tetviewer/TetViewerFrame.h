/*
 * =====================================================================================
 *
 *       Filename:  TetViewerFrame.h
 *
 *        Version:  1.0
 *        Created:  11/17/10 16:38:15
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  TET_VIEWER_FRAME_INC
#define  TET_VIEWER_FRAME_INC

#include <QApplication>
#include <QLabel>

#if QT_VERSION >= 0x040000
#include "ui_tetviewer.h"
#else
#error Only Qt with the version higher than 4.0 is supported right now
#endif

#include <deformable/ModeData.h>

#include <geometry/FixVtxTetMesh.hpp>

class TetViewerFrame : public QMainWindow, private Ui_TetViewerFrame
{
    friend class TetViewerCanvas;

    Q_OBJECT

    public slots:
        void open();
        void load_modes();
        void export_bin_tet();
        void export_abaqus_tet();
        void check_useless_vtx();
        void update_active_mode();
        void update_mode_displacement();

    public:
        typedef FixVtxTetMesh<double>     TMesh;

        TetViewerFrame():mesh_(NULL) 
        {
            setupUi(this);
            canvas->parent_ = this;

            QObject::connect(actionOpen, SIGNAL(triggered()), this, SLOT(open()));
            QObject::connect(actionLoadModes, SIGNAL(triggered()), this, SLOT(load_modes()));
            QObject::connect(actionWireframe, SIGNAL(toggled(bool)), canvas,
                             SLOT(toggle_wireframe(bool)));
            QObject::connect(actionMeshInfo, SIGNAL(toggled(bool)), canvas,
                             SLOT(toggle_show_info(bool)));
            QObject::connect(actionBinaryTetFormat, SIGNAL(triggered()), this,
                             SLOT(export_bin_tet()));
            QObject::connect(actionAbaqusTetFormat, SIGNAL(triggered()), this,
                             SLOT(export_abaqus_tet()));
            QObject::connect(actionCheckUselessVertex, SIGNAL(triggered()), this,
                             SLOT(check_useless_vtx()));
            QObject::connect(modeIndex, SIGNAL(valueChanged(int)), this,
                             SLOT(update_active_mode()));
            QObject::connect(modalCoordinate, SIGNAL(valueChanged(int)), this,
                             SLOT(update_mode_displacement()));
            QObject::connect(modeScale, SIGNAL(valueChanged(double)), this,
                             SLOT(update_mode_displacement()));
            QObject::connect(objectDensity, SIGNAL(valueChanged(double)), this,
                             SLOT(update_active_mode()));
            statusbar->addWidget(new QLabel("  Press \"H\" for help   ", statusbar));
        }
        ~TetViewerFrame()
        {   delete mesh_; }
        inline const ModeData &modeData() const {return modeData_;}
        inline REAL density() const {return objectDensity->value();}

    protected:

    private:
        void update_mesh(TMesh* msh);
        void update_normals();

    private:
        TMesh*                  mesh_;
        std::vector<Point3d>    vtx_;       // vertices in the tet mesh
        std::vector<Vector3d>   nml_;       // vertex normal. if the vertex is not on the surface
                                            // its normal is undefined.

        ModeData                modeData_;
        int                     activeMode_;

};

#endif   /* ----- #ifndef TET_VIEWER_FRAME_INC  ----- */

