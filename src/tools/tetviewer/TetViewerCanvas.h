/*
 * =====================================================================================
 *
 *       Filename:  TetViewerCanvas.h
 *
 *        Version:  1.0
 *        Created:  11/17/10 16:41:25
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  TET_VIEWER_CANVAS_INC
#define  TET_VIEWER_CANVAS_INC

#include <QGLViewer/qglviewer.h>

class TetViewerFrame;

class TetViewerCanvas : public QGLViewer
{
    friend class TetViewerFrame;

    Q_OBJECT

    public slots:
        void toggle_wireframe(bool wf);
        void toggle_show_info(bool si);

    public:
        TetViewerCanvas(QWidget* parent):QGLViewer(parent), 
                wireframe_(false), showMshInfo_(false),
                modeInfo_index_(0)
        { }

        void update_mode_info(int modeIndex, double modeFrequency)
        {
            modeInfo_index_ = modeIndex;
            modeInfo_frequency_ = modeFrequency;
        }

    protected:
        void init();
        void draw();
        void keyPressEvent(QKeyEvent* e);

    private:
        void init_gl();
        void draw_obj();
        void show_mesh_info();
        void show_mode_info();

    private:
        bool            wireframe_;
        bool            showMshInfo_;

        int             modeInfo_index_;
        double          modeInfo_frequency_;

        TetViewerFrame* parent_;
};

#endif   /* ----- #ifndef TET_VIEWER_CANVAS_INC  ----- */
