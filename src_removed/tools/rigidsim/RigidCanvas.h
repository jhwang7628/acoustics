#ifndef UI_RIGIDSIM_CANVAS_H
#   define UI_RIGIDSIM_CANVAS_H

#include <QGLViewer/qglviewer.h>

class XDemo;
class RigidSim;

class RigidCanvas : public QGLViewer
{
    friend class RigidSim;

    public:
        RigidCanvas(QWidget* parent):QGLViewer(parent), 
                pdemo_(NULL), wireframe_(false)
        { }

    protected:
        void draw();
        void init();
        void keyPressEvent(QKeyEvent* e);

    private:
        void draw_ground() const;

    private:
        XDemo*      pdemo_;
        bool        wireframe_;
};

#endif
