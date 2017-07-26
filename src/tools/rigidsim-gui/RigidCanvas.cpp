#include "RigidCanvas.h"
#include <QKeyEvent>
#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include "demo/demo.h"

// ----------------------------------------------------------------------------
void RigidCanvas::init()
{
#ifdef __linux
    int dummy = 0;
    glutInit(&dummy, NULL);
#endif

    restoreStateFromFile();
    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);
}

void RigidCanvas::draw()
{
#ifndef NO_DRAW
    if ( pdemo_ ) pdemo_->draw();
    draw_ground();
#endif
}

void RigidCanvas::draw_ground() const
{
    const float GD_SIZE = 0.01;
    float step = GD_SIZE * 10;

    glColor3f(0.7, 0.7, 0.7);

    float d = step;
    for(int i = 0;i < 100;++ i, d += step)
    {
        glBegin(GL_LINE_LOOP);
        glVertex3f(-d, 0, -d);
        glVertex3f( d, 0, -d);
        glVertex3f( d, 0,  d);
        glVertex3f(-d, 0,  d);
        glEnd();
    }
}

void RigidCanvas::keyPressEvent(QKeyEvent* e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    const int key = e->key();

    if ( key == Qt::Key_W && modifiers == Qt::NoButton )
    {
        wireframe_ = !wireframe_;
        if ( wireframe_ )
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode(GL_FRONT, GL_FILL);
        updateGL();
    }
    else if ((e->key() == Qt::Key_T) && (modifiers == Qt::NoButton)) {
      if (camera()->type() == qglviewer::Camera::ORTHOGRAPHIC){
        camera()->setType(qglviewer::Camera::PERSPECTIVE);
        std::cout << "camera: perspective\n"; 
      }
      else {
        camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
        std::cout << "camera: orthographic\n"; 
      }
    }
    else if ( pdemo_ && pdemo_->key_pressed(e) ) 
        updateGL();
    else
        QGLViewer::keyPressEvent(e);
}
