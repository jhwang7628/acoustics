/*
 * =====================================================================================
 *
 *       Filename:  TetViewerCanvas.cpp
 *
 *        Version:  1.0
 *        Created:  11/17/10 17:13:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "TetViewerCanvas.h"

#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include <QKeyEvent>
#include "TetViewerFrame.h"
#include "geometry/TetMesh.hpp"

using namespace std;

void TetViewerCanvas::init()
{
    init_gl();

    setKeyDescription(Qt::Key_W, "Toggles wire frame display");
    setKeyDescription(Qt::Key_I, "Toggles simulation information");
    setKeyDescription(Qt::Key_I, "Advance single step");

    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);
}

void TetViewerCanvas::init_gl()
{
#ifdef __linux
    int dummy = 0;
    glutInit(&dummy, NULL);
#endif
    const GLfloat GLOBAL_AMBIENT[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    const GLfloat SPECULAR_COLOR[] = { 0.6f, 0.6f, 0.6f, 1.0 };

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, GLOBAL_AMBIENT);
    const GLfloat ambientLight[] = { 0.f, 0.f, 0.f, 1.0f };
    const GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    const GLfloat specularLight[] = { 1.f, 1.f, 1.f, 1.0f };
    const GLfloat position[] = { -0.5f, 1.0f, 0.4f, 1.0f };

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COLOR);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    // antialiasing
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void TetViewerCanvas::toggle_wireframe(bool wf)
{
    wireframe_ = wf;
    updateGL();
}
void TetViewerCanvas::toggle_show_info(bool si)
{
    showMshInfo_ = si;
    updateGL();
}

void TetViewerCanvas::keyPressEvent(QKeyEvent* e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    
    if ( e->key() == Qt::Key_W && modifiers == Qt::NoButton ) {
        toggle_wireframe(!wireframe_);
    } else if ( e->key() == Qt::Key_I && modifiers == Qt::NoButton ) {
        printf("=========== Mesh Statistics ===========\n");
        printf(" # of vertices:          %d\n", (int)parent_->mesh_->num_vertices());
        printf(" # of tetrahedron:       %d\n", (int)parent_->mesh_->num_tets());
        printf(" # of surface triangles: %d\n", (int)parent_->mesh_->num_surface_tgls());
        printf(" # of free vertices:     %d\n", (int)parent_->mesh_->num_free_vertices());
        printf(" # of fixed vertices:    %d\n", (int)parent_->mesh_->num_fixed_vertices());
        printf("=======================================\n");
    } else {
        QGLViewer::keyPressEvent(e);
    }
}

void TetViewerCanvas::show_mesh_info()
{
    glColor3f(1.f, 1.f, 1.f);
    drawText(width()-180, 20, QString("# of vtx:       %1").arg(
                parent_->mesh_->num_vertices(), 0));
    drawText(width()-180, 40, QString("# of tets:      %1").arg(
                parent_->mesh_->num_tets(), 0));
    drawText(width()-180, 60, QString("# of surf tgl:  %1").arg(
                parent_->mesh_->num_surface_tgls(), 0));
    drawText(width()-180, 80, QString("# of free vtx:  %1").arg(
                parent_->mesh_->num_free_vertices(), 0));
    drawText(width()-180, 100, QString("# of fixed vtx: %1").arg(
                parent_->mesh_->num_fixed_vertices(), 0));
}

void TetViewerCanvas::show_mode_info()
{
    if ( modeInfo_index_ > 0 ) {
        glColor3f(1.f, 1.f, 1.f);
        drawText(width()-180, 550, QString("Mode %1: %2 Hz").arg(
                modeInfo_index_ ).arg( modeInfo_frequency_));
    }
}

void TetViewerCanvas::draw()
{
    if ( !parent_->mesh_ ) return;

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(0.f, 1.f, 0.f);
    glEnable(GL_LIGHTING);
    draw_obj();

    if ( wireframe_ )
    {
        glLineWidth(2.);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0.8f, 0.8f, 0.8f);
        glDisable(GL_LIGHTING);
        draw_obj();
    }

    if ( showMshInfo_ ) show_mesh_info();

    show_mode_info();
}

void TetViewerCanvas::draw_obj()
{
    const vector<Tuple3ui>& tgl = parent_->mesh_->surface_indices();
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(&(parent_->vtx_[0])));
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)(&(parent_->nml_[0])));
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
}

