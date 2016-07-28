#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMap>
#include <QCursor>
#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <geometry/Point3.hpp>

using namespace qglviewer;

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
SetAllKeyDescriptions()
{
    //setKeyDescription(Qt::Key_I, "Toggle impulses display"); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
init()
{
    restoreStateFromFile();
    RestoreDefaultDrawOptions(); 
    glDisable(GL_LIGHTING);
    glPointSize(3.0);
    setGridIsDrawn();
    SetAllKeyDescriptions();
    setAnimationPeriod(40); // in milliseconds

    std::cout << "\n>> Press key 'h' for help, 'esc' for exit.\n\n";
}

//##############################################################################
// This method shouldn't be called directly. Call 'updateGL()' instead. See
// documentation for QGLViewer
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
draw()
{
    DrawMesh(); 
    glColor3f(0.6f, 0.6f, 0.6f); 
    drawText(10, height()-20, _message); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
drawWithNames()
{
//    //// draw rigid sound object mesh
//    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
//    const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
//    const int N_vertices = vertices.size(); 
//    const REAL ballSize = 3E-4;
//
//    // draw points
//    glPointSize(3.0); 
//    for (int v_idx=0; v_idx<N_vertices; ++v_idx)
//    {
//        const Point3<REAL> &vertex = vertices.at(v_idx); 
//        glColor3f(0.6f, 0.6f, 0.6f); 
//        glPushMatrix();
//        glTranslatef(vertex.x, vertex.y, vertex.z); 
//        glPushName(v_idx);
//        GL_Wrapper::DrawSphere(ballSize, 3, 3);
//        glPopName();
//        glPopMatrix();
//    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawMesh()
{
//    bool isDrawModes = false; 
//    //bool isDrawModes = true;
//    //if (_drawModes < 0 || _vertexValues.size() == 0) 
//    //    isDrawModes = false;
//
//    // draw rigid sound object mesh
//    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
//    const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
//    const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
//    const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
//    const int N_triangles = triangles.size(); 
//    const REAL offsetEpsilon = 1E-5;
//    const REAL ballSize = 3E-4;
//
//    // draw points only if selected
//    if (selectedName() != -1)
//    {
//        const int v_idx = selectedName(); 
//        // offset the points to make it more visible
//        const Point3<REAL> &vertex = vertices.at(v_idx); 
//        //const Vector3<REAL> &normal = normals.at(v_idx); 
//        glPointSize(10.0);
//        glColor3f(0.0f, 1.0f, 0.0f); 
//        glPushMatrix(); 
//        glTranslatef(vertex.x, vertex.y, vertex.z); 
//        GL_Wrapper::DrawSphere(ballSize, 3, 3);
//        glPopMatrix(); 
//    }
//
//    // draw edges of the triangles
//    if (_wireframe == 0 || _wireframe == 1)
//    {
//        glLineWidth(1.0f); 
//        glBegin(GL_LINES); 
//        glColor3f(0.6f, 0.6f, 0.6f); 
//        for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
//        {
//            const Tuple3ui &triangle = triangles.at(t_idx); 
//            Point3<REAL> x = vertices.at(triangle.x); 
//            Point3<REAL> y = vertices.at(triangle.y); 
//            Point3<REAL> z = vertices.at(triangle.z); 
//            const Vector3<REAL> &normal_x = normals.at(triangle.x); 
//            const Vector3<REAL> &normal_y = normals.at(triangle.y); 
//            const Vector3<REAL> &normal_z = normals.at(triangle.z); 
//            x = x + normal_x.normalized() * offsetEpsilon; 
//            y = y + normal_y.normalized() * offsetEpsilon; 
//            z = z + normal_z.normalized() * offsetEpsilon; 
//            glVertex3f(x.x, x.y, x.z); 
//            glVertex3f(y.x, y.y, y.z); 
//
//            glVertex3f(y.x, y.y, y.z); 
//            glVertex3f(z.x, z.y, z.z); 
//
//            glVertex3f(z.x, z.y, z.z); 
//            glVertex3f(x.x, x.y, x.z); 
//        }
//        glEnd(); 
//    }
//
//    // draw triangles
//    if (_wireframe == 0 || _wireframe == 2)
//    {
//        glBegin(GL_TRIANGLES); 
//        glColor3f(1, 1, 1); 
//        for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
//        {
//            const Tuple3ui &triangle = triangles.at(t_idx); 
//            const Point3<REAL> &x = vertices.at(triangle.x); 
//            const Point3<REAL> &y = vertices.at(triangle.y); 
//            const Point3<REAL> &z = vertices.at(triangle.z); 
//            if (isDrawModes)
//            {
//                const REAL xValue = _vertexValues(triangle.x);
//                const REAL yValue = _vertexValues(triangle.y);
//                const REAL zValue = _vertexValues(triangle.z);
//                glColor3f(xValue, 0, -xValue);
//                glVertex3f(x.x, x.y, x.z); 
//                glColor3f(yValue, 0, -xValue);
//                glVertex3f(y.x, y.y, y.z); 
//                glColor3f(zValue, 0, -xValue);
//                glVertex3f(z.x, z.y, z.z); 
//            }
//            else
//            {
//                glVertex3f(x.x, x.y, x.z); 
//                glVertex3f(y.x, y.y, y.z); 
//                glVertex3f(z.x, z.y, z.z); 
//            }
//        }
//        glEnd(); 
//    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
animate()
{
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
keyPressEvent(QKeyEvent *e)
{
//    const Qt::KeyboardModifiers modifiers = e->modifiers(); 
//    bool optionsChanged = false;
    bool handled = false; 
//    if ((e->key() == Qt::Key_I) && (modifiers == Qt::NoButton)) {
//        _drawImpulse = !_drawImpulse; 
//        optionsChanged = true;}
//    if ((e->key() == Qt::Key_W) && (modifiers == Qt::NoButton)) {
//        _wireframe = (_wireframe+1)%3; 
//        optionsChanged = true;}
//    else if ((e->key() == Qt::Key_BracketLeft) && (modifiers == Qt::NoButton)) {
//        if (!animationIsStarted())
//            DrawOneFrameBackward(); }
//    else if ((e->key() == Qt::Key_BracketRight) && (modifiers == Qt::NoButton)) {
//        if (!animationIsStarted())
//            DrawOneFrameForward(); }
//    else if ((e->key() == Qt::Key_M) && (modifiers == Qt::ShiftModifier)){
//        _drawModes--; 
//        UpdateVertexValues();
//        optionsChanged = true;}
//    else if ((e->key() == Qt::Key_M) && (modifiers == Qt::NoButton)){
//        _drawModes++; 
//        UpdateVertexValues();
//        optionsChanged = true;}
//    else if ((e->key() == Qt::Key_N) && (modifiers == Qt::NoButton)){
//        _rigidSoundObject->AdvanceModalODESolvers(1);
//        UpdateVertexValues();}
//    else if ((e->key() == Qt::Key_R) && (modifiers == Qt::NoButton)){
//        StepODEAndStoreResults();}
//    else if ((e->key() == Qt::Key_F) && (modifiers == Qt::NoButton)){
//        PrintAllFrequencies(std::cout);}
//    else {
//        handled = false;}
//
//
//    // still enable the default qglviewer event handling but this function has
//    // priority
    if (!handled)
        QGLViewer::keyPressEvent(e);
//    if (optionsChanged) 
//        PrintDrawOptions();
    updateGL();
}

//##############################################################################
//##############################################################################
QString FDTD_AcousticSimulator_Viewer::
helpString() const
{
    QString text("<h2>Modal FDTD_AcousticSimulator_Viewer</h2>");
    text += "Used for debugging FDTD_AcousticSimulator and associated classes";
    return text;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
postSelection(const QPoint &point)
{
    _messageSelection = QString("Vertex ID= " + QString::number(selectedName()));
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
RestoreDefaultDrawOptions()
{
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
PrintDrawOptions()
{
    std::cout << "\n"
              << "Draw Options \n"
              << "------------\n"
              << "\n"; 
}

