#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMap>
#include <QCursor>
#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <geometry/Point3.hpp>
#include <ui/ModalViewer.h>
#include <utils/STL_Wrapper.h>

using namespace qglviewer;
using namespace std;

//##############################################################################
//##############################################################################
void ModalViewer::
SetAllKeyDescriptions()
{
    setKeyDescription(Qt::Key_I, "Toggle impulses display"); 
    setKeyDescription(Qt::Key_W, "Toggle wireframe-only display"); 
    setKeyDescription(Qt::Key_BracketLeft, "Previous frame (when no animation)"); 
    setKeyDescription(Qt::Key_BracketRight, "Next frame (when no animation)"); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
init()
{
    restoreStateFromFile();
    RestoreDefaultDrawOptions(); 
    glDisable(GL_LIGHTING);
    glPointSize(3.0);
    setGridIsDrawn();
    SetAllKeyDescriptions();
    PrepareImpulses(); 
    setAnimationPeriod(40); // in milliseconds

    std::cout << "\n>> Press key 'h' for help, 'esc' for exit.\n\n";
}

//##############################################################################
// This method shouldn't be called directly. Call 'updateGL()' instead. See
// documentation for QGLViewer
//##############################################################################
void ModalViewer::
draw()
{
    DrawMesh(); 
    if (_drawImpulse)
        DrawImpulses(); 
    glColor3f(0.6f, 0.6f, 0.6f); 
    if (_displayMessage)
        drawText(10, height()-20, _message); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
DrawMesh()
{
    // draw rigid sound object mesh
    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
    const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
    const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
    const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
    const int N_vertices = vertices.size(); 
    const int N_triangles = triangles.size(); 
    const int N_normals = normals.size(); 
    const REAL offsetEpsilon = 1E-5;

    // draw points
    glPointSize(3.0); 
    glBegin(GL_POINTS);
    glColor3f(0.6f, 0.6f, 0.6f); 
    for (int v_idx=0; v_idx<N_vertices; ++v_idx)
    {
        // offset the points to make it more visible
        const Point3<REAL> &vertex = vertices.at(v_idx); 
        const Vector3<REAL> &normal = normals.at(v_idx); 
        Point3<REAL> offsetVertex = vertex + normal.normalized() * offsetEpsilon;
        glVertex3f(offsetVertex.x, offsetVertex.y, offsetVertex.z);
    }
    glEnd();

    // draw edges of the triangles
    glLineWidth(1.0f); 
    glBegin(GL_LINES); 
    glColor3f(0.6f, 0.6f, 0.6f); 
    for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
    {
        const Tuple3ui &triangle = triangles.at(t_idx); 
        Point3<REAL> x = vertices.at(triangle.x); 
        Point3<REAL> y = vertices.at(triangle.y); 
        Point3<REAL> z = vertices.at(triangle.z); 
        const Vector3<REAL> &normal_x = normals.at(triangle.x); 
        const Vector3<REAL> &normal_y = normals.at(triangle.y); 
        const Vector3<REAL> &normal_z = normals.at(triangle.z); 
        x = x + normal_x.normalized() * offsetEpsilon; 
        y = y + normal_y.normalized() * offsetEpsilon; 
        z = z + normal_z.normalized() * offsetEpsilon; 
        glVertex3f(x.x, x.y, x.z); 
        glVertex3f(y.x, y.y, y.z); 

        glVertex3f(y.x, y.y, y.z); 
        glVertex3f(z.x, z.y, z.z); 

        glVertex3f(z.x, z.y, z.z); 
        glVertex3f(x.x, x.y, x.z); 
    }
    glEnd(); 

    // draw triangles
    if (!_wireframe)
    {
        glBegin(GL_TRIANGLES); 
        glColor3f(1, 1, 1); 
        for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
        {
            const Tuple3ui &triangle = triangles.at(t_idx); 
            const Point3<REAL> &x = vertices.at(triangle.x); 
            const Point3<REAL> &y = vertices.at(triangle.y); 
            const Point3<REAL> &z = vertices.at(triangle.z); 
            glVertex3f(x.x, x.y, x.z); 
            glVertex3f(y.x, y.y, y.z); 
            glVertex3f(z.x, z.y, z.z); 
        }
        glEnd(); 
    }
}

//##############################################################################
//##############################################################################
void ModalViewer::
DrawImpulses()
{
    const REAL impulseScaling = 0.01;
    const int N_frames = _rigidSoundObject->Size(); 
    _currentImpulseFrame = _currentFrame % N_frames; 
    // get impulse from object
    REAL timestamp; 
    int vertexID; 
    Vector3d impulse; 
    _rigidSoundObject->GetImpulse(_currentImpulseFrame, timestamp, vertexID, impulse); 

    // draw impulse vector
    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
    const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
    const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
    const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
    glLineWidth(3.0f); 
    glBegin(GL_LINES); 
    const Point3<REAL> &vertexEnd = vertices.at(vertexID); 
    glColor3f(1.0f, 1.0f, 0.0f);
    Point3<REAL> vertexBegin = vertexEnd - impulse * impulseScaling;
    glVertex3f(vertexBegin.x, vertexBegin.y, vertexBegin.z); 
    glVertex3f(vertexEnd.x, vertexEnd.y, vertexEnd.z); 
    glEnd(); 

    // draw impulse applied vertex
    glPointSize(10.0); 
    glBegin(GL_POINTS); 
    glColor3f(1.0f, 0.0f, 0.0f); 
    glVertex3f(vertexEnd.x, vertexEnd.y, vertexEnd.z); 
    glEnd(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
animate()
{
    if (_drawImpulse)
        DrawImpulses(); 
    _currentFrame ++; 
    PrintFrameInfo(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
keyPressEvent(QKeyEvent *e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers(); 
    bool optionsChanged = false;
    if ((e->key() == Qt::Key_I) && (modifiers == Qt::NoButton)) {
        _drawImpulse = !_drawImpulse; 
        optionsChanged = true;}
    if ((e->key() == Qt::Key_W) && (modifiers == Qt::NoButton)) {
        _wireframe = !_wireframe; 
        optionsChanged = true;}
    else if ((e->key() == Qt::Key_BracketLeft) && (modifiers == Qt::NoButton)) {
        if (!animationIsStarted())
            DrawOneFrameBackward(); }
    else if ((e->key() == Qt::Key_BracketRight) && (modifiers == Qt::NoButton)) {
        if (!animationIsStarted())
            DrawOneFrameForward(); }

    // still enable the default qglviewer event handling
    QGLViewer::keyPressEvent(e);
    if (optionsChanged) 
        PrintDrawOptions();
    updateGL();
}

//##############################################################################
//##############################################################################
QString ModalViewer::
helpString() const
{
    QString text("<h2>Modal ModalViewer</h2>");
    text += "Used for debugging FDTD_RigidSoundObject and associated classes";
    return text;
}

//##############################################################################
//##############################################################################
void ModalViewer::
DrawOneFrameForward()
{
    _currentFrame ++; 
    updateGL();
    PrintFrameInfo(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
DrawOneFrameBackward()
{
    _currentFrame --; 
    updateGL();
    PrintFrameInfo(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
PrepareImpulses()
{
    const std::string impulseFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/modalImpulses.txt"); 
    ImpulseSeriesReader reader(impulseFile); 
    std::shared_ptr<ImpulseSeriesObject> objectPtr = std::static_pointer_cast<ImpulseSeriesObject>(_rigidSoundObject); 
    reader.LoadImpulses(0, objectPtr); 
    _rigidSoundObject->GetImpulseRange(_impulseRange.start, _impulseRange.stop); 
    std::cout << "Impulses Read:\n"
              << " Number of impulses: " << _rigidSoundObject->Size() << "\n"
              << " Time range of impulses: [" << _impulseRange.start << ", " << _impulseRange.stop << "]\n"
              << "\n";
}

//##############################################################################
//##############################################################################
void ModalViewer::
RestoreDefaultDrawOptions()
{
    _drawImpulse = false; 
    _wireframe = false; 
    _displayMessage = true;
}

//##############################################################################
//##############################################################################
void ModalViewer::
PrintDrawOptions()
{
    std::cout << "\n"
              << "Draw Options \n"
              << "------------\n"
              << " Draw Impulse       : " << _drawImpulse << "\n"
              << " Draw Text Info     : " << _displayMessage << "\n"
              << " Draw Wireframe only: " << _wireframe << "\n"
              << "\n"; 
}

//##############################################################################
//##############################################################################
void ModalViewer::
PrintFrameInfo()
{
    //const std::string frameInfo("Current Frame: " + std::to_string(_currentFrame)); 
    //_message += QString::fromStdString(frameInfo); 
    _message = QString("");
    _message += "Current Frame: " + QString::number(_currentFrame) + "; "; 
    _message += "Current Impulse Frame: " + QString::number(_currentImpulseFrame); 
}
