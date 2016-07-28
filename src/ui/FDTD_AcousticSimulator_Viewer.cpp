#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMap>
#include <QCursor>
#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/FDTD_Objects.h>
#include <wavesolver/PML_WaveSolver_Settings.h>
#include <utils/GL_Wrapper.h>
#include <config.h>

using namespace qglviewer;

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
SetAllKeyDescriptions()
{
    setKeyDescription(Qt::Key_W, "Toggle wireframe-only display"); 
    setKeyDescription(Qt::Key_B, "Toggle simulation box display"); 
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
    //setGridIsDrawn();
    SetAllKeyDescriptions();
    setAnimationPeriod(40); // in milliseconds

    init_gl();
    std::cout << "\n>> Press key 'h' for help, 'esc' for exit.\n\n";
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
init_gl()
{
    // fetch objects to set up for colors, material for gl rendering
    const int N_objects = _simulator->GetSceneObjects()->GetRigidSoundObjects().size(); 
    _objectColors.resize(N_objects); 
    for (int obj_idx=0; obj_idx<N_objects; ++obj_idx)
    {
        const float x = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) *0.6f + 0.4f;
        const float y = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) *0.6f + 0.4f;
        const float z = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) *0.6f + 0.4f;
        _objectColors[obj_idx] = Vector3f(x, y, z);
    }

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

    if (_drawBox)
        DrawBox(); 
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
    const auto &sceneObjects = _simulator->GetSceneObjects(); 
    const auto &rigidSoundObjects = sceneObjects->GetRigidSoundObjects();
    const int N_objects = rigidSoundObjects.size(); 
    for (int obj_idx=0; obj_idx<N_objects; ++obj_idx)
    {
        bool isDrawModes = false; 
        //bool isDrawModes = true;
        //if (_drawModes < 0 || _vertexValues.size() == 0) 
        //    isDrawModes = false;

        // draw rigid sound object mesh
        auto &object = rigidSoundObjects.at(obj_idx); 
        std::shared_ptr<TriangleMesh<REAL> > meshPtr = object->GetMeshPtr();
        const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
        const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
        const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
        const int N_triangles = triangles.size(); 
        const REAL offsetEpsilon = 1E-5;
        const REAL ballSize = 3E-4;

        // draw points only if selected
        if (selectedName() != -1)
        {
            const int v_idx = selectedName(); 
            // offset the points to make it more visible
            const Point3<REAL> &vertex = vertices.at(v_idx); 
            //const Vector3<REAL> &normal = normals.at(v_idx); 
            glPointSize(10.0);
            glColor3f(0.0f, 1.0f, 0.0f); 
            glPushMatrix(); 
            glTranslatef(vertex.x, vertex.y, vertex.z); 
            GL_Wrapper::DrawSphere(ballSize, 3, 3);
            glPopMatrix(); 
        }

        // draw edges of the triangles
        if (_wireframe == 0 || _wireframe == 1)
        {
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
        }

        // draw triangles
        glEnable(GL_LIGHTING);
        if (_wireframe == 0 || _wireframe == 2)
        {
            glBegin(GL_TRIANGLES); 
            const auto &color = _objectColors.at(obj_idx); 
            glColor3f(color.x, color.y, color.z); 
            for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
            {
                const Tuple3ui &triangle = triangles.at(t_idx); 
                const Point3<REAL> &x = vertices.at(triangle.x); 
                const Point3<REAL> &y = vertices.at(triangle.y); 
                const Point3<REAL> &z = vertices.at(triangle.z); 
                const Vector3<REAL> &nx = normals.at(triangle.x); 
                const Vector3<REAL> &ny = normals.at(triangle.y); 
                const Vector3<REAL> &nz = normals.at(triangle.z); 
                if (isDrawModes)
                {
                    //const REAL xValue = _vertexValues(triangle.x);
                    //const REAL yValue = _vertexValues(triangle.y);
                    //const REAL zValue = _vertexValues(triangle.z);
                    //glColor3f(xValue, 0, -xValue);
                    //glVertex3f(x.x, x.y, x.z); 
                    //glColor3f(yValue, 0, -xValue);
                    //glVertex3f(y.x, y.y, y.z); 
                    //glColor3f(zValue, 0, -xValue);
                    //glVertex3f(z.x, z.y, z.z); 
                }
                else
                {
                    glNormal3f(nx.x, nx.y, nx.z); 
                    glVertex3f(x.x, x.y, x.z); 
                    glNormal3f(ny.x, ny.y, ny.z); 
                    glVertex3f(y.x, y.y, y.z); 
                    glNormal3f(nz.x, nz.y, nz.z); 
                    glVertex3f(z.x, z.y, z.z); 
                }
            }
            glEnd(); 
        }
        glDisable(GL_LIGHTING);
    }
}

//##############################################################################
// Draw simulation box
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawBox()
{
    const auto &settings = _simulator->GetSolverSettings(); 
    const REAL cellSize = settings->cellSize; 
    const int division = settings->cellDivisions; 
    const REAL halfLength = (REAL)division*cellSize / 2.0; 
    const double minBound[3] = {-halfLength, -halfLength, -halfLength}; 
    const double maxBound[3] = { halfLength,  halfLength,  halfLength}; 
    glColor3f(1.0f, 1.0f, 1.0f);
    GL_Wrapper::DrawWireBox(&minBound[0], &maxBound[0]); 
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
    const Qt::KeyboardModifiers modifiers = e->modifiers(); 
    bool optionsChanged = false;
    bool handled = true; 
    if ((e->key() == Qt::Key_W) && (modifiers == Qt::NoButton)) {
        _wireframe = (_wireframe+1)%3; 
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_B) && (modifiers == Qt::NoButton)) {
        _drawBox = !_drawBox; 
        optionsChanged = true;
    }
    else {
        handled = false;
    }

    // still enable the default qglviewer event handling but this function has
    // priority
    if (!handled)
        QGLViewer::keyPressEvent(e);
    if (optionsChanged) 
        PrintDrawOptions();
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
    _wireframe = 2;
    _drawBox = true; 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
PrintDrawOptions()
{
    std::cout << "\n"
              << "Draw Options \n"
              << "------------\n"
              << " Draw simulation box: " << _drawBox << "\n"
              << " Draw wireframe only: " << _wireframe << "\n"
              << "\n"; 
}

