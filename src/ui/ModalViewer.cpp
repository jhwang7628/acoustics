#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMap>
#include <QCursor>
#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <macros.h>
#include <geometry/Point3.hpp>
#include <ui/ModalViewer.h>
#include <utils/IO/IO.h>
#include <utils/STL_Wrapper.h>
#include <utils/GL_Wrapper.h>
#include <parser/ImpulseResponseParser.h>

using namespace qglviewer;
using namespace std;

//##############################################################################
//##############################################################################
void ModalViewer::
SetAllKeyDescriptions()
{
    setKeyDescription(Qt::Key_I, "Toggle impulses display"); 
    setKeyDescription(Qt::Key_W, "Toggle wireframe-only display"); 
    setKeyDescription(Qt::Key_BracketLeft, "Previous impulse frame (when no animation)"); 
    setKeyDescription(Qt::Key_BracketRight, "Next impulse frame (when no animation)"); 
    setKeyDescription(Qt::Key_M, "Next mode"); 
    setKeyDescription(Qt::ShiftModifier+Qt::Key_M, "Previous mode"); 
    setKeyDescription(Qt::Key_N, "Step Modal ODE"); 
    setKeyDescription(Qt::Key_F, "Print all modal frequencies"); 
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
    PrepareModes(); 
    _rigidSoundObject->FDTD_RigidSoundObject::Initialize(); 
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
    {
        drawText(10, height()-20, _message); 
        drawText(10, height()-40, _messageSelection);
    }
}

//##############################################################################
//##############################################################################
void ModalViewer::
drawWithNames()
{
    // draw rigid sound object mesh
    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
    const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
    const int N_vertices = vertices.size(); 
    const REAL ballSize = 3E-4;

    // draw points
    glPointSize(3.0); 
    for (int v_idx=0; v_idx<N_vertices; ++v_idx)
    {
        const Point3<REAL> &vertex = vertices.at(v_idx); 
        glColor3f(0.6f, 0.6f, 0.6f); 
        glPushMatrix();
        glTranslatef(vertex.x, vertex.y, vertex.z); 
        glPushName(v_idx);
        GL_Wrapper::DrawSphere(ballSize, 3, 3);
        glPopName();
        glPopMatrix();
    }
}

//##############################################################################
//##############################################################################
void ModalViewer::
DrawMesh()
{
    bool isDrawModes = true;
    if (_drawModes < 0 || _vertexValues.size() == 0) 
        isDrawModes = false;

    // draw rigid sound object mesh
    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
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
    if (_wireframe == 0 || _wireframe == 2)
    {
        glBegin(GL_TRIANGLES); 
        glColor3f(1, 1, 1); 
        for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
        {
            const Tuple3ui &triangle = triangles.at(t_idx); 
            const Point3<REAL> &x = vertices.at(triangle.x); 
            const Point3<REAL> &y = vertices.at(triangle.y); 
            const Point3<REAL> &z = vertices.at(triangle.z); 
            if (isDrawModes)
            {
                const REAL xValue = _vertexValues(triangle.x);
                const REAL yValue = _vertexValues(triangle.y);
                const REAL zValue = _vertexValues(triangle.z);
                glColor3f(xValue, 0, -xValue);
                glVertex3f(x.x, x.y, x.z); 
                glColor3f(yValue, 0, -xValue);
                glVertex3f(y.x, y.y, y.z); 
                glColor3f(zValue, 0, -xValue);
                glVertex3f(z.x, z.y, z.z); 
            }
            else
            {
                glVertex3f(x.x, x.y, x.z); 
                glVertex3f(y.x, y.y, y.z); 
                glVertex3f(z.x, z.y, z.z); 
            }
        }
        glEnd(); 
    }
}

//##############################################################################
//##############################################################################
void ModalViewer::
DrawImpulses()
{
    const REAL impulseScaling = 0.1 * _timeStepSize;
    const int N_frames = _rigidSoundObject->N_Impulses(); 
    _currentImpulseFrame = _currentFrame % N_frames; 
    // get impulse from object
    const REAL timeStart = CurrentTime(); 
    const REAL timeStop  = timeStart + _timeStepSize; 
    std::vector<ImpulseSeriesObject::ImpactRecord> impactRecords; 
    _rigidSoundObject->GetForces(timeStart, timeStop, impactRecords); 
    if (impactRecords.size() > 0) 
    {
        const int N_impacts = impactRecords.size(); 
        for (int imp_idx=0; imp_idx<N_impacts; ++imp_idx)
        {
            const auto &record = impactRecords.at(imp_idx); 
            const int &vertexID = record.appliedVertex; 
            const Vector3d &impulse = record.impactVector; 

            //_rigidSoundObject->GetImpulse(_currentImpulseFrame, timestamp, vertexID, impulse); 

            // draw impulse vector
            std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
            const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
            //const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
            //const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
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
    }
}

//##############################################################################
//##############################################################################
void ModalViewer::
animate()
{
    DrawOneFrameForward(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
keyPressEvent(QKeyEvent *e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers(); 
    bool optionsChanged = false;
    bool handled = true; 
    if ((e->key() == Qt::Key_I) && (modifiers == Qt::NoButton)) {
        _drawImpulse = !_drawImpulse; 
        optionsChanged = true;}
    if ((e->key() == Qt::Key_W) && (modifiers == Qt::NoButton)) {
        _wireframe = (_wireframe+1)%3; 
        optionsChanged = true;}
    else if ((e->key() == Qt::Key_BracketLeft) && (modifiers == Qt::NoButton)) {
        if (!animationIsStarted())
            DrawOneFrameBackward(); }
    else if ((e->key() == Qt::Key_BracketRight) && (modifiers == Qt::NoButton)) {
        if (!animationIsStarted())
            DrawOneFrameForward(); }
    else if ((e->key() == Qt::Key_M) && (modifiers == Qt::ShiftModifier)){
        _drawModes--; 
        UpdateVertexValues();
        optionsChanged = true;}
    else if ((e->key() == Qt::Key_M) && (modifiers == Qt::NoButton)){
        _drawModes++; 
        UpdateVertexValues();
        optionsChanged = true;}
    else if ((e->key() == Qt::Key_N) && (modifiers == Qt::NoButton)){
        _rigidSoundObject->AdvanceModalODESolvers(1);
        UpdateVertexValues();}
    else if ((e->key() == Qt::Key_R) && (modifiers == Qt::NoButton)){
        StepODEAndStoreResults();}
    else if ((e->key() == Qt::Key_F) && (modifiers == Qt::NoButton)){
        PrintAllFrequencies(std::cout);}
    else {
        handled = false;}


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
postSelection(const QPoint &point)
{
    _messageSelection = QString("Vertex ID= " + QString::number(selectedName()));
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
    _currentFrame = std::max<int>(_currentFrame, 0);
    updateGL();
    PrintFrameInfo(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
UpdateVertexValues()
{
    // color vertex using absolute modal values normalized
    if (_drawModes >= 0 && _drawModes < _rigidSoundObject->N_Modes())
    {
        // draw the modal displacement excited by impulse
        _rigidSoundObject->GetModalDisplacement(_drawModes, _vertexValues);

        //// draw the modes directly
        //_rigidSoundObject->GetVertexModeValuesNormalized(_drawModes, _vertexValues); 
        PrintFrameInfo(); 
    }
}

//##############################################################################
//##############################################################################
void ModalViewer::
PrepareImpulses()
{
    const std::string impulseFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/modalImpulses.txt"); 
    const std::string rigidsimConfigFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/default.cfg");
    ImpulseSeriesReader reader(impulseFile, rigidsimConfigFile); 
    std::shared_ptr<ImpulseSeriesObject> objectPtr = std::static_pointer_cast<ImpulseSeriesObject>(_rigidSoundObject); 
    reader.LoadImpulses(0, objectPtr); 
    _rigidSoundObject->GetRangeOfImpulses(_impulseRange.start, _impulseRange.stop); 
    _timeStepSize = _rigidSoundObject->GetRigidsimTimeStepSize();  // set the time step size always the same as rigidsim

    std::cout << "Impulses Read:\n"
              << " Number of impulses: " << _rigidSoundObject->N_Impulses() << "\n"
              << " Time range of impulses: [" << _impulseRange.start << ", " << _impulseRange.stop << "]\n"
              << "\n";
}

//##############################################################################
//##############################################################################
void ModalViewer::
PrepareModes()
{
    const std::string modeFile("/home/jui-hsien/code/acoustics/work/plate/plate.modes"); 
    const std::string parseFile("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml"); 
    ImpulseResponseParser parser(parseFile); 
    ModalMaterialList materials; 
    parser.GetModalMaterials(materials); 
    auto materialPtr = materials.at(0);

    _ODEStepSize = _timeStepSize / 40.0; 
    std::cout << "ODEStepSize = " << _ODEStepSize << std::endl;
    _rigidSoundObject->ModalAnalysisObject::Initialize(_ODEStepSize, modeFile, materialPtr); 
    _rigidSoundObject->InitializeModeVectors(); 
}

//##############################################################################
//##############################################################################
void ModalViewer::
RestoreDefaultDrawOptions()
{
    _drawImpulse = false; 
    _wireframe = 2;
    _displayMessage = true;
    _drawModes = 0; 
}

//##############################################################################
// Assume starts from time = 0
//##############################################################################
void ModalViewer::
StepODEAndStoreResults()
{
    assert(EQUAL_FLOATS(_rigidSoundObject->GetODESolverTime(), 0.0)); 
    REAL timeStop; 
    std::string outFile_displacement, outFile_q; 
    std::cout << "\nEnter time to stop (s): ";
    std::cin >> timeStop; 
    std::cout << "\nEnter file path to store vertex displacements: ";
    std::cin >> outFile_displacement; 
    std::cout << "\nEnter file path to store time-history of q: ";
    std::cin >> outFile_q; 
    const int N_steps = (timeStop - 0.0) / _ODEStepSize; 

    if (IO::ExistFile(outFile_displacement) || IO::ExistFile(outFile_q))
    {
        char yn = 'N'; 
        std::cout << "Either vertex displacement or q file exists. Overwrite both? [y/N] "; 
        std::cin >> yn; 
        if (yn == 'N')
            return; 
    } 

    std::ofstream of1(outFile_displacement.c_str()); 
    std::ofstream of2(outFile_q.c_str()); 
    std::cout << "\n\n" << N_steps << " steps are needed to advance ODE solvers. Store vertex displacements to " << outFile_displacement << "; store q to " << outFile_q << std::endl;
    _rigidSoundObject->AdvanceModalODESolvers(N_steps, of1, of2); 
    of1.close(); 
    of2.close(); 
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
              << " Draw Modes         : " << _drawModes << "\n"
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
    _message += "Current Time: " + QString::number(CurrentTime()) + "; "; 
    _message += "Current Impulse Frame: " + QString::number(_currentImpulseFrame) + "; ";
    _message += "Current Mode Frame: " + QString::number(_drawModes) + "(" + QString::number(_rigidSoundObject->GetModeFrequency(_drawModes)) + " Hz); "; 
    _message += "Current Modal ODE time: " + QString::number(_rigidSoundObject->GetODESolverTime()) + "; "; 
}

//##############################################################################
//##############################################################################
void ModalViewer::
PrintAllFrequencies(std::ostream &os)
{
    const int N_modes = _rigidSoundObject->N_Modes();
    for (int mode_idx=0; mode_idx<N_modes; ++mode_idx)
        os << "Mode " << mode_idx << ": " << _rigidSoundObject->GetModeFrequency(mode_idx) << " Hz" << std::endl;
}
