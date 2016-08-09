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
#include <iostream>

using namespace qglviewer;

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
SetAllKeyDescriptions()
{
    setKeyDescription(Qt::Key_W, "Toggle wireframe-only display"); 
    setKeyDescription(Qt::Key_B, "Toggle simulation box display"); 
    setKeyDescription(Qt::Key_R, "Run Simulator in the background of GL"); 
    setKeyDescription(Qt::Key_P, "Draw a sphere at given position"); 
    setKeyDescription(Qt::Key_C, "Clear all debug sphere"); 
    setKeyDescription(Qt::Key_N, "Draw an arrow"); 
    setKeyDescription(Qt::Key_N, "Draw arrows from file <x, y, z, nx, ny, nz>"); 
    setKeyDescription(Qt::Key_Y, "Draw slice for data display"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_C, "Clear all debug arrows"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_P, "Debug: draw failed reflections arrows"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_W, "Toggle slice grid lines"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_Y, "Toggle slice data pointer"); 
    setKeyDescription(Qt::ControlModifier + Qt::Key_P, "Draw arrows from file"); 
    setKeyDescription(Qt::AltModifier + Qt::Key_P, "Draw spheres from file"); 
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

    DrawLights();
}

//##############################################################################
// This method shouldn't be called directly. Call 'updateGL()' instead. See
// documentation for QGLViewer
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
draw()
{
    DrawMesh(); 
    DrawListeningPoints();
    DrawDebugCin();
    if (_sliceDataPointer == 0)
    {
        DrawSlices(0); 
    }
    else if (_sliceDataPointer == 1)
    {
        DrawSlices(1); 
    }
    else
    {
        DrawSlices(0); 
        DrawSlices(1); 
    }
    glColor3f(1.0, 1.0, 1.0);
    drawText(10, height()-20, _message); 

    glLineWidth(3.0f);
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

        // get transformations
        Eigen::Vector3d translation = object->GetTranslation(); 
        Eigen::Vector3d rotationAxis;
        REAL rotationAngle;
        object->GetRotationDegree(rotationAngle, rotationAxis); 
        glPushMatrix(); 
        glTranslated(translation[0], translation[1], translation[2]); 
        glRotated(rotationAngle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);

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

        glPopMatrix();
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
// Draw simulation listening point
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawListeningPoints()
{
    const auto &settings = _simulator->GetSolverSettings(); 
    const auto &points = settings->listeningPoints; 
    const int N_points = points.size(); 
    for (int pt_idx=0; pt_idx<N_points; ++pt_idx)
    {
        const Vector3d &vertex = points.at(pt_idx); 
        glPushMatrix();
        glTranslatef(vertex.x, vertex.y, vertex.z); 
        glColor3f(0.9f, 0.1f, 0.1f);
        GL_Wrapper::DrawSphere(5E-3, 10, 10); 
        glPopMatrix(); 
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawLights()
{
    glEnable(GL_LIGHT0);
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
    glEnable(GL_COLOR_MATERIAL);
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawSlices(const int &dataPointer)
{
    if (_sliceCin.size() == 0)
        return;

    const int N_slices = _sliceCin.size();
    std::vector<Eigen::MatrixXd> dataSlices(N_slices); 
    REAL maxCoeffSlices = std::numeric_limits<REAL>::min(); 
    REAL minCoeffSlices = std::numeric_limits<REAL>::max(); 
    // first fetch all data and compute min, max
    for (int s_idx=0; s_idx<N_slices; ++s_idx) 
    {
        const auto &slice = _sliceCin.at(s_idx); 
        Eigen::MatrixXd &data = dataSlices.at(s_idx); 
        if (dataPointer == 0)
            _simulator->GetSolver()->FetchPressureData(slice.samples, data);
        else if (dataPointer == 1)
            _simulator->GetSolver()->FetchPressureCellType(slice.samples, data);

        maxCoeffSlices = max(maxCoeffSlices, data.maxCoeff()); 
        minCoeffSlices = min(minCoeffSlices, data.minCoeff()); 
    }

    minCoeffSlices = (minCoeffSlices==maxCoeffSlices ?  maxCoeffSlices-EPS : minCoeffSlices);
    if (dataPointer == 0)
        _sliceColorMap->set_interpolation_range(minCoeffSlices, maxCoeffSlices); 
    else if (dataPointer == 1) 
        _sliceColorMap->set_interpolation_range(-1.0, 1.0); 

    // draw slices 
    for (int s_idx=0; s_idx<N_slices; ++s_idx) 
    {
        const auto &slice = _sliceCin.at(s_idx); 
        const Eigen::MatrixXd &data = dataSlices.at(s_idx); 
        if (_sliceWireframe == 0 || _sliceWireframe == 1)
        {
            glLineWidth(1.0f);
            glColor3f(0.4f, 0.4f, 0.4f); 
            const int N_gridLines = slice.gridLines.size()/2; 
            for (int l_idx=0; l_idx<N_gridLines; ++l_idx)
            {
                glBegin(GL_LINES);
                const Vector3d &vertex_0 = slice.gridLines.at(l_idx*2); 
                const Vector3d &vertex_1 = slice.gridLines.at(l_idx*2+1); 
                glVertex3f(vertex_0.x, vertex_0.y, vertex_0.z); 
                glVertex3f(vertex_1.x, vertex_1.y, vertex_1.z); 
                glEnd();
            }
        }

        if (_sliceWireframe == 0 || _sliceWireframe == 2)
        {
            glBegin(GL_QUADS);
            const int divisions = slice.N_sample_per_dim; 
            for (int dim_0_idx=0; dim_0_idx<divisions-1; ++dim_0_idx)
                for (int dim_1_idx=0; dim_1_idx<divisions-1; ++dim_1_idx)
                {
                    const int idx_0_0 =  dim_0_idx     *divisions + dim_1_idx; 
                    const int idx_0_1 =  dim_0_idx     *divisions + dim_1_idx + 1; 
                    const int idx_1_1 = (dim_0_idx + 1)*divisions + dim_1_idx; 
                    const int idx_1_0 = (dim_0_idx + 1)*divisions + dim_1_idx + 1; 
                    const Tuple3f c_0_0 = _sliceColorMap->get_interpolated_color(data(idx_0_0)); 
                    const Tuple3f c_0_1 = _sliceColorMap->get_interpolated_color(data(idx_0_1)); 
                    const Tuple3f c_1_0 = _sliceColorMap->get_interpolated_color(data(idx_1_0)); 
                    const Tuple3f c_1_1 = _sliceColorMap->get_interpolated_color(data(idx_1_1)); 
                    const Vector3d &vertex_0_0 = slice.samples.at(idx_0_0); 
                    const Vector3d &vertex_0_1 = slice.samples.at(idx_0_1); 
                    const Vector3d &vertex_1_0 = slice.samples.at(idx_1_0); 
                    const Vector3d &vertex_1_1 = slice.samples.at(idx_1_1); 
                    glColor3f(c_0_0.x, c_0_0.y, c_0_0.z); 
                    glVertex3f(vertex_0_0.x, vertex_0_0.y, vertex_0_0.z); 
                    glColor3f(c_0_1.x, c_0_1.y, c_0_1.z); 
                    glVertex3f(vertex_0_1.x, vertex_0_1.y, vertex_0_1.z); 
                    glColor3f(c_1_0.x, c_1_0.y, c_1_0.z); 
                    glVertex3f(vertex_1_0.x, vertex_1_0.y, vertex_1_0.z); 
                    glColor3f(c_1_1.x, c_1_1.y, c_1_1.z); 
                    glVertex3f(vertex_1_1.x, vertex_1_1.y, vertex_1_1.z); 
                }
            glEnd();
        }
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawDebugCin()
{
    glEnable(GL_LIGHTING);
    // debug sphere
    for (size_t sph_idx=0; sph_idx<_sphereCin.size(); ++sph_idx)
    {
        const auto &sphere = _sphereCin.at(sph_idx); 
        const REAL &x = sphere.origin.x; 
        const REAL &y = sphere.origin.y; 
        const REAL &z = sphere.origin.z; 
        glPushMatrix();
        glTranslatef(x, y, z); 
        glColor3f(0.0f, 1.0f, 0.0f); 
        GL_Wrapper::DrawSphere(5E-4 * sphere.scale, 10, 10);
        glPopMatrix();
    }
    glDisable(GL_LIGHTING);
   
    const REAL arrowScale = 1;
    // debug arrows
    for (size_t arr_idx=0; arr_idx<_arrowCin.size(); ++arr_idx)
    {
        glLineWidth(5.0f);
        const REAL &start_x = _arrowCin.at(arr_idx).start.x; 
        const REAL &start_y = _arrowCin.at(arr_idx).start.y; 
        const REAL &start_z = _arrowCin.at(arr_idx).start.z; 
        const REAL &stop_x = _arrowCin.at(arr_idx).start.x + _arrowCin.at(arr_idx).normal.x*arrowScale; 
        const REAL &stop_y = _arrowCin.at(arr_idx).start.y + _arrowCin.at(arr_idx).normal.y*arrowScale; 
        const REAL &stop_z = _arrowCin.at(arr_idx).start.z + _arrowCin.at(arr_idx).normal.z*arrowScale; 
        glBegin(GL_LINES);
        glColor3f(0.0f, 0.0f, 1.0f); 
        glVertex3f(start_x, start_y, start_z); 
        glVertex3f(stop_x, stop_y, stop_z); 
        glEnd();

        glPushMatrix();
        glTranslatef(start_x, start_y, start_z); 
        glColor3f(1.0f, 1.0f, 0.0f); 
        GL_Wrapper::DrawSphere(5E-4, 10, 10);
        glPopMatrix();
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
animate()
{
    DrawOneFrameForward();
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
        _wireframe = (_wireframe+1)%4; 
        optionsChanged = true;
    }
    if ((e->key() == Qt::Key_G) && (modifiers == Qt::ShiftModifier)) {
        _sliceWireframe = (_sliceWireframe+1)%4; 
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_B) && (modifiers == Qt::NoButton)) {
        _drawBox = !_drawBox; 
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_BracketRight) && (modifiers == Qt::NoButton)) {
        if (!animationIsStarted())
            DrawOneFrameForward(); 
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::NoButton)) {
            Sphere sphere;
            std::cout << "Sphere <x, y, z, scale>: " << std::flush; 
            std::cin >> sphere.origin.x >> sphere.origin.y >> sphere.origin.z >> sphere.scale; 
            _sphereCin.push_back(sphere); 
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::ShiftModifier)) {
        // write to disk 
        std::string filename; 
        std::cout << "Debug: failed reflection arrow filename: " << std::flush; 
        std::cin >> filename; 
        auto &sceneObjects = _simulator->GetSceneObjects(); 
        const int N = sceneObjects->N();
        sceneObjects->WriteFailedReflections(filename);

        // push all arrows to debug draw
        for (int obj_idx=0; obj_idx<N; ++obj_idx)
        {
            const std::string objFilename = filename + sceneObjects->GetMeshName(obj_idx); 
            std::ifstream inFile(objFilename.c_str()); 
            Vector3f x, n; 
            while(inFile >> x.x >> x.y >> x.z >> n.x >> n.y >> n.z)
            {
                Arrow arrow; 
                arrow.start = x; 
                arrow.normal = n; 
                _arrowCin.push_back(arrow); 
            }
        }
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::ControlModifier)) {
        std::string filename; 
        std::cout << "To-draw arrow filename: " << std::flush; 
        std::cin >> filename; 
        std::ifstream inFile(filename.c_str()); 
        Vector3f x, n; 
        while(inFile >> x.x >> x.y >> x.z >> n.x >> n.y >> n.z)
        {
            Arrow arrow; 
            arrow.start = x; 
            arrow.normal = n; 
            _arrowCin.push_back(arrow); 
        }
        std::cout << " " << _arrowCin.size() << " arrows read and ready to draw.\n" << std::flush;
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::AltModifier)) {
        std::string filename; 
        std::cout << "To-draw sphere filename: " << std::flush; 
        std::cin >> filename; 
        REAL scale; 
        std::cout << "To-draw sphere scale: " << std::flush; 
        std::cin >> scale; 
        std::ifstream inFile(filename.c_str()); 
        Sphere sphere;
        Vector3f &x = sphere.origin;
        while(inFile >> x.x >> x.y >> x.z)
        {
            sphere.scale = scale; 
            _sphereCin.push_back(sphere); 
        }
        std::cout << " " << _sphereCin.size() << " spheres read and ready to draw.\n" << std::flush;
    }
    else if ((e->key() == Qt::Key_N) && (modifiers == Qt::NoButton)) {
            Vector3f x, n; 
            std::cout << "Arrow start location <x, y, z>: " << std::flush; 
            std::cin >> x.x >> x.y >> x.z; 
            std::cout << "Arrow normal <x, y, z>: " << std::flush; 
            std::cin >> n.x >> n.y >> n.z; 
            Arrow arrow; 
            arrow.start = x; 
            arrow.normal = n; 
            _arrowCin.push_back(arrow); 
    }
    else if ((e->key() == Qt::Key_M) && (modifiers == Qt::NoButton)) {
        std::string filename; 
        std::cout << "Read debug arrow files: " << std::flush; 
        std::cin >> filename; 
        std::ifstream inFile(filename.c_str()); 
        if (inFile)
        {
            Vector3f x, n; 
            while (inFile >> x.x >> x.y >> x.z >> n.x >> n.y >> n.z)
            {
                Arrow arrow; 
                arrow.start = x; 
                arrow.normal = n; 
                _arrowCin.push_back(arrow); 
            }
        }
    }
    else if ((e->key() == Qt::Key_Y) && (modifiers == Qt::NoButton)) {
        Slice slice; 
        std::cout << "Slice dim: " << std::flush; 
        std::cin >> slice.dim;
        std::cout << "Slice origin: " << std::flush; 
        std::cin >> slice.origin.x >> slice.origin.y >> slice.origin.z; 
        ConstructSliceSamples(slice);
        _sliceCin.push_back(slice); 
    }
    else if ((e->key() == Qt::Key_Y) && (modifiers == Qt::ShiftModifier)) {
        _sliceDataPointer = (_sliceDataPointer + 1)%2; 
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_C) && (modifiers == Qt::NoButton)) {
            _sphereCin.clear(); 
    }
    else if ((e->key() == Qt::Key_C) && (modifiers == Qt::ShiftModifier)) {
            _arrowCin.clear(); 
    }
    else if ((e->key() == Qt::Key_R) && (modifiers == Qt::NoButton)) {
        DrawOneFrameForward(); 
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
ConstructSliceSamples(Slice &slice)
{
    const int &dim = slice.dim; 
    const Vector3d &origin = slice.origin; 
    Vector3Array &samples = slice.samples; 
    Vector3Array &gridLines = slice.gridLines; 

    const auto &settings = _simulator->GetSolverSettings(); 
    const REAL cellSize = settings->cellSize; 
    const int division = settings->cellDivisions; 
    const REAL halfLength = (REAL)division*cellSize / 2.0; 
    const REAL minBound = -halfLength + cellSize/2.0; 

    const int dim_0 = (dim + 1) % 3; 
    const int dim_1 = (dim + 2) % 3; 
    for (int dim_0_idx=0; dim_0_idx<division; ++dim_0_idx)
        for (int dim_1_idx=0; dim_1_idx<division; ++dim_1_idx)
        {
            Vector3d sample; 
            sample(dim) = origin(dim); 
            sample(dim_0) = minBound + cellSize*(REAL)dim_0_idx; 
            sample(dim_1) = minBound + cellSize*(REAL)dim_1_idx; 
            samples.push_back(sample); 
        }
    slice.N_sample_per_dim = division; 
    slice.minBound = minBound; 
    slice.maxBound = minBound + (REAL)(division-1)*cellSize; 


    // horizontal grid lines
    const REAL xStart = slice.minBound - cellSize/2.0; 
    const REAL xStop =  slice.maxBound + cellSize/2.0; 
    for (int dim_0_idx=0; dim_0_idx<division+1; ++dim_0_idx)
    {
        Vector3d start; 
        Vector3d stop; 
        start(dim) = origin(dim); 
        stop(dim) = origin(dim); 
        start(dim_0) = xStart + cellSize*(REAL)dim_0_idx;
        stop(dim_0) = xStart + cellSize*(REAL)dim_0_idx;
        start(dim_1) = xStart; 
        stop(dim_1) = xStop; 
        gridLines.push_back(start); 
        gridLines.push_back(stop); 
    }

    // vertical grid lines
    const REAL yStart = slice.minBound - cellSize/2.0; 
    const REAL yStop =  slice.maxBound + cellSize/2.0; 
    for (int dim_1_idx=0; dim_1_idx<division+1; ++dim_1_idx)
    {
        Vector3d start; 
        Vector3d stop; 
        start(dim) = origin(dim); 
        stop(dim) = origin(dim); 
        start(dim_0) = yStart; 
        stop(dim_0) = yStop; 
        start(dim_1) = yStart + cellSize*(REAL)dim_1_idx;
        stop(dim_1) = yStart + cellSize*(REAL)dim_1_idx;
        gridLines.push_back(start); 
        gridLines.push_back(stop); 
    }

    // make colormap for this slice.
    if (!_sliceColorMap)
        _sliceColorMap = std::make_shared<JetColorMap>(); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawOneFrameForward()
{
    _currentFrame++;
    PrintFrameInfo();
    //_simulator->TestAnimateObjects(150); 
    _simulator->RunForSteps(1); 
    updateGL(); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
RestoreDefaultDrawOptions()
{
    _wireframe = 2;
    _sliceWireframe = 2;
    _drawBox = true; 
    _sliceDataPointer = 0; 
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
              << " Draw slice wireframe only: " << _sliceWireframe << "\n"
              << " Draw slice data pointer: " << _sliceDataPointer << "\n"
              << "\n"; 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
PrintFrameInfo()
{
    //const std::string frameInfo("Current Frame: " + std::to_string(_currentFrame)); 
    _message = QString("");
    _message += "Current Frame: " + QString::number(_currentFrame) + "; "; 
    _message += "max= " + QString::number(_drawAbsMax) + "; "; 
}
