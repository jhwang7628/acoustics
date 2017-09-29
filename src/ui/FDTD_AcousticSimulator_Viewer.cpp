#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMap>
#include <QCursor>
#include <qimage.h>
#include <complex>
#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/FDTD_Objects.h>
#include <modal_model/KirchhoffIntegralSolver.h> 
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
    setKeyDescription(Qt::Key_T, "Toggle perspective/orthogonal view"); 
    setKeyDescription(Qt::Key_G, "Draw ground"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_H, "Draw hashed cells"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_S, "Save slice data to file"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_C, "Clear all debug arrows"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_F, "Debug: execute some debug function"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_P, "Debug: draw debug arrows stored in FDTD_RigidObject class"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_W, "Toggle slice grid lines"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_D, "Change slice division (default: 80)"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_R, "Run half-step"); 
    setKeyDescription(Qt::ControlModifier + Qt::Key_C, "Clear all slices"); 
    setKeyDescription(Qt::ControlModifier + Qt::Key_P, "Draw arrows from file"); 
    setKeyDescription(Qt::ControlModifier + Qt::Key_F, "Toggle fixed colormap"); 
    setKeyDescription(Qt::ControlModifier + Qt::Key_R, "Read FBem solutions to BEM solver"); 
    setKeyDescription(Qt::ShiftModifier + Qt::Key_Y, "Toggle slice data pointer forward"); 
    setKeyDescription(Qt::ControlModifier + Qt::Key_Y, "Toggle slice data pointer backward"); 
    setKeyDescription(Qt::AltModifier + Qt::Key_P, "Draw spheres from file"); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
init()
{
    restoreStateFromFile();
    RestoreDefaultDrawOptions(); 
    setSceneRadius(5.0);
    glDisable(GL_LIGHTING);
    glPointSize(3.0);
    //setGridIsDrawn();
    setBackgroundColor(QColor(102,153,255)); 
    SetAllKeyDescriptions();

    setAnimationPeriod(40); // in milliseconds
    init_gl();

    // determine if this program is run remotely
    char *vDISPLAY = getenv("DISPLAY"); 
    if (strcmp(vDISPLAY, ":0") != 0) 
    {
        _remoteConnection = true;
        std::cout << "Remote connection detected.\n";
        RestoreDefaultDrawOptions();
    }

    std::cout << "\n>> Press key 'h' for help, 'esc' for exit.\n\n";
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
init_gl()
{
    // fetch objects to set up for colors, material for gl rendering
    const int N_objects = _simWorld->GetSceneObjects()->N(); 
    _objectColors.resize(N_objects); 
    for (int obj_idx=0; obj_idx<N_objects; ++obj_idx)
    {
        const float x = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) *0.6f + 0.4f;
        const float y = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) *0.6f + 0.4f;
        const float z = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) *0.6f + 0.4f;
        _objectColors[obj_idx] = Vector3f(x, y, z);
    }

    // antialiasing
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Enable GL textures
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    LoadGround(); 

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
    DrawImpulses();
    DrawSelection(); 
    DrawDebugCin();
    if (_sliceWireframe.count() != 0)
        DrawSlices(_sliceDataPointer); 
    DrawGround();
    if (_drawHashedCells)
        DrawHashedCells();

    // for some reason drawText causes bug if run remotely
    if (!_remoteConnection)
    {
        glColor3f(1.0, 1.0, 1.0);
        drawText(10, height()-20, _message); 
        drawText(10, height()-40, _messageColormap); 
    }

    glLineWidth(3.0f);
    if (_drawBoxLis)
    {
        DrawListeningPoints();
        DrawBox(); 
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
drawWithNames()
{
    const REAL ballSize = _solverSettings->cellSize/1.9; 
    // draw cell centroid near the slices 
    for (auto &slice : _sliceCin)
    {
        const REAL offset = slice.origin[slice.dim]; 
        slice.intersectingUnit->simulator->GetSolver()->SampleAxisAlignedSlice(slice.dim, offset, slice.cells); 
        for (const auto &cell : slice.cells) 
        {
            const Vector3d &vertex = cell.centroidPosition; 
            glPushMatrix(); 
            glTranslatef(vertex.x, vertex.y, vertex.z); 
            glPushName(cell.index); 
            GL_Wrapper::DrawSphere(ballSize, 3, 3);
            glPopName(); 
            glPopMatrix(); 
        }
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawMesh()
{
    const auto &sceneObjects = _simWorld->GetSceneObjects(); 
    const auto &rigidObjects = sceneObjects->GetRigidObjects();
    const int N_objects = rigidObjects.size(); 
    for (int obj_idx=0; obj_idx<N_objects; ++obj_idx)
    {
        // draw rigid sound object mesh
        auto &object = rigidObjects.at(obj_idx); 
        std::shared_ptr<TriangleMesh<REAL> > meshPtr = object->GetMeshPtr();
        const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
        const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
        const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
        const int N_triangles = triangles.size(); 
        const int N_vertices = vertices.size(); 
        const REAL offsetEpsilon = 1E-5;

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
                x = object->ObjectToWorldPoint(x); 
                y = object->ObjectToWorldPoint(y); 
                z = object->ObjectToWorldPoint(z); 
                glVertex3f(x.x, x.y, x.z); 
                glVertex3f(y.x, y.y, y.z); 

                glVertex3f(y.x, y.y, y.z); 
                glVertex3f(z.x, z.y, z.z); 

                glVertex3f(z.x, z.y, z.z); 
                glVertex3f(x.x, x.y, x.z); 
            }
            glEnd(); 
        }

        // get curvatures 
        const std::vector<REAL> *meanCurvatures = (_meshDataPointer==1 ? meshPtr->mean_curvatures() : nullptr);
        std::shared_ptr<ColorMap> curvatureColorMap; 
        if (meanCurvatures) 
        {
            curvatureColorMap = std::make_shared<JetColorMap>();
            REAL maxCurvature = std::numeric_limits<REAL>::min(); 
            REAL minCurvature = std::numeric_limits<REAL>::max(); 
            for (int v_idx=0; v_idx<N_vertices; ++v_idx)
            {
                maxCurvature = std::max<REAL>(maxCurvature, meanCurvatures->at(v_idx)); 
                minCurvature = std::min<REAL>(minCurvature, meanCurvatures->at(v_idx)); 
            }
            curvatureColorMap->set_interpolation_range(minCurvature, maxCurvature); 
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
                Point3<REAL> x = vertices.at(triangle.x); 
                Point3<REAL> y = vertices.at(triangle.y); 
                Point3<REAL> z = vertices.at(triangle.z); 
                Vector3<REAL> nx = normals.at(triangle.x); 
                Vector3<REAL> ny = normals.at(triangle.y); 
                Vector3<REAL> nz = normals.at(triangle.z); 
                x = object->ObjectToWorldPoint(x); 
                y = object->ObjectToWorldPoint(y); 
                z = object->ObjectToWorldPoint(z); 
                nx = object->ObjectToWorldVector(nx); 
                ny = object->ObjectToWorldVector(ny); 
                nz = object->ObjectToWorldVector(nz); 
                if (_meshDataPointer == 0)
                {
                    glNormal3f(nx.x, nx.y, nx.z); 
                    glVertex3f(x.x, x.y, x.z); 
                    glNormal3f(ny.x, ny.y, ny.z); 
                    glVertex3f(y.x, y.y, y.z); 
                    glNormal3f(nz.x, nz.y, nz.z); 
                    glVertex3f(z.x, z.y, z.z); 
                }
                else
                {
                    const REAL xCurvature = meanCurvatures->at(triangle.x);
                    const REAL yCurvature = meanCurvatures->at(triangle.y);
                    const REAL zCurvature = meanCurvatures->at(triangle.z);
                    const Tuple3f cx = curvatureColorMap->get_interpolated_color(xCurvature); 
                    const Tuple3f cy = curvatureColorMap->get_interpolated_color(yCurvature); 
                    const Tuple3f cz = curvatureColorMap->get_interpolated_color(zCurvature); 
                    glNormal3f(nx.x, nx.y, nx.z); 
                    glColor3f(cx.x, cx.y, cx.z);
                    glVertex3f(x.x, x.y, x.z); 
                    glNormal3f(ny.x, ny.y, ny.z); 
                    glColor3f(cy.x, cy.y, cy.z);
                    glVertex3f(y.x, y.y, y.z); 
                    glNormal3f(nz.x, nz.y, nz.z); 
                    glColor3f(cz.x, cz.y, cz.z);
                    glVertex3f(z.x, z.y, z.z); 
                }
            }
            glEnd(); 
        }

        // draw rasterized cells
        if (_wireframe == 3)
        {
            const auto &simUnits = _simWorld->GetActiveSimUnits();
            for (const auto &unit : simUnits)
            {
                auto &grid = unit->simulator->GetGrid(); 
                const auto &gcMap = grid.GetGhostCells(); 
                glBegin(GL_QUADS); 
                for (const auto &m : gcMap)
                {
                    const int &cell_idx = m.second->ownerCell; 
                    const int &nmldir = abs(m.second->topology) - 1; // 0:x; 1:y; 2:z
                    const Vector3d position = (grid.pressureFieldPosition(m.second->ownerCell)
                                              +grid.pressureFieldPosition(m.second->neighbourCell))/2.0; 
                    Vector3d lowerCorner = position;
                    Vector3d upperCorner = position; 
                    lowerCorner[(nmldir+1)%3] -= 0.5*_solverSettings->cellSize; 
                    lowerCorner[(nmldir+2)%3] -= 0.5*_solverSettings->cellSize; 
                    upperCorner[(nmldir+1)%3] += 0.5*_solverSettings->cellSize; 
                    upperCorner[(nmldir+2)%3] += 0.5*_solverSettings->cellSize; 
                    Vector3d nml(0,0,0);
                    nml[nmldir] = (m.second->topology > 0 ? 1.0 : -1.0); 
                    Tuple3f color(1.,0.,0.);
                    if (_sliceColorMap)
                        color = _sliceColorMap->get_interpolated_color(m.second->pressure);
                    glNormal3f(nml[0],nml[1],nml[2]);
                    glColor3f(color[0],color[1],color[2]); 
                    Vector3d lu = lowerCorner; lu[(nmldir+1)%3] = upperCorner[(nmldir+1)%3];
                    Vector3d ul = lowerCorner; ul[(nmldir+2)%3] = upperCorner[(nmldir+2)%3];
                    glVertex3f(lowerCorner[0],lowerCorner[1],lowerCorner[2]); 
                    glVertex3f(lu[0],lu[1],lu[2]); 
                    glVertex3f(upperCorner[0],upperCorner[1],upperCorner[2]); 
                    glVertex3f(ul[0],ul[1],ul[2]); 
                    //const auto &cell_idx = m.second->ownerCell; 
                    //MAC_Grid::Cell cell; 
                    //unit->simulator->GetSolver()->FetchCell(cell_idx, cell); 
                    //glColor3f(1.0f, 0.0f, 0.0f); 
                    //GL_Wrapper::DrawBox(&(cell.lowerCorner.x), &(cell.upperCorner.x)); 
                }
                glEnd();
                const auto &bgcMap = unit->simulator->GetGrid().GetBoundaryGhostCells(); 
                for (const auto &m : bgcMap)
                {
                    const auto &cell_idx = m.second->ownerCell; 
                    MAC_Grid::Cell cell; 
                    unit->simulator->GetSolver()->FetchCell(cell_idx, cell); 
                    glColor3f(1.0f, 0.0f, 0.0f); 
                    GL_Wrapper::DrawBox(&(cell.lowerCorner.x), &(cell.upperCorner.x)); 
                }
            }
        }
        glDisable(GL_LIGHTING);

        // draw bounding box only 
        if (_wireframe == 4)
        {
            // this draws bounding box in world space
            const REAL maxValue = std::numeric_limits<REAL>::max(); 
            const REAL minValue = std::numeric_limits<REAL>::lowest(); 
            Point3d low(maxValue, maxValue, maxValue); 
            Point3d  up(minValue, minValue, minValue); 
            for (const auto &v : vertices)
            {
                const Vector3<REAL> vv = object->ObjectToWorldPoint(v); 
                low.x = min(low.x, vv.x); 
                low.y = min(low.y, vv.y); 
                low.z = min(low.z, vv.z); 
                 up.x = max( up.x, vv.x); 
                 up.y = max( up.y, vv.y); 
                 up.z = max( up.z, vv.z); 
            }
            const auto &color = _objectColors.at(obj_idx); 
            glColor3f(color.x, color.y, color.z); 
            GL_Wrapper::DrawWireBox(&(low.x), &(up.x)); 
        }
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawImpulses()
{
    const REAL time = _simWorld->GetWorldTime();
    const auto &sceneObjects = _simWorld->GetSceneObjects(); 
    for (int obj_idx=0; obj_idx<sceneObjects->N(); ++obj_idx)
    {
        const auto &object = sceneObjects->GetPtr(obj_idx); 
        std::shared_ptr<TriangleMesh<REAL> > meshPtr = object->GetMeshPtr();
        const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
        std::vector<ImpulseSeriesObject::ImpactRecord> records; 
        std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object)
            ->GetImpulseWithinSupport(time, records); 
        for (const auto &imp : records) 
        {
            if (imp.supportLength < SMALL_NUM)
                continue;
            const auto &color = _objectColors.at(obj_idx); 
            Point3<REAL> vertexEnd   = vertices.at(imp.appliedVertex); 
            Point3<REAL> vertexBegin = vertexEnd - imp.impactVector * _drawImpulseScaling;
            vertexEnd = object->ObjectToWorldPoint(vertexEnd); 
            vertexBegin = object->ObjectToWorldPoint(vertexBegin); 
            glLineWidth(3.0f); 
            glBegin(GL_LINES); 
            glColor3f(color.x, color.y, color.z); 
            glVertex3f(vertexBegin.x, vertexBegin.y, vertexBegin.z); 
            glVertex3f(vertexEnd.x, vertexEnd.y, vertexEnd.z); 
            glEnd(); 
            glPushMatrix();
            glTranslatef(vertexEnd.x, vertexEnd.y, vertexEnd.z); 
            GL_Wrapper::DrawSphere(0.5E-3, 10, 10); 
            glPopMatrix(); 
        }
    }
}

//##############################################################################
// Draw simulation box
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawBox()
{
    const auto bboxs = _simWorld->GetSolverBBoxs(); 
    for (const auto &p : bboxs)
    {
        auto &unit = p.first; 
        auto newBoxCenter = SimWorld::rasterizer.rasterize(
                unit->unitCenter);  
        auto oldBoxCenter = SimWorld::rasterizer.rasterize(
                unit->BoundingBoxCenter());
        if (!newBoxCenter.equals(oldBoxCenter))
        {
            unit->boxCenterChanged = true;
        }
        const Vector3d &minBound = p.second.minBound(); 
        const Vector3d &maxBound = p.second.maxBound(); 
        glColor3f(1.0f, 1.0f, 1.0f);
        GL_Wrapper::DrawWireBox(&minBound[0], &maxBound[0]); 
        if (_sceneBox) 
        {
            glColor3f(1.0f, 0.0f, 1.0f);
            GL_Wrapper::DrawWireBox(&(_sceneBox->minBound()[0]), 
                                    &(_sceneBox->maxBound()[0])); 
        }

        // draw field center
        glEnable(GL_LIGHTING);
        glPushMatrix();
        const auto &center = unit->BoundingBoxCenter();
        glTranslatef(center.x, 
                     center.y, 
                     center.z); 
        glColor3f(0.9f, 0.9f, 0.1f);
        GL_Wrapper::DrawSphere(8E-4, 30, 30); 
        glPopMatrix(); 
        glDisable(GL_LIGHTING);
    }
}

//##############################################################################
// Load image for ground
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
LoadGround()
{
    QImage glImg = QGLWidget::convertToGLFormat(QImage(textureMeta.name));
    glTexImage2D(GL_TEXTURE_2D, 0, 4, glImg.width(), glImg.height(), 0,
            GL_RGBA, GL_UNSIGNED_BYTE, glImg.bits());
    textureMeta.groundLoaded = true; 
}

//##############################################################################
// Draw simulation listening point
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawGround()
{
    if (_drawGround == 0) 
    {
        return; 
    }
    else if (_drawGround == 1)
    {
        const float GD_SIZE = 0.01;
        const float step = GD_SIZE * 10;
        float d = step;
        glColor3f(0.7, 0.7, 0.7);
        for(int i=0; i<20; ++i, d+=step)
        {
            glBegin(GL_LINE_LOOP);
            glVertex3f(-d, 0, -d);
            glVertex3f( d, 0, -d);
            glVertex3f( d, 0,  d);
            glVertex3f(-d, 0,  d);
            glEnd();
        }
    }
    else if (textureMeta.groundLoaded && _drawGround == 2)
    {
        const float floorsize = 2.0f; 
        const float uvbound = floorsize/textureMeta.blockSizeUV; 
        glEnable(GL_TEXTURE_2D);
        // Display the quad
        glColor3f(1.0, 1.0, 1.0);
        glNormal3f(0.0, 1.0, 0.0);
        glBegin(GL_QUADS);
        glTexCoord2f(-uvbound, -uvbound);	glVertex3f(-floorsize, 0., -floorsize);
        glTexCoord2f(-uvbound,  uvbound);	glVertex3f(-floorsize, 0.,  floorsize);
        glTexCoord2f( uvbound,  uvbound);	glVertex3f( floorsize, 0.,  floorsize);
        glTexCoord2f( uvbound, -uvbound);	glVertex3f( floorsize, 0., -floorsize);
        glEnd();
        glDisable(GL_TEXTURE_2D);
    } 
}

//##############################################################################
// Draw simulation listening point
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawListeningPoints()
{
    const auto &settings = _solverSettings; 
    const auto &points = settings->listeningPoints; 
    const int N_points = points.size(); 
    glEnable(GL_LIGHTING);
    for (int pt_idx=0; pt_idx<N_points; ++pt_idx)
    {
        const Vector3d &vertex = points.at(pt_idx); 
        glPushMatrix();
        glTranslatef(vertex.x, vertex.y, vertex.z); 
        glColor3f(0.9f, 0.1f, 0.1f);
        GL_Wrapper::DrawSphere(1E-2, 30, 30); 
        glPopMatrix(); 
    }
    glDisable(GL_LIGHTING);
    auto units = _simWorld->GetActiveSimUnits(); 
    glPointSize(1.0); 
    glColor3f(0.9f, 0.1f, 0.9f);
    glBegin(GL_POINTS);
    for (const auto &unit : units)
    {
        const auto &spks = unit->listen->speakers; 
        for (const auto &spk : spks)
            glVertex3f(spk.x,spk.y,spk.z); 
    }
    glEnd();
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawSelection()
{
    if (selectedName() != -1 && _listenedUnit)
    {
        // draw box
        const int cell_idx = selectedName(); 
        MAC_Grid::Cell cell; 
        auto &simulator = _listenedUnit->simulator;
        simulator->GetSolver()->FetchCell(cell_idx, cell); 
        glColor3f(1.0f, 1.0f, 1.0f); 
        GL_Wrapper::DrawWireBox(&(cell.lowerCorner.x), &(cell.upperCorner.x)); 

        // draw hashed triangles
        const auto search = simulator->GetSolver()->GetGrid().GetFVMetaData().cellMap.find(cell_idx); 
        if (search != simulator->GetSolver()->GetGrid().GetFVMetaData().cellMap.end())
        {
            for (std::vector<MAC_Grid::TriangleIdentifier>::const_iterator it=(search->second)->begin(); it!=(search->second)->end(); ++it)
            {
                const auto &object = simulator->GetSceneObjects()->GetPtr(it->objectID);
                const std::vector<Point3<REAL> > &vertices = object->GetMeshPtr()->vertices(); 
                const Tuple3ui &triangle = object->GetMeshPtr()->triangle_ids(it->triangleID);
                const Vector3f &color = _objectColors.at(it->objectID); 
                Point3<REAL> x = vertices.at(triangle.x); 
                Point3<REAL> y = vertices.at(triangle.y); 
                Point3<REAL> z = vertices.at(triangle.z); 
                x = object->ObjectToWorldPoint(x); 
                y = object->ObjectToWorldPoint(y); 
                z = object->ObjectToWorldPoint(z); 
                glColor3f(color.x, color.y, color.z); 
                glBegin(GL_TRIANGLES); 
                glVertex3f(x.x, x.y, x.z); 
                glVertex3f(y.x, y.y, y.z); 
                glVertex3f(z.x, z.y, z.z); 
                glEnd(); 
                glPushMatrix(); 
                glTranslatef(it->centroid.x, it->centroid.y, it->centroid.z); 
                glColor3f(0.2f, 0.2f, 0.2f); 
                GL_Wrapper::DrawSphere(3E-4, 6, 6);
                glPopMatrix(); 
            }
        }

#ifdef USE_FV
        // draw ghost cell related stuff
        const auto gc = simulator->GetSolver()->GetGrid().GetGhostCell(cell_idx);
        typedef std::vector<MAC_Grid::GhostCell::VolumeSamples>::iterator Iterator_VS;
        typedef std::vector<MAC_Grid::GhostCell::BoundarySamples>::iterator Iterator_BS;
        if (gc)
        {
            //for (Iterator_BS sp=gc->boundarySamples.begin(); sp!=gc->boundarySamples.end(); ++sp)
            for (Iterator_VS sp=gc->volumeSamples.begin(); sp!=gc->volumeSamples.end(); ++sp)
            {
                const Vector3d &pos = sp->position; 
                glPushMatrix(); 
                glTranslatef(pos.x, pos.y, pos.z); 
                if (sp->isBulk)
                    glColor3f(0.0f, 1.0f, 0.0f); 
                else
                    glColor3f(0.2f, 0.2f, 0.2f); 
                GL_Wrapper::DrawSphere(1E-4, 5, 5);
                glPopMatrix(); 
            }
        }
#endif
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
    // draw slices 
    ComputeAndCacheSliceData(dataPointer); // do-nothing if cached
    for (int s_idx=0; s_idx<N_slices; ++s_idx) 
    {
        auto &slice = _sliceCin.at(s_idx); 
        const Eigen::MatrixXd &data = slice.data; 
        if (_sliceWireframe[0])
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

        if (_sliceWireframe[1] && data.rows()!=0 && data.cols()!=0)
        {
            glBegin(GL_QUADS);
            const Tuple3i divisions = slice.N_sample_per_dim; 
            const int dim_0 = (slice.dim + 1) %3; 
            const int dim_1 = (slice.dim + 2) %3; 

            for (int dim_0_idx=0; dim_0_idx<divisions[dim_0]-1; ++dim_0_idx)
                for (int dim_1_idx=0; dim_1_idx<divisions[dim_1]-1; ++dim_1_idx)
                {
                    const int idx_0_0 =  dim_0_idx     *divisions[dim_1] + dim_1_idx; 
                    const int idx_0_1 =  dim_0_idx     *divisions[dim_1] + dim_1_idx + 1; 
                    const int idx_1_1 = (dim_0_idx + 1)*divisions[dim_1] + dim_1_idx; 
                    const int idx_1_0 = (dim_0_idx + 1)*divisions[dim_1] + dim_1_idx + 1; 
                    const Tuple3f c_0_0 = _sliceColorMap->get_interpolated_color(data(idx_0_0, 0)); 
                    const Tuple3f c_0_1 = _sliceColorMap->get_interpolated_color(data(idx_0_1, 0)); 
                    const Tuple3f c_1_0 = _sliceColorMap->get_interpolated_color(data(idx_1_0, 0)); 
                    const Tuple3f c_1_1 = _sliceColorMap->get_interpolated_color(data(idx_1_1, 0)); 
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
    glBegin(GL_LINES);
    for (size_t arr_idx=0; arr_idx<_arrowCin.size(); ++arr_idx)
    {
        glLineWidth(5.0f);
        const REAL &start_x = _arrowCin.at(arr_idx).start.x; 
        const REAL &start_y = _arrowCin.at(arr_idx).start.y; 
        const REAL &start_z = _arrowCin.at(arr_idx).start.z; 
        const REAL &stop_x = _arrowCin.at(arr_idx).start.x + _arrowCin.at(arr_idx).normal.x*arrowScale; 
        const REAL &stop_y = _arrowCin.at(arr_idx).start.y + _arrowCin.at(arr_idx).normal.y*arrowScale; 
        const REAL &stop_z = _arrowCin.at(arr_idx).start.z + _arrowCin.at(arr_idx).normal.z*arrowScale; 
        glColor3f(1.0f, 1.0f, 0.0f); 
        glVertex3f(start_x, start_y, start_z); 
        glVertex3f(stop_x, stop_y, stop_z); 

        // end spheres
        //glEnable(GL_LIGHTING);
        //glPushMatrix();
        //glTranslatef(start_x, start_y, start_z); 
        //glColor3f(1.0f, 1.0f, 0.0f); 
        //GL_Wrapper::DrawSphere(2.5E-4, 10, 10);
        //glPopMatrix();
        //glDisable(GL_LIGHTING);
    }
    glEnd();
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawHashedCells()
{
    // FIXME debug
    //typedef std::unordered_map<int, std::shared_ptr<std::vector<MAC_Grid::TriangleIdentifier> > >::const_iterator Iterator; 
    //const MAC_Grid::FVMetaData &fvMetaData = _simulator->GetSolver()->GetGrid().GetFVMetaData();
    //const auto &cellMap = fvMetaData.cellMap; 
    //for (Iterator it=cellMap.begin(); it!=cellMap.end(); ++it)
    //{
    //    MAC_Grid::Cell cell; 
    //    _simulator->GetSolver()->FetchCell(it->first, cell); 
    //    const Vector3f &color = _objectColors.at(it->second->at(0).objectID);
    //    glColor3f(color.x, color.y, color.z);
    //    GL_Wrapper::DrawWireBox(&(cell.lowerCorner.x), &(cell.upperCorner.x)); 
    //}
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
        _wireframe = (_wireframe+1)%6; 
        optionsChanged = true;
    }
    if ((e->key() == Qt::Key_W) && (modifiers == Qt::ControlModifier)) {
        _meshDataPointer = (_meshDataPointer+1)%2; 
        optionsChanged = true;
    }
    if ((e->key() == Qt::Key_D) && (modifiers == Qt::ShiftModifier)) {
        std::cout << "Input slice division: " << std::flush; 
        PRE_CIN_CLEAR;
        std::cin >> _sliceDivision; 
        POST_CIN_CLEAR;
    }
    if ((e->key() == Qt::Key_H) && (modifiers == Qt::ShiftModifier)) {
        _drawHashedCells = !_drawHashedCells;
        optionsChanged = true;
    }
    if ((e->key() == Qt::Key_W) && (modifiers == Qt::ShiftModifier)) {
        if (_sliceWireframe.count() == _sliceWireframe.size())
            _sliceWireframe.reset(); 
        else
            _sliceWireframe = std::bitset<2>(_sliceWireframe.to_ulong()+1);
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_B) && (modifiers == Qt::NoButton)) {
        _drawBoxLis = !_drawBoxLis; 
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_BracketRight) && (modifiers == Qt::NoButton)) {
        if (!animationIsStarted())
            DrawOneFrameForward(); 
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::NoButton)) {
            Sphere sphere;
            std::cout << "Sphere <x, y, z, scale>: " << std::flush; 
            PRE_CIN_CLEAR; 
            std::cin >> sphere.origin.x >> sphere.origin.y >> sphere.origin.z >> sphere.scale; 
            POST_CIN_CLEAR;
            _sphereCin.push_back(sphere); 
    }
    else if ((e->key() == Qt::Key_T) && (modifiers == Qt::NoButton)) {
      if (camera()->type() == Camera::ORTHOGRAPHIC){
        camera()->setType(Camera::PERSPECTIVE);
        std::cout << "camera: perspective\n"; 
      }
      else {
        camera()->setType(Camera::ORTHOGRAPHIC);
        std::cout << "camera: orthographic\n"; 
      }
    }
    else if ((e->key() == Qt::Key_G) && (modifiers == Qt::NoButton)) {
        _drawGround = (_drawGround + 1)%3; 
        optionsChanged = true;
    }
    else if ((e->key() == Qt::Key_S) && (modifiers == Qt::ControlModifier)) 
    {
        setSnapshotFormat("JPEG");
        setSnapshotQuality(80);
        setSnapshotFileName("frames/test");
        _takeSnapshots = !_takeSnapshots;
        std::cout << "takeSnapshots: " << std::boolalpha << _takeSnapshots << std::endl;
    }
    else if ((e->key() == Qt::Key_S) && (modifiers == Qt::ShiftModifier)) {
        std::string filename; 
        std::cout << "Enter filename for saving all slices data: " << std::flush; 
        PRE_CIN_CLEAR;
        std::cin >> filename; 
        POST_CIN_CLEAR;
        std::ofstream of(filename.c_str()); 
        if (of) 
        {
            std::cout << " Writing all slice data to file: " << filename << "\n"; 
            const int N_slices = _sliceCin.size(); 
            of << std::setprecision(16); 
            of << "# <number slices> <_sliceDataPointer>\n"; 
            of << N_slices << " " << _sliceDataPointer << "\n"; 
            for (int s_idx=0; s_idx<N_slices; ++s_idx)
            {
                auto &slice = _sliceCin.at(s_idx); 
                auto &data = slice.data;
                auto &samples = slice.samples; 
                const int N_probes = data.rows(); 
                const int N_dimension = data.cols(); 
                of << "# <slice index> <number of data probes on slice> <data dimension per probe> \n"; 
                of << s_idx << " " << N_probes << " " << N_dimension << "\n"; 
                of << "# <position_x> <position_y> <position_z> <data_1> <data_2> ... <data_dim>\n"; 
                for (int p_idx=0; p_idx<N_probes; ++p_idx)
                {
                    auto &vertex = samples.at(p_idx); 
                    of << vertex.x << " " << vertex.y << " " << vertex.z << " ";
                    for (int d_idx=0; d_idx<N_dimension; ++d_idx)
                    {
                        of << data(p_idx, d_idx) << " ";
                    }
                    of << std::endl; 
                }
            }
            of.close(); 
            std::cout << " Write complete." << std::endl;
        }
    }
    else if ((e->key() == Qt::Key_F) && (modifiers == Qt::ShiftModifier)) {
        // print all velocity BC
        std::cout << "Debug: execute some debug function" << std::endl;
        auto &sceneObjects = _simWorld->GetSceneObjects(); 
        const auto &object = sceneObjects->GetPtr(0); 
        std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object)
            ->PrintAllVelocity("allVelocityFDTD.txt", 0);
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::ShiftModifier)) {
        // write to disk 
        std::string filename; 
        std::cout << "Debug: failed reflection arrow filename: " << std::flush; 
        PRE_CIN_CLEAR;
        std::cin >> filename; 
        POST_CIN_CLEAR;
        Push_Back_ReflectionArrows(filename); 
    }
    else if ((e->key() == Qt::Key_F) && (modifiers == Qt::ControlModifier)) {
        PRE_CIN_CLEAR; 
        _fixedSliceColorMapRange = !_fixedSliceColorMapRange; 
        if (_fixedSliceColorMapRange) 
        {
            REAL cMin, cMax; 
            std::cout << "Input colormap range <cmin> <cmax>: " << std::flush; 
            while (std::cin >> cMin >> cMax)
            {
                if (cMin < cMax)
                    break; 
                std::cout << "Invalid range. Input colormap range <cmin> <cmax>: " << std::flush; 
            }
            _sliceColorMapRange.x = cMin; 
            _sliceColorMapRange.y = cMax;  
        }
        ResetSliceColormap(); 
        POST_CIN_CLEAR;
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::ControlModifier)) {
        std::string filename; 
        std::cout << "To-draw arrow filename: " << std::flush; 
        PRE_CIN_CLEAR;
        std::cin >> filename; 
        POST_CIN_CLEAR;
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
    else if ((e->key() == Qt::Key_C) && (modifiers == Qt::ControlModifier)) {
        _sliceCin.clear(); 
    }
    else if ((e->key() == Qt::Key_P) && (modifiers == Qt::AltModifier)) {
        PRE_CIN_CLEAR;
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
        POST_CIN_CLEAR;
    }
    else if ((e->key() == Qt::Key_N) && (modifiers == Qt::NoButton)) {
            PRE_CIN_CLEAR;
            Vector3f x, n; 
            std::cout << "Arrow start location <x, y, z>: " << std::flush; 
            std::cin >> x.x >> x.y >> x.z; 
            std::cout << "Arrow normal <x, y, z>: " << std::flush; 
            std::cin >> n.x >> n.y >> n.z; 
            Arrow arrow; 
            arrow.start = x; 
            arrow.normal = n; 
            _arrowCin.push_back(arrow); 
            POST_CIN_CLEAR;
    }
    else if ((e->key() == Qt::Key_M) && (modifiers == Qt::NoButton)) {
        PRE_CIN_CLEAR;
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
        POST_CIN_CLEAR;
    }
    else if ((e->key() == Qt::Key_Y) && (modifiers == Qt::NoButton)) {
        int dim; 
        Vector3d origin; 
        PRE_CIN_CLEAR;
        std::cout << "Slice dim: " << std::flush; 
        std::cin >> dim;
        std::cout << "Slice origin: " << std::flush; 
        std::cin >> origin.x >> origin.y >> origin.z; 
        POST_CIN_CLEAR; 
        AddSlice(dim, origin[dim]);
    }
    else if ((e->key() == Qt::Key_Y) && (modifiers == Qt::ShiftModifier)) {
        _sliceDataPointer = (_sliceDataPointer + 1)%8; 
        optionsChanged = true;
        SetAllSliceDataReady(false); 
    }
    else if ((e->key() == Qt::Key_Y) && (modifiers == Qt::ControlModifier)) {
        _sliceDataPointer = std::max<int>(_sliceDataPointer - 1, 0); 
        optionsChanged = true;
        SetAllSliceDataReady(false); 
    }
    else if ((e->key() == Qt::Key_C) && (modifiers == Qt::NoButton)) {
            _sphereCin.clear(); 
    }
    else if ((e->key() == Qt::Key_C) && (modifiers == Qt::ShiftModifier)) {
            _arrowCin.clear(); 
    }
    else if ((e->key() == Qt::Key_R) && (modifiers == Qt::ControlModifier)) {
        PRE_CIN_CLEAR; 
        std::string inputFile, outputFile; 
        std::cout << "FBem input file: " << std::flush;
        std::cin >> inputFile; 
        std::cout << "FBem output file: " << std::flush;
        std::cin >> outputFile; 
        std::cout << std::endl;

        REAL frequency; 
        std::cout << "Frequency for this mode (Hz): " << std::flush;
        std::cin >> frequency; 
        std::cout << std::endl;
        POST_CIN_CLEAR; 

        _bemSolver->AddFBemSolution(inputFile, outputFile, 2.0*M_PI*frequency); 

        // always points to the latest mode
        _bemModePointer = _bemSolver->N_Modes()-1;  
        SetAllSliceDataReady(false); 
    }
    else if ((e->key() == Qt::Key_R) && (modifiers == Qt::NoButton)) {
        DrawOneFrameForward(); 
    }
    else if ((e->key() == Qt::Key_Right) && (modifiers == Qt::NoButton)) {
        MoveSceneCenter(0, _solverSettings->cellSize*0.9); 
    }
    else if ((e->key() == Qt::Key_Left) && (modifiers == Qt::NoButton)) {
        MoveSceneCenter(0, -_solverSettings->cellSize*0.9); 
    }
    else if ((e->key() == Qt::Key_Up) && (modifiers == Qt::NoButton)) {
        MoveSceneCenter(2, _solverSettings->cellSize*0.9); 
    }
    else if ((e->key() == Qt::Key_Down) && (modifiers == Qt::NoButton)) {
        MoveSceneCenter(2, -_solverSettings->cellSize*0.9); 
    }
    else {
        handled = false; 
    }

    // still enable the default qglviewer event handling but this function has
    // priority
    if (!handled)
        QGLViewer::keyPressEvent(e);
    if (optionsChanged) 
    {
        PrintDrawOptions();
        PrintFrameInfo(); 
    }
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
    if (selectedName() != -1) 
    {
        const int cell_idx = selectedName(); 
        bool found; 
        qglviewer::Vec selectedPoint = camera()->pointUnderPixel(point, found);
        Vector3d x(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
        _listenedUnit.reset();
        for (auto &s : _sliceCin)
        {
            if (s.intersectingUnit->GetBoundingBox().isInside(x))
            {
                _listenedUnit = s.intersectingUnit;
            }
        }
        if (_listenedUnit)
        {
            MAC_Grid::Cell cell; 
            _listenedUnit->simulator->GetSolver()->FetchCell(cell_idx, cell); 
            std::cout << cell << std::endl; 
            _listenedCell = cell; 
        }
    }
    else
    {
        _listenedCell.index = -1; // reset 
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
InitializeBEMSolver()
{
    //FIXME debug: 1) always set mesh id to 0; 2) hard-coded frequency and read
    //files
    if (!_bemSolver)
        _bemSolver = std::make_shared<KirchhoffIntegralSolver>(); 
    int bemMeshID = 0; 
    //std::cout << "Input BEM solution mesh id in the simulator: " << std::flush; 
    //std::cin >> bemMeshID; 
    std::shared_ptr<TriangleMesh<REAL> > bemMesh = _simWorld->GetSceneObjects()->GetPtr(bemMeshID)->GetMeshPtr();
    _bemSolver->SetMesh(bemMesh);
    std::cout << " Set BEM solution corresponding mesh to " << _simWorld->GetSceneObjects()->GetMeshName(bemMeshID) << std::endl;

    const REAL frequency = 1020.01;
    const std::string inputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/input-0_0.txt"); 
    const std::string outputFile("/home/jui-hsien/code/acoustics/work/plate_drop_long/fastbem/ret-0_0.txt"); 
    _bemSolver->AddFBemSolution(inputFile, outputFile, 2.0*M_PI*frequency); 
    _bemModePointer = _bemSolver->N_Modes()-1;  
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
AddSlice(const int &dim, const REAL &offset)
{
    const auto &simUnits = _simWorld->GetActiveSimUnits(); 
    for (const auto &unit : simUnits)
    {
        const BoundingBox &bbox = unit->GetBoundingBox(); 
        if (offset>bbox.axismax(dim) && offset<bbox.axismin(dim))
            continue; // this plane not intersecting the box
        Slice slice; 
        slice.dim = dim; 
        slice.origin.set(0, 0, 0);
        slice.origin[dim] += offset; 
        slice.dataReady = false;
        slice.intersectingUnit = unit; 
        ConstructSliceSamples(slice);
        _sliceCin.push_back(slice); 
    }
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
    samples.clear(); 
    gridLines.clear(); 

    const auto &settings = _solverSettings; 
    const auto &simulator = slice.intersectingUnit->simulator; 
    const Tuple3i &divisions = simulator->GetGrid().pressureFieldDivisions(); 
    const REAL cellSize = settings->cellSize; 

    //const REAL halfLength = (REAL)division*cellSize / 2.0; 
    //const REAL minBound = -halfLength + cellSize/2.0; 
    const BoundingBox pBBox = simulator->GetGrid().PressureBoundingBox(); 
    slice.minBound = pBBox.minBound() + 0.5*cellSize;
    slice.maxBound = pBBox.maxBound() - 0.5*cellSize; 
    slice.N_sample_per_dim = divisions; 

    const int dim_0 = (dim + 1) % 3; 
    const int dim_1 = (dim + 2) % 3; 
    for (int dim_0_idx=0; dim_0_idx<divisions[dim_0]; ++dim_0_idx)
        for (int dim_1_idx=0; dim_1_idx<divisions[dim_1]; ++dim_1_idx)
        {
            Vector3d sample; 
            sample(dim) = origin(dim); 
            sample(dim_0) = slice.minBound[dim_0] + cellSize*(REAL)dim_0_idx; 
            sample(dim_1) = slice.minBound[dim_1] + cellSize*(REAL)dim_1_idx; 
            samples.push_back(sample); 
        }

    // horizontal grid lines
    const REAL xStart = slice.minBound[dim_0] - 0.5*cellSize; 
    const REAL xStop =  slice.maxBound[dim_0] + 0.5*cellSize; 
    const REAL yStart = slice.minBound[dim_1] - 0.5*cellSize; 
    const REAL yStop =  slice.maxBound[dim_1] + 0.5*cellSize; 
    for (int dim_0_idx=0; dim_0_idx<divisions[dim_0]+1; ++dim_0_idx)
    {
        Vector3d start; 
        Vector3d stop; 
        start(dim) = origin(dim); 
        stop(dim) = origin(dim); 
        start(dim_0) = xStart + cellSize*(REAL)dim_0_idx;
        stop(dim_0) = xStart + cellSize*(REAL)dim_0_idx;
        start(dim_1) = yStart; 
        stop(dim_1) = yStop; 
        gridLines.push_back(start); 
        gridLines.push_back(stop); 
    }

    // vertical grid lines
    for (int dim_1_idx=0; dim_1_idx<divisions[dim_1]+1; ++dim_1_idx)
    {
        Vector3d start; 
        Vector3d stop; 
        start(dim) = origin(dim); 
        stop(dim) = origin(dim); 
        start(dim_0) = xStart; 
        stop(dim_0) = xStop; 
        start(dim_1) = yStart + cellSize*(REAL)dim_1_idx;
        stop(dim_1) = yStart + cellSize*(REAL)dim_1_idx;
        gridLines.push_back(start); 
        gridLines.push_back(stop); 
    }

    // make colormap for this slice.
    if (!_sliceColorMap)
        //_sliceColorMap = std::make_shared<JetColorMap>(); 
        _sliceColorMap = std::make_shared<DipoleColorMap>(); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
ComputeAndCacheSliceData(const int &dataPointer)
{
    bool updateColormap = false; 
    const int N_slices = _sliceCin.size();
    REAL maxCoeffSlices = std::numeric_limits<REAL>::min(); 
    REAL minCoeffSlices = std::numeric_limits<REAL>::max(); 
    // first fetch all data and compute min, max
    for (int s_idx=0; s_idx<N_slices; ++s_idx) 
    {
        auto &slice = _sliceCin.at(s_idx); 
        auto &unit = slice.intersectingUnit; 
        if (unit->boxCenterChanged)
        {
            ConstructSliceSamples(slice); 
            slice.dataReady = false; 
            unit->boxCenterChanged = false; // reset
        }
        if (slice.dataReady)
            continue; 
        Eigen::MatrixXd &data = slice.data; 
        data.setZero();
        auto simulator = slice.intersectingUnit->simulator; 
        if (dataPointer == 0)
        {
            simulator->GetSolver()->FetchPressureData(slice.samples, data);
        } 
        else if (dataPointer == 1)
        {
            // 0.0: bulk cell; 1: solid cell; -1: ghost cell
            simulator->GetSolver()->FetchPressureCellType(slice.samples, data, _sceneBox);
        } 
        else if (dataPointer == 2)
        {
            simulator->GetSolver()->FetchVelocityData(slice.samples, 0, data);
        }
        else if (dataPointer == 3)
        {
            simulator->GetSolver()->FetchVelocityData(slice.samples, 1, data);
        }
        else if (dataPointer == 4)
        {
            simulator->GetSolver()->FetchVelocityData(slice.samples, 2, data);
        }
        else if (dataPointer == 5)
        {
            simulator->GetSolver()->FetchPressureData(slice.samples, data, 0);
        } 
        else if (dataPointer == 6)
        {
            simulator->GetSolver()->FetchPressureData(slice.samples, data, 1);
        } 
        else if (dataPointer == 7)
        {
            simulator->GetSolver()->FetchPressureData(slice.samples, data, 2);
        } 
        else if (dataPointer == 8 || dataPointer == 9) 
        {
            const int N_samples = slice.samples.size(); 
            data.resize(N_samples, 1);
            int count = 0;
            for (int d_idx=0; d_idx<N_samples; ++d_idx)
            {
                if (dataPointer == 8) 
                {
                    const std::complex<REAL> transferValue = _bemSolver->Solve(_bemModePointer, slice.samples.at(d_idx));
                    data(d_idx, 0) = std::abs(transferValue);
                }
                else if (dataPointer == 9) 
                {
                    // if distance > threashold, computes transfer residual
                    REAL transferResidual; 
                    const REAL distance = simulator->GetSceneObjects()->LowestObjectDistance(slice.samples.at(d_idx)); 
                    if (distance > 0.005)
                        _bemSolver->TestSolver(_bemModePointer, _bemSolver->GetMode_k(_bemModePointer), slice.samples.at(d_idx), transferResidual);
                    else
                        transferResidual = 0.0; 
                    data(d_idx, 0) = transferResidual;
                }
                count ++; 
                std::cout << "\r" << (REAL)count / (REAL)N_samples * 100.0 << "\% completed" << std::flush;
            }
            std::cout << std::endl;
        }
        else 
        {
            throw std::runtime_error("**ERROR** dataPointer out of range");
        }

        maxCoeffSlices = max(maxCoeffSlices, data.maxCoeff()); 
        minCoeffSlices = min(minCoeffSlices, data.minCoeff()); 
        updateColormap = true; 
        slice.dataReady = true;
    }

    if (updateColormap && !_fixedSliceColorMapRange)
    {
        minCoeffSlices = (minCoeffSlices==maxCoeffSlices ?  maxCoeffSlices-EPS : minCoeffSlices);
        if (dataPointer == 1) 
        {
            _sliceColorMap->set_interpolation_range(-1.0, 1.0); 
            _sliceColorMapRange.x = -1; 
            _sliceColorMapRange.y =  1; 
        }
        else
        {
            _sliceColorMap->set_interpolation_range(minCoeffSlices, maxCoeffSlices); 
            _sliceColorMapRange.x = minCoeffSlices; 
            _sliceColorMapRange.y = maxCoeffSlices; 
        }
    }
    _messageColormap = "Colormap range = [" + QString::number(_sliceColorMapRange.x) + ", " + QString::number(_sliceColorMapRange.y) + "]"; 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
DrawOneFrameForward()
{
    if (_previewSpeed == 0)
    {
        _currentFrame++;
        _simWorld->StepWorld(); 
        if (_listenedCell.index >= 0 && _listenedUnit)
        {
            _listenedUnit->simulator->GetSolver()->FetchCell(_listenedCell.index, _listenedCell); 
            std::cout << _listenedCell << std::endl;
        }
#if DEBUG_WRITE_REFLECTION_ARROWS_INTERVAL > 0
        if (_currentFrame % DEBUG_WRITE_REFLECTION_ARROWS_INTERVAL == 0)
            Push_Back_ReflectionArrows("a"); 
#endif
        SetAllSliceDataReady(false); 
    }
    else
    {
        _currentFrame += _previewSpeed;
        _simWorld->PreviewStepping(_previewSpeed);
    }
    if (_takeSnapshots)
        saveSnapshot(true, true);
    PrintFrameInfo();
    updateGL(); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
RestoreDefaultDrawOptions()
{
    _wireframe = (_remoteConnection ? 4 : 2); 
    _sliceWireframe.reset(); 
    _sliceWireframe.set(1); // draw face only
    _drawBoxLis = true; 
    _drawGround = (_remoteConnection ? 0 : 2); 
    _drawHashedCells = false;
    _sliceDataPointer = 0; 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
PrintDrawOptions()
{
    std::cout << "\n" << std::boolalpha
              << "Draw Options \n"
              << "------------\n"
              << " Draw simulation box      : " << _drawBoxLis << "\n"
              << " Draw ground              : " << _drawGround << "\n"
              << " Draw wireframe only      : " << _wireframe << "\n"
              << " Draw hashed cells        : " << _drawHashedCells << "\n"
              << " Draw slice wireframe only: " << _sliceWireframe.to_ulong() << "\n"
              << " Draw slice data pointer  : " << _sliceDataPointer << "\n"
              << "\n"; 
    std::cout << std::flush;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
PrintFrameInfo()
{
    const int N_dataType = 10; 
    static const QString dataString[] = {"p_full", "cell_id", "v_x", "v_y", "v_z", "p_x", "p_y", "p_z", "freq_transfer", "freq_transfer_residual"}; 
    int p=0;
    for (; p<N_dataType; ++p)
        if (p == _sliceDataPointer)
            break; 
    //const std::string frameInfo("Current Frame: " + std::to_string(_currentFrame)); 
    _message = QString("");
    _message += "Current Frame: " + QString::number(_currentFrame) + "; "; 
    _message += "Current Time: " + QString::number(_simWorld->GetWorldTime()) + "; "; 
    _message += "Current Data: " + dataString[p] + "; ";
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
Push_Back_ReflectionArrows(const std::string &filename)
{
    auto &sceneObjects = _simWorld->GetSceneObjects(); 
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

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Viewer::
MoveSceneCenter(const int &dim, const double &displacement) 
{
    // FIXME debug
    ////if (!_sceneBox) 
    ////{
    ////    const BoundingBox bbox = _simulator->GetGrid().PressureBoundingBox();
    ////    _sceneBox = new BoundingBox(bbox.minBound(), bbox.maxBound()); 
    ////}
    ////(_sceneBox->minBound())[dim] += displacement; 
    ////(_sceneBox->maxBound())[dim] += displacement; 
    ////std::cout << _sceneBox->center() << std::endl; 
    ////_simulator->GetGrid().UpdatePML(*_sceneBox); 
    //Vector3d move; 
    //move[dim] = displacement;
    //const bool changed = _simulator->MoveSimBox(move); 
    //if (changed)
    //{
    //    // recompute slices
    //    for (auto &slice : _sliceCin)
    //    {
    //        slice.samples.clear(); 
    //        slice.gridLines.clear(); 
    //        slice.cells.clear(); 
    //        ConstructSliceSamples(slice); 
    //        slice.dataReady = false; 
    //    }
    //    updateGL();
    //}
}
