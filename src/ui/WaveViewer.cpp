//////////////////////////////////////////////////////////////////////
// WaveViewer.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "WaveViewer.h"

#include <geometry/TriangleMesh.hpp>

#include <utils/IO.h>

#include <GL/glut.h>

//////////////////////////////////////////////////////////////////////
// WaveViewer implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
WaveViewer::WaveViewer( Solver &solver )
    : _solver( solver ),
      _acceleration( NULL ),
      _minDrawBound( 0, 0, 0 ),
      _maxDrawBound( solver.fieldDivisions()[ 0 ],
                     solver.fieldDivisions()[ 1 ],
                     solver.fieldDivisions()[ 2 ] ),
      _drawColourMax( 1.0 ),
      _wireframe( false ),
      _drawMesh( false ),
      _drawReceivers( false ),
      _drawGhostCells( false ),
      _drawInterfacialCells( false ),
      _drawField( 0 ),
      _drawPressure( _solver.N() )
{
    for ( int obj_idx = 0; obj_idx < _solver.meshes().size(); obj_idx++ )
    {
        printf( "(WaveViewer) Solver mesh %d has %d vertices and %d triangles\n",
                obj_idx, (int)_solver.meshes()[ obj_idx ]->vertices().size(),
                (int)_solver.meshes()[ obj_idx ]->triangles().size() );
    }

    Vector3d sceneCenter = solver.sceneCenter();
    qglviewer::Vec c( sceneCenter[ 0 ], sceneCenter[ 1 ], sceneCenter[ 2 ] );
    setSceneCenter( c );
    setSceneRadius( solver.fieldDiameter() );

    glEnable( GL_COLOR_MATERIAL );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
WaveViewer::~WaveViewer()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::setDrawRange( const Tuple3i &minBound,
                               const Tuple3i &maxBound )
{
    const Tuple3i             &divs = _solver.fieldDivisions();

    _minDrawBound = minBound;
    _maxDrawBound = maxBound;

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        _minDrawBound[ dimension ] = max( _minDrawBound[ dimension ], 0 );
        _maxDrawBound[ dimension ] = min( _maxDrawBound[ dimension ],
                divs[ dimension ] );
    }

    updateGL();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::setDrawColourRange( REAL colourRange )
{
    _drawColourMax = colourRange;

    updateGL();
}

//////////////////////////////////////////////////////////////////////
// Draws the desired range of the current finite difference grid
//////////////////////////////////////////////////////////////////////
void WaveViewer::draw()
{
    Tuple3i                    drawSize = _maxDrawBound - _minDrawBound;

    if ( ( drawSize[ 0 ] < 2 && drawSize[ 1 ] < 2 )
            || ( drawSize[ 0 ] < 2 && drawSize[ 2 ] < 2 )
            || ( drawSize[ 1 ] < 2 && drawSize[ 2 ] < 2 ) )
    {
        return;
    }

    drawGridWalls( GRID_X );
    drawGridWalls( GRID_Y );
    drawGridWalls( GRID_Z );

    if ( _drawMesh )
    {
        drawMesh();
    }

    if ( _drawReceivers )
    {
        drawReceivers();
    }

#if 0
    if ( _drawGhostCells )
    {
        drawGhostCells();
    }
#endif

#if 0
    if ( _drawInterfacialCells )
    {
        drawInterfacialCells();
    }
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::init()
{
    restoreStateFromFile();

    QColor bgColor( 0.05, 0.05, 0.05 );
    setBackgroundColor( bgColor );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::animate()
{
    if ( !_solver.stepSystem( *_acceleration ) )
    {
        // Disable animation if told to by the solver
        toggleAnimation();
    }

    QGLViewer::animate();
}

//////////////////////////////////////////////////////////////////////
// Define shortcut for stepping system
//////////////////////////////////////////////////////////////////////
void WaveViewer::keyPressEvent( QKeyEvent *e )
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    bool                        handled = false;

    if ( e->key() == Qt::Key_S )
    {
        _solver.stepSystem( *_acceleration );

        handled = true;

        updateGL();
    }
    else if ( e->key() == Qt::Key_W )
    {
        _wireframe = !_wireframe;

        if ( _wireframe )
        {
            glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
        }
        else
        {
            glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        }

        handled = true;

        updateGL();
    }
    else if ( e->key() == Qt::Key_N )
    {
        _solver.writeWaveOutput();
    }
    else if ( e->key() == Qt::Key_D )
    {
        _drawMesh = !_drawMesh;

        handled = true;

        updateGL();
    }
    else if ( e->key() == Qt::Key_L )
    {
        _drawReceivers = !_drawReceivers;

        handled = true;

        updateGL();
    }
    else if ( e->key() == Qt::Key_A )
    {
        toggleAnimation();

        handled = true;
    }
    else if ( e->key() == Qt::Key_G )
    {
        // Draw ghost cells
        _drawGhostCells = !_drawGhostCells;

        handled = true;
    }
    else if ( e->key() == Qt::Key_I )
    {
        // Draw interfacial cells
        _drawInterfacialCells = !_drawInterfacialCells;

        handled = true;
    }
    else if ( e->key() >= Qt::Key_1 && e->key() <= Qt::Key_6 )
    {
        _drawField = e->key() - Qt::Key_1;

        if ( _drawField >= _solver.N() )
        {
            _drawField = _solver.N() - 1;
        }

        handled = true;

        updateGL();
    }

    if ( !handled )
    {
        QGLViewer::keyPressEvent( e );
    }
}

//////////////////////////////////////////////////////////////////////
// Draws walls of the finite difference grid along one axis
//////////////////////////////////////////////////////////////////////
void WaveViewer::drawGridWalls( GridDirection direction )
{
    Tuple3i                    drawSize = _maxDrawBound - _minDrawBound;

    int                        plane1 = _minDrawBound[ direction ];
    int                        plane2 = _maxDrawBound[ direction ] - 1;

    int                        uDirection;
    int                        vDirection;

    Tuple3i                    vertexIndex;
    Vector3d                   vertexPosition;
    Vector3d                   normal;

    if ( drawSize[ ( direction + 1 ) % NUM_DIRECTIONS ] < 2
            || drawSize[ ( direction + 2 ) % NUM_DIRECTIONS ] < 2 )
    {
        return;
    }

    uDirection = ( direction + 1 ) % NUM_DIRECTIONS;
    vDirection = ( direction + 2 ) % NUM_DIRECTIONS;

    vertexIndex[ direction ] = plane1;

    normal[ direction ] = -1.0;

    drawPlane( direction, uDirection, vDirection, plane1, normal );

    if ( plane1 >= plane2 )
        return;

    vertexIndex[ direction ] = plane2;
    normal[ direction ] = 1.0;

    drawPlane( direction, uDirection, vDirection, plane2, normal );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::drawPlane( int planeDirection, int uDirection, int vDirection,
                            int planeIndex, const Vector3d &normal )
{
    Tuple3i                    vertexIndex;
    Vector3d                   vertexPosition;
    Vector3d                   colour;

    float                      planeColour[ 4 ] = { 0.2, 0.2, 0.2,  1.0 };
    float                      specColour[ 4 ] = { 0.2, 0.2, 0.2, 1.0 };

    vertexIndex[ planeDirection ] = planeIndex;

    glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, planeColour );
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specColour );
    glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 2.0 );
    for ( int i = _minDrawBound[ uDirection ];
            i < _maxDrawBound[ uDirection ] - 1; i++ )
    {
        glBegin( GL_QUAD_STRIP );
        for ( int j = _minDrawBound[ vDirection ];
                j < _maxDrawBound[ vDirection ]; j++ )
        {
            vertexIndex[ uDirection ] = i;
            vertexIndex[ vDirection ] = j;

            colour = computeVertexColour( vertexIndex );

            //glColor3f( 0.8, 0.0, 0.3 );
            glColor3f( colour[ 0 ], colour[ 1 ], colour[ 2 ] );
            glNormal3f( normal[ 0 ], normal[ 1 ], normal[ 2 ] );

            vertexPosition = _solver.fieldPosition( vertexIndex );
            glVertex3f( vertexPosition[ 0 ],
                    vertexPosition[ 1 ],
                    vertexPosition[ 2 ] );

            vertexIndex[ uDirection ] += 1;

            colour = computeVertexColour( vertexIndex );

            glColor3f( colour[ 0 ], colour[ 1 ], colour[ 2 ] );

            vertexPosition = _solver.fieldPosition( vertexIndex );
            glVertex3f( vertexPosition[ 0 ],
                    vertexPosition[ 1 ],
                    vertexPosition[ 2 ] );
        }
        glEnd();
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::drawMesh()
{
#if 0
    const TriMesh             &mesh = _solver.mesh();
#endif
    const vector<const TriMesh *>   &meshes = _solver.meshes();

    float                            meshColour[4] = { 0.5, 0.0, 0.5, 1.0 };
    float                            specColour[ 4 ] = { 0.1, 0.1, 0.1, 1.0 };

    for ( int mesh_idx = 0; mesh_idx < meshes.size(); mesh_idx++ )
    {
        const TriMesh               &mesh = *( meshes[ mesh_idx ] );
        const vector<Point3d>       &vertices = mesh.vertices();
        const vector<Tuple3ui>      &triangles = mesh.triangles();
        const vector<Vector3d>      &normals = mesh.normals();

        glBegin( GL_TRIANGLES );
        glColor3fv( meshColour );
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, meshColour );
        glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specColour );
        glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 10.0 );
        for ( int i = 0; i < triangles.size(); i++ )
        {
#if 0
            //glColor3f( 0.5, 0.0, 0.5 );
            const Vector3d          &n = triangles[ i ]->getNormal();

            for ( int j = 0; j < 3; j++ )
            {
                const Vector3d        &v = triangles[ i ]->getX( j );

                glNormal3f( n[ 0 ], n[ 1 ], n[ 2 ] );
                glVertex3f( v[ 0 ], v[ 1 ], v[ 2 ] );
            }
#endif
            for ( int j = 0; j < 3; j++) {
                const Point3d       &v = vertices[ triangles[ i ][ j ] ];
                const Vector3d      &n = normals[ triangles[ i ][ j ] ];

                glNormal3f( n[ 0 ], n[ 1 ], n[ 2 ] );
                glVertex3f( v[ 0 ], v[ 1 ], v[ 2 ] );
            }
        }
        glEnd();
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveViewer::drawReceivers()
{
    const Vector3Array        *receivers = _solver.listeningPositions();

    REAL                       radius = _solver.fieldDiameter() * 0.01;

    if ( !receivers )
    {
        return;
    }

    float                      sphereColour[4] = { 0.0, 0.75, 0.75, 1.0 };
    float                      specColour[ 4 ] = { 0.1, 0.1, 0.1, 1.0 };

    for ( int listen_idx = 0; listen_idx < receivers->size(); listen_idx++ )
    {
        const Vector3d          &pos = receivers->at( listen_idx );

        //glColor3f( 0.0, 0.75, 0.75 );
        glPushMatrix();
        glTranslated( pos[ 0 ], pos[ 1 ], pos[ 2 ] );
        glColor3fv( sphereColour );
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, sphereColour );
        glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specColour );
        glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 10.0 );
        glutSolidSphere( radius, 10, 20 );
        glPopMatrix();
    }
}

//////////////////////////////////////////////////////////////////////
// Computes a red/green (negative/positive) colour for the
// given pressure value
//////////////////////////////////////////////////////////////////////
Vector3d WaveViewer::computeVertexColour( const Tuple3i &index )
{
    REAL                       pressure;

    Vector3d                   colour( 0.2, 0.2, 0.2 );

    REAL                       absDiff;

    vertexPressure( index, _drawPressure );

    pressure = _drawPressure( _drawField );

    absDiff = abs( pressure );
    absDiff = min( absDiff, _drawColourMax );
    absDiff /= _drawColourMax;
    absDiff *= 0.8;

    if ( pressure < 0 )
    {
        colour[ 0 ] += absDiff;
    }
    else
    {
        colour[ 1 ] += absDiff;
    }

    return colour;
}

//////////////////////////////////////////////////////////////////////
// Get the current pressure associated with the given vertex 
//////////////////////////////////////////////////////////////////////
void WaveViewer::vertexPressure( const Tuple3i &index, VECTOR &pressure )
{
    _solver.vertexPressure( index, pressure );
}

//////////////////////////////////////////////////////////////////////
// WaveWindow implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
WaveWindow::WaveWindow( Solver &solver )
    : _xLabel( QApplication::translate( "xlabel", "Grid range (x):" ) ),
      _yLabel( QApplication::translate( "xlabel", "Grid range (y):" ) ),
      _zLabel( QApplication::translate( "xlabel", "Grid range (z):" ) ),
      _absMaxLabel( QApplication::translate( "maxlabel", "Drawing maximum" ) ),
      _viewer( solver )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
WaveWindow::~WaveWindow()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
QWidget *WaveWindow::createWindow()
{
    QWidget                   *window = new QWidget();

    const Tuple3i             &divs = _viewer.fieldDivisions();

    _xMinBox.setRange( 0, divs[ 0 ] );
    _yMinBox.setRange( 0, divs[ 1 ] );
    _zMinBox.setRange( 0, divs[ 2 ] );

    _xMaxBox.setRange( 0, divs[ 0 ] );
    _yMaxBox.setRange( 0, divs[ 1 ] );
    _zMaxBox.setRange( 0, divs[ 2 ] );

    _xMaxBox.setValue( divs[ 0 ] );
    _yMaxBox.setValue( divs[ 1 ] );
    _zMaxBox.setValue( divs[ 2 ] );

    // Connect signals
    connect( &_xMinBox, SIGNAL( valueChanged( int ) ), this,
            SLOT( rangeValueChanged( int ) ) );
    connect( &_xMaxBox, SIGNAL( valueChanged( int ) ), this,
            SLOT( rangeValueChanged( int ) ) );
    connect( &_yMinBox, SIGNAL( valueChanged( int ) ), this,
            SLOT( rangeValueChanged( int ) ) );
    connect( &_yMaxBox, SIGNAL( valueChanged( int ) ), this,
            SLOT( rangeValueChanged( int ) ) );
    connect( &_zMinBox, SIGNAL( valueChanged( int ) ), this,
            SLOT( rangeValueChanged( int ) ) );
    connect( &_zMaxBox, SIGNAL( valueChanged( int ) ), this,
            SLOT( rangeValueChanged( int ) ) );

    _xControlLayout.addWidget( &_xLabel );
    _xControlLayout.addWidget( &_xMinBox );
    _xControlLayout.addWidget( &_xMaxBox );

    _yControlLayout.addWidget( &_yLabel );
    _yControlLayout.addWidget( &_yMinBox );
    _yControlLayout.addWidget( &_yMaxBox );

    _zControlLayout.addWidget( &_zLabel );
    _zControlLayout.addWidget( &_zMinBox );
    _zControlLayout.addWidget( &_zMaxBox );

    _absMaxBox.setRange( 1e-10, 1e6 );
    _absMaxBox.setValue( 1.0 );
    _absMaxBox.setDecimals( 10 );

    connect( &_absMaxBox, SIGNAL( valueChanged( double ) ), this,
            SLOT( drawRangeValueChanged( double ) ) );

    _drawingLayout.addWidget( &_absMaxLabel );
    _drawingLayout.addWidget( &_absMaxBox );

    _controlLayout.addLayout( &_xControlLayout );
    _controlLayout.addLayout( &_yControlLayout );
    _controlLayout.addLayout( &_zControlLayout );
    _controlLayout.addLayout( &_drawingLayout );

    _layout.addLayout( &_controlLayout );
    _layout.addWidget( &_viewer );

    window->setLayout( &_layout );
    window->resize( 1024, 768 );
    window->setWindowTitle( QApplication::translate( "waveviewer",
                "Wave Viewer" ) );

    return window;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveWindow::rangeValueChanged( int x )
{
    Tuple3i                    minBound( _xMinBox.value(),
            _yMinBox.value(),
            _zMinBox.value() );

    Tuple3i                    maxBound( _xMaxBox.value(),
            _yMaxBox.value(),
            _zMaxBox.value() );

    _viewer.setDrawRange( minBound, maxBound );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveWindow::drawRangeValueChanged( double x )
{
    _viewer.setDrawColourRange( _absMaxBox.value() );
}