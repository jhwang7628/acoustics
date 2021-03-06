//////////////////////////////////////////////////////////////////////
// WaveViewer.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef WAVE_VIEWER_H
#define WAVE_VIEWER_H

#include <QtGui>

#if QT_VERSION >= 0x050000
    #include <QtWidgets>
#endif

#include <QGLViewer/qglviewer.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <geometry/SimplePointSet.hpp>

#include <utils/Evaluator.h>

#include <wavesolver/WaveSolver.h>

#include <TYPES.h>
#include <string>

#include <sndgen/player/SineWavePlayer.h>
#include <multipole/MultipolePlayer.h>

typedef struct ExtraGLData
{
    std::string test; 
    SimplePointSet<REAL> * pointset1;
    SimplePointSet<REAL> * pointset2;
} ExtraGLData; 

//////////////////////////////////////////////////////////////////////
// WaveViewer class
//
// OpenGL viewer for visualizing wave phenomena on a finite
// difference grid
//////////////////////////////////////////////////////////////////////
class WaveViewer : public QGLViewer {
    private:
        // Convenience definition
        typedef TriangleMesh<REAL>  TriMesh;
        typedef SimplePointSet<REAL> * PointSet; 

    public:
        WaveViewer( Solver &solver );

        WaveViewer( Solver &solver, ExtraGLData * extraGLData ); 

        // Destructor
        virtual ~WaveViewer();

        inline const Tuple3i    &fieldDivisions()
        {
            return _solver.fieldDivisions();
        }

        void                     setAccelerationFunction(
                BoundaryEvaluator *acceleration )
        {
            _acceleration = acceleration;
        }

        void setDrawRange( const Tuple3i &minBound, const Tuple3i &maxBound );

        void setDrawColourRange( REAL colourRange );

    protected:
        // Draws the desired range of the current finite difference grid
        virtual void draw();
        virtual void init();

        virtual void animate();

        // Define shortcut for stepping system
        virtual void keyPressEvent( QKeyEvent *e );

    private:
        // Draws walls of the finite difference grid along one axis
        void drawGridWalls( GridDirection direction );

        void drawPlane( int planeDirection, int uDirection, int vDirection,
                        int planeIndex, const Vector3d &normal );

        void drawMesh();
        void drawExtraGLData(); 
        void drawReceivers();

#if 0
        void drawGhostCells();
        void drawInterfacialCells();
#endif

        // Computes a red/green (negative/positive) colour for the
        // given pressure value
        Vector3d computeVertexColour( const Tuple3i &index );

        // Get the current pressure associated with the given vertex
        void vertexPressure( const Tuple3i &index, VECTOR &pressure );
        virtual void postSelection(const QPoint& point);
    private:
        Solver                  &_solver;
        ExtraGLData             *_extraGLData; 

        BoundaryEvaluator       *_acceleration;

        Tuple3i                  _minDrawBound;
        Tuple3i                  _maxDrawBound;

        REAL                     _drawColourMax;

        bool                     _wireframe;

        bool                     _drawMesh;
        bool                     _drawReceivers;

        bool                     _drawGhostCells;
        bool                     _drawInterfacialCells;

        int                      _drawExtraGLData; 

        int                      _drawField;
        VECTOR                   _drawPressure;

        qglviewer::Vec orig, dir, selectedPoint;
        REAL _amplitude;
        REAL _frequency;
        REAL _phase;
        REAL _estAmplitude;
        REAL _estPhase;


};

//////////////////////////////////////////////////////////////////////
// WaveWindow class
//
// Main application window storing a wave viewer and some controls
//////////////////////////////////////////////////////////////////////
class WaveWindow : public QObject {
    Q_OBJECT

    public:
        WaveWindow( Solver &solver );
        WaveWindow( Solver &solver , ExtraGLData* extraGLData); 

        // Destructor
        virtual ~WaveWindow();

        QWidget *createWindow();

        void                     setAccelerationFunction(
                BoundaryEvaluator *acceleration )
        {
            _viewer.setAccelerationFunction( acceleration );
        }

    protected:

        public slots:
            void rangeValueChanged( int x );

        void drawRangeValueChanged( double x );

    private:
        //////////////////////////////////////////////////////////////////////
        // UI widgits
        //////////////////////////////////////////////////////////////////////

        // Main layout
        QHBoxLayout              _layout;

        // Layout for the control pane
        QVBoxLayout              _controlLayout;

        // Individual control layouts
        QHBoxLayout              _xControlLayout;
        QHBoxLayout              _yControlLayout;
        QHBoxLayout              _zControlLayout;

        // Spin boxes for controlling drawing region, and labels
        QSpinBox                 _xMinBox;
        QSpinBox                 _xMaxBox;
        QSpinBox                 _yMinBox;
        QSpinBox                 _yMaxBox;
        QSpinBox                 _zMinBox;
        QSpinBox                 _zMaxBox;

        QLabel                   _xLabel;
        QLabel                   _yLabel;
        QLabel                   _zLabel;

        // For controlling wave drawing
        QHBoxLayout              _drawingLayout;
        QDoubleSpinBox           _absMaxBox;
        QLabel                   _absMaxLabel;

        // OpenGL drawing region
        WaveViewer               _viewer;

};

#endif
