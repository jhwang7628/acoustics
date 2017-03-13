//////////////////////////////////////////////////////////////////////
// RigidMesh.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "RigidMesh.h"

#include <utils/IO.h>
#include <utils/STLUtil.h>
#include <utils/trace.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
RigidMesh::RigidMesh( const TriangleMesh<REAL> &mesh, const string &rigidFilePrefix,
                      REAL density )
    : _mesh( mesh ),
      _density( density )
{
    char                       buf[ 1024 ]; 
    FloatArray                 massData;

    sprintf( buf, "%s_inertia.3matrix", rigidFilePrefix.c_str() );
    _inertia.read( buf, true /* File is in row-major format */ );
    _inertia *= density;
    _inertiaInverse = _inertia.inverse();

    sprintf( buf, "%s_centerOfMass.3vector", rigidFilePrefix.c_str() );
    _centerOfMass.read( buf );

    sprintf( buf, "%s_mass.vector", rigidFilePrefix.c_str() );
    readVector( buf, massData );
    TRACE_ASSERT( massData.size() == 1 );
    _mass = massData[ 0 ];
    _mass *= density;

    printf( "Rigid mesh has mass %e\n", _mass );
    printf( "Rigid mesh has volume %e\n", _mass / density );

    massData.clear();
    sprintf( buf, "%s_materialProperties.vector", rigidFilePrefix.c_str() );
    readVector( buf, massData );
    TRACE_ASSERT( massData.size() == 2 );
    _youngsModulus = massData[ 0 ];
    _poissonRatio = massData[ 1 ];
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
RigidMesh::~RigidMesh()
{
}

//////////////////////////////////////////////////////////////////////
// Figures out an inverse effective mass for the given
// position on the surface of the mesh, and a force applied in the
// given direction
//
//      Let x0 be the center of mass and I_x be the cross product
//      matrix for the vector (x - x0).
//      Then the inverse effective mass is:
//        (1/m) + f^T * I_x^{T} * M^{-1} * I_x * f
//      where m is the object mass and M is the inertia matrix
//////////////////////////////////////////////////////////////////////
REAL RigidMesh::inverseEffectiveMass( const Point3d &x,
                                      const Vector3d &forceDirection ) const
{
    Matrix3d                   I_x;
    Vector3d                   forceDirectionNormalized;
    Vector3d                   tmp;
    REAL                       invEffectiveMass = 0.0;

    I_x = Matrix3d::crossProductMatrix( x - _centerOfMass );

    forceDirectionNormalized = forceDirection;
    forceDirectionNormalized.normalize();

    tmp = forceDirectionNormalized;
    tmp = I_x * tmp;
    tmp = _inertiaInverse * tmp;
    tmp = I_x.transpose() * tmp;

    invEffectiveMass = forceDirectionNormalized.dotProduct( tmp );
    invEffectiveMass += 1.0 / _mass;

    return invEffectiveMass;
}

//////////////////////////////////////////////////////////////////////
// Returns acceleration directional magnitudes, given a force
// applied in the given direction to the given point
//////////////////////////////////////////////////////////////////////
void RigidMesh::accelerationCoefficients(
        const Point3d &x, const Vector3d &force,
        Vector3d &translationalAcceleration,
        Vector3d &rotationalAcceleration ) const
{
    Vector3d                   r = x - _centerOfMass;
    Vector3d                   torque = r.crossProduct( force );

    translationalAcceleration = force / _mass;

    rotationalAcceleration = _inertiaInverse * torque;
}

//////////////////////////////////////////////////////////////////////
// Fits this mesh to an ellipsoid in which the z-axis has unit length
//////////////////////////////////////////////////////////////////////
void RigidMesh::fitEllipsoid( ObjectMeasure measure )
{
    REAL                       a, b, c;
    REAL                       V;
    REAL                       A;
    REAL                       ellipseArea;
    REAL                       volumeScaling = 0.;

    const REAL                 p = 1.6075;

    // To start, we need the moment of inertia for this object
#if 0
    _inertia.symmetricEigensystem( _inertiaMoments, _inertiaAxes );
#endif
    MATRIX::symmetricEigensystem( _inertia, _inertiaMoments, _inertiaAxes );

    if ( measure == MATCH_RIGID ) {
        // Simultaneously match the rigid mass and moment inertia of
        // the proxy ellipsoid, assuming that the ellipsoid has the same
        // density as the input mesh
        REAL                     C1, C2, C3;
        Vector3d                 d = _inertiaMoments;

        d *= 15.0 / ( 4.0 * _density * M_PI );

        C1 = d[ 1 ] - d[ 0 ] + d[ 2 ];
        C1 = C1 / 2.0;

        C2 = d[ 1 ] / C1 - 1.0;
        C3 = d[ 2 ] / C1 - 1.0;

        a = d[ 0 ];
        a /= sqrt( C2 );
        a /= sqrt( C3 );
        a /= C2 + C3;
        a = pow( a, 0.2 );

        b = a * sqrt( C3 );
        c = a * sqrt( C2 );

#if 0
        REAL ellipseMass = _density * 4.0 * M_PI * a * b * c / 3.0;
        Vector3d ellipseMoments( b * b + c * c, a * a + c * c, a * a + b * b );
        ellipseMoments *= ellipseMass / 5.0;
        cout << SDUMP( ellipseMass ) << endl;
        cout << SDUMP( ellipseMoments ) << endl;
        cout << SDUMP( _mass ) << endl;
        cout << SDUMP( _inertiaMoments ) << endl;
        char tmp;
        cout << "Enter a character: ";
        cin >> tmp;
        cout << endl;
#endif
    } else {
        // Neglecting mass scaling for a moment, the moments of inertia for
        // the three axes of the matching ellipse (major, intermediate, minor)
        // are:
        //    (b^2 + c^2)
        //    (a^2 + c^2)
        //    (a^2 + b^2)
        //
        // Solve for a; a^2 = ( d_3 + d_2 - d_1 ) / 2
        a = _inertiaMoments[ 2 ] + _inertiaMoments[ 1 ] - _inertiaMoments[ 0 ];
        a /= 2;
        a = sqrt( a );

        // Solve for b; b^2 = d_3 - a^2
        b = sqrt( _inertiaMoments[ 2 ] - a * a );

        // Solve for c; c^2 = d_2 - a^2
        c = sqrt( _inertiaMoments[ 1 ] - a * a );

        // These quantities should be decreasing
        TRACE_ASSERT( a >= b && b >= c, "Invalid ellipse parameters" );

        if ( measure == VOLUME ) {
            // Make the ellipse volume and mesh volume match
            V = _mass / _density;

            volumeScaling = 3.0 * V / ( 4.0 * M_PI * a * b * c );
            volumeScaling = pow( volumeScaling, 1.0 / 3.0 );
        } else if ( measure == SURFACE_AREA ) {
#if 0
            A = _mesh.surfaceArea();
#endif
            A = _mesh.total_area();

            ellipseArea = pow( a, p ) * pow( b, p );
            ellipseArea += pow( a, p ) * pow( c, p );
            ellipseArea += pow( b, p ) * pow( c, p );
            ellipseArea /= 3.0;
            ellipseArea = pow( ellipseArea, 1.0 / p );
            ellipseArea *= 4.0 * M_PI;

            volumeScaling = A / ellipseArea;
            volumeScaling = pow( volumeScaling, 1.0 / 2.0 );
        }

        // Rescale ellipse axes
        a *= volumeScaling;
        b *= volumeScaling;
        c *= volumeScaling;

#if 0
        REAL ellipseMass = _density * 4.0 * M_PI * a * b * c / 3.0;
        Vector3d ellipseMoments( b * b + c * c, a * a + c * c, a * a + b * b );
        ellipseMoments *= ellipseMass / 5.0;
        cout << SDUMP( ellipseMass ) << endl;
        cout << SDUMP( ellipseMoments ) << endl;
        cout << SDUMP( _mass ) << endl;
        cout << SDUMP( _inertiaMoments ) << endl;
        char tmp;
        cout << "Enter a character: ";
        cin >> tmp;
        cout << endl;
#endif
    }

    // Figure out the scaling necessary to make the major axis have
    // unit length (ie. a = 0.5)
    _matchingScale = a / 0.5;

    // Figure out scaling of the other axes in the unscaled mesh
    _xScale = c / _matchingScale;
    _yScale = b / _matchingScale;

    _xScale /= 0.5;
    _yScale /= 0.5;

    TRACE_ASSERT( _xScale <= _yScale, "Invalid ellipse scales" );

    printf( "Generated fitting ellipse with axes (%e, %e, %e) and scales "
            "(%e, %e)\n", a, b, c, _xScale, _yScale );

#if 0
    // FIXME: debugging
    {
        MATRIX                   inertiaCopy( _inertia );
        MATRIX                   axisCopy( _inertiaAxes );
        VECTOR                   momentCopy( _inertiaMoments );

        inertiaCopy.write( "inertia.matrix" );
        axisCopy.write( "axes.matrix" );
        momentCopy.write( "moments.vector" );

        abort();
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Resets the fit ellipsoid for this mesh and determines a new scaling
//////////////////////////////////////////////////////////////////////
void RigidMesh::setEllipsoid( REAL xScale, REAL yScale, ObjectMeasure measure )
{
    REAL                       a, b, c;
    REAL                       V;
    REAL                       A;
    REAL                       ellipseArea;
    REAL                       volumeScaling;

    const REAL                 p = 1.6075;

    a = 1.0;
    b = yScale;
    c = xScale;

    if ( measure == VOLUME ) {
        // Make the ellipse volume and mesh volume match
        V = _mass / _density;

        volumeScaling = 3.0 * V / ( 4.0 * M_PI * a * b * c );
        volumeScaling = pow( volumeScaling, 1.0 / 3.0 );
    } else if ( measure == SURFACE_AREA ) {
#if 0
        A = _mesh.surfaceArea();
#endif
        A = _mesh.total_area();

        ellipseArea = pow( a, p ) * pow( b, p );
        ellipseArea += pow( a, p ) * pow( c, p );
        ellipseArea += pow( b, p ) * pow( c, p );
        ellipseArea /= 3.0;
        ellipseArea = pow( ellipseArea, 1.0 / p );
        ellipseArea *= 4.0 * M_PI;

        volumeScaling = A / ellipseArea;
        volumeScaling = pow( volumeScaling, 1.0 / 2.0 );
    }

    // Rescale ellipse axes
    a *= volumeScaling;
    b *= volumeScaling;
    c *= volumeScaling;

    // Figure out the scaling necessary to make the major axis have
    // unit length (ie. a = 0.5)
    _matchingScale = a / 0.5;

    // Figure out scaling of the other axes in the unscaled mesh
    _xScale = c / _matchingScale;
    _yScale = b / _matchingScale;

    _xScale /= 0.5;
    _yScale /= 0.5;

    TRACE_ASSERT( _xScale <= _yScale, "Invalid ellipse scales" );

    printf( "Adjusting fitting ellipse with axes (%e, %e, %e) and scales "
            "(%e, %e)\n", a, b, c, _xScale, _yScale );
}
