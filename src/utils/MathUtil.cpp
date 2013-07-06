//////////////////////////////////////////////////////////////////////
// MathUtil.cpp
//////////////////////////////////////////////////////////////////////

#include "MathUtil.h"

//////////////////////////////////////////////////////////////////////
// Generates a set of points on a sphere with the given center and
// radius, uniformly sampled in angle space (but not uniformly sampled
// on the sphere itself)
//////////////////////////////////////////////////////////////////////
int MathUtil::GenerateSpherePoints( const Vector3d &center, REAL radius,
                                    int resolution,
                                    Vector3Array &points )
{
    REAL theta, phi;
    REAL cos_theta, sin_theta;

    // The azimuth resolution is 2 * resolution, so the number of points
    // on the sphere is given by (resolution) * (2 * resolution) + 2
    // where the extra 2 points come from the north and south poles.
    int azimuthResolution = 2 * resolution;
    int sphereResolution = azimuthResolution * resolution + 2;

    REAL azimuthDivision = 2.0 * M_PI / (REAL)azimuthResolution;
    REAL zenithDivision = M_PI / (REAL)( resolution + 1 );

    // Add the north pole
    points.push_back( center + Vector3d( 0.0, 0.0, radius ) );

    for (int t = 1; t <= resolution; t++)
    {
        theta = zenithDivision * (REAL)t;
        cos_theta = cos( theta );
        sin_theta = sin( theta );

        for (int p = 0; p < azimuthResolution; p++)
        {
            phi = azimuthDivision * (REAL)p;

            Vector3d offset( radius * sin_theta * cos( phi ),
                             radius * sin_theta * sin( phi ),
                             radius * cos_theta );

            points.push_back( center + offset );
        }
    }

    // Add the south pole
    points.push_back( center + Vector3d( 0.0, 0.0, -radius ) );

    return sphereResolution;
}

//////////////////////////////////////////////////////////////////////
// Generates a set of points on the positive octant of a sphere.
//////////////////////////////////////////////////////////////////////
int MathUtil::GenerateOctantPoints( const Vector3d &center, REAL radius,
                                    int resolution,
                                    Vector3Array &points )
{
    REAL theta, phi;
    REAL cos_theta, sin_theta;

    int azimuthResolution = resolution + 1;
    int sphereResolution = azimuthResolution * resolution + 1;

    REAL azimuthDivision = M_PI / 2.0 / (REAL)( azimuthResolution - 1 );
    REAL zenithDivision = M_PI / 2.0 / (REAL)( resolution );
#if 0
    int azimuthResolution = (int)ceil( (REAL)resolution / 2.0 );
    int sphereResolution = azimuthResolution * resolution + 1;

    REAL azimuthDivision = M_PI / 2.0 / (REAL)( azimuthResolution - 1 );
    REAL zenithDivision = M_PI / 2.0 / (REAL)( resolution );
#endif

    // Add the north pole
    points.push_back( center + Vector3d( 0.0, 0.0, radius ) );

    for (int t = 1; t <= resolution; t++)
    {
        theta = zenithDivision * (REAL)t;
        cos_theta = cos( theta );
        sin_theta = sin( theta );

        for (int p = 0; p < azimuthResolution; p++)
        {
            phi = azimuthDivision * (REAL)p;

            Vector3d offset( radius * sin_theta * cos( phi ),
                    radius * sin_theta * sin( phi ),
                    radius * cos_theta );

            points.push_back( center + offset );
        }
    }

    return sphereResolution;
}

//////////////////////////////////////////////////////////////////////
// Converts cartesian coordinates to spherical coordinates
// with respect to the origin
//////////////////////////////////////////////////////////////////////
Vector3d MathUtil::spherical_coordinates( const Vector3d &x )
{
    Vector3d r;

    r[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    r[0] = sqrt( r[0] );

    r[1] = acos( x[2] / r[0] );

    r[2] = atan2( x[1], x[0] );

    return r;
}
