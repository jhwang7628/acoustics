//////////////////////////////////////////////////////////////////////
// MathUtil.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <TYPES.h>

#include <linearalgebra/Vector3.hpp>

#include "MERSENNETWISTER.h"

#include <math.h>
#include <vector>

//////////////////////////////////////////////////////////////////////
// Some random utility functions that I needed a place for
//////////////////////////////////////////////////////////////////////
class MathUtil {
    public:
        //////////////////////////////////////////////////////////////////////
        // Computes the correct mod function (ie. handles negative numbers)
        // (eg. i % j)
        //////////////////////////////////////////////////////////////////////
        static inline int mod_fixed( int i, int j )
        {
            int mod_result = i % j;

            if ( mod_result < 0 )
            {
                mod_result += j;
            }

            return mod_result;
        }

        static inline void swap( REAL &x, REAL &y )
        {
            REAL                   tmp = x;

            x = y;
            y = tmp;
        }

        //////////////////////////////////////////////////////////////////////
        // Computes a uniformly randomly distributed point on a sphere
        // of the given radius, centered at the origin
        //////////////////////////////////////////////////////////////////////
        static inline Vector3d UniformRandomSpherePoint( MERSENNETWISTER &generator, REAL r )
        {
#if 0
            REAL                       z = generator.rand( 2.0 * r ) - r;
            REAL                       t = generator.rand( 2.0 * M_PI );
            REAL                       rDisc = sqrt( r * r - z * z );

            return Vector3d( rDisc * cos( t ), rDisc * sin( t ), z );
#endif
            REAL                       u = generator.rand( 1.0 );
            REAL                       v = generator.rand( 1.0 );

            REAL                       theta = 2.0 * M_PI * u;
            REAL                       phi = acos( 2.0 * v - 1 );

            return Vector3d( r * sin( phi ) * cos( theta ),
                             r * sin( phi ) * sin( theta ),
                             r * cos( phi ) );
        }

        //////////////////////////////////////////////////////////////////////
        // Generates a set of points on a sphere with the given center and
        // radius, uniformly sampled in angle space (but not uniformly sampled
        // on the sphere itself)
        //////////////////////////////////////////////////////////////////////
        static int GenerateSpherePoints( const Vector3d &center, REAL radius,
                                         int resolution,
                                         Vector3Array &points );

        // Generates a set of points on the positive octant of a sphere.
        static int GenerateOctantPoints( const Vector3d &center, REAL radius,
                                         int resolution,
                                         Vector3Array &points );

        // Converts cartesian coordinates to spherical coordinates
        // with respect to the origin
        static Vector3d spherical_coordinates( const Vector3d &x );

};

#endif
