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

#ifndef M_PI
#   define M_PI        3.14159265358979323846264338327950288   /* pi */
#endif

#ifndef M_1_PI
#   define M_1_PI      0.318309886183790671537767526745028724  /* 1/pi */
#endif

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

        template <typename T>
        static inline T M_NEG(T a)
        {   return -a; }

        template <typename T>
        static inline T M_MAX(const T a, const T b)
        {   return a > b ? a : b; }

        template <typename T>
        static inline T M_MIN(const T a, const T b)
        {   return a < b ? a : b; }

        template <typename T>
        static inline T M_DEG2RAD(T x)
        {   return x * M_PI / 180.0; }

        template <typename T>
        static inline T M_RAD2DEG(T x)
        {   return x * 180.0 * M_1_PI; }

        template <typename T>
        static inline T M_ABS(T a)
        {   return a > (T)0 ? a : -a; }

        template <typename T>
        static inline T M_TRI(T a)
        {   return a * a * a; }

        template <typename T>
        static inline T M_SQR(T x)
        {   return x*x; }

        template <typename T>
        static inline int M_SIGN(T x)
        {   return x > (T)0 ? 1 : (x < (T)0 ? -1 : 0); }

        static inline int gcd(int a, int b)
        {
                for(int c;b;c=a,a=b,b=c%b) ;
                    return a;
        }

        static inline int lcm(int a, int b)
        {
                return a/gcd(a,b)*b;
        }

        template <typename T>
        static inline T clamp(T a, T minv, T maxv)
        {
                return a <= minv ? minv : (a > maxv ? maxv : a);
        }

        template <typename T>
        static inline T gaussian(T v, T mu, T d)
        {
                return exp(-M_SQR((v - mu) / d));
        }

        template <typename T>
        static inline T gaussian_normalized(T v, T mu, T d)
        {
                return exp(-M_SQR((v - mu) / d)) / sqrt(d);
        }

};

#endif
