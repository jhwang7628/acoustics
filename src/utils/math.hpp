#ifndef DIFF_DEFINE
/******************************************************************************
 *  File: math.hpp
 *
 *  This file is part of isostuffer
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#endif /* ! DIFF_DEFINE */
/*
 * math.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef MATH_UTILS_H
#   define MATH_UTILS_H

#ifdef DIFF_DEFINE
#include <math.h>
#include "utils/macros.h"

#endif /* DIFF_DEFINE */
#ifndef M_PI
#   define M_PI        3.14159265358979323846264338327950288   /* pi */
#endif

#ifndef M_1_PI
#   define M_1_PI      0.318309886183790671537767526745028724  /* 1/pi */
#endif

template <typename T> inline T M_NEG(T a)
{   return -a; }

template <typename T> inline T M_MAX(const T a, const T b)
{   return a > b ? a : b; }

template <typename T> inline T M_MIN(const T a, const T b)
{   return a < b ? a : b; }

template <typename T> inline T M_DEG2RAD(T x)
{   return x * M_PI / 180.0; }

template <typename T> inline T M_RAD2DEG(T x)
{   return x * 180.0 * M_1_PI; }

template <typename T> inline T M_ABS(T a)
{   return a > (T)0 ? a : -a; }

template <typename T> inline T M_TRI(T a)
{   return a * a * a; }

template <typename T> inline T M_SQR(T x)
{   return x*x; }

template <typename T> inline int M_SIGN(T x)
{   return x > (T)0 ? 1 : (x < (T)0 ? -1 : 0); }

inline int gcd(int a, int b)
{
    for(int c;b;c=a,a=b,b=c%b) ;
    return a;
}

inline int lcm(int a, int b)
{
    return a/gcd(a,b)*b;
}

template <typename T>
inline T clamp(T a, T minv, T maxv)
{
    return a <= minv ? minv : (a > maxv ? maxv : a);
}

template <typename T>
inline T gaussian(T v, T mu, T d)
{
    return exp(-M_SQR((v - mu) / d));
#ifdef DIFF_DEFINE
}

template <typename T>
inline T cube_root(T x);

template <>
inline float cube_root<float>(float x)
{
    float fr, r;
    int   ex, shx;
    bool neg = false;

    if ( unlikely(x==0.) ) return 0.;
    /* Argument reduction */
    if ( x < 0.F ) 
    {
        x *= -1.F;
        neg = true;
    }
    fr = frexp(x, &ex);

    /* separate into mantissa and exponent */
    shx = ex % 3;
    if (shx > 0) shx -= 3; /* compute shx such that (ex - shx) is divisible by 3 */
    ex = (ex - shx) / 3;
    /* exponent of cube root */
    fr = ldexp(fr, shx);
    /* 0.125 <= fr < 1.0 */

#ifdef NO_ITERATE_CUBE_ROOT
    /* Use quartic rational polynomial with error < 2^(-24) */
    fr = ((((45.2548339756803022511987494 * fr +
          192.2798368355061050458134625) * fr +
          119.1654824285581628956914143) * fr +
          13.43250139086239872172837314) * fr +
          0.1636161226585754240958355063) /
         ((((14.80884093219134573786480845 * fr +
          151.9714051044435648658557668) * fr +
          168.5254414101568283957668343) * fr +
          33.9905941350215598754191872) * fr +
          1.0);
    r = ldexp(fr, ex); /* 24 bits of precision */
#else 
    /* Compute seed with a quadratic qpproximation */
    fr = (-0.46946116F * fr + 1.072302F) * fr + 0.3812513F;     /* 0.5<=fr<1 */
    r = ldexp(fr, ex);
    /* 6 bits of precision */
    /* Newton-Raphson iterations */
    r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r*r);    /* 12 bits */
    r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r*r);    /* 24 bits */
#endif
    return neg ? -r : r;
}

template <>
inline double cube_root<double>(double x)
{
    double fr, r;
    int    ex, shx;
    bool   neg = false;

    if ( unlikely(x==0.) ) return 0.;
    /* Argument reduction */
    if ( x < 0. ) 
    {
        x *= -1.;
        neg = true;
    }
    fr = frexp(x, &ex);

    /* separate into mantissa and exponent */
    shx = ex % 3;
    if (shx > 0) shx -= 3; /* compute shx such that (ex - shx) is divisible by 3 */
    ex = (ex - shx) / 3;
    /* exponent of cube root */
    fr = ldexp(fr, shx);
    /* 0.125 <= fr < 1.0 */

    /* Use quartic rational polynomial with error < 2^(-24) */
    fr = ((((45.2548339756803022511987494 * fr +
          192.2798368355061050458134625) * fr +
          119.1654824285581628956914143) * fr +
          13.43250139086239872172837314) * fr +
          0.1636161226585754240958355063) /
         ((((14.80884093219134573786480845 * fr +
          151.9714051044435648658557668) * fr +
          168.5254414101568283957668343) * fr +
          33.9905941350215598754191872) * fr +
          1.0);
    r = ldexp(fr, ex); /* 24 bits of precision */
    /* Newton-Raphson iterations */
    r = (2.0/3.0) * r + (1.0/3.0) * x / (r*r);    /* 48 bits */
    r = (2.0/3.0) * r + (1.0/3.0) * x / (r*r);    /* 96 bits */
    return neg ? -r : r;
#endif /* DIFF_DEFINE */
}

#endif
