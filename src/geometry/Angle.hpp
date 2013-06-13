/*
 * =====================================================================================
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
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  Angle.hpp
 *
 *    Description:  functions related to angles in 3D space
 *
 *        Version:  1.0
 *        Created:  02/24/2011 03:03:21 PM
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef GEOMETRY_ANGLE_INC
#   define GEOMETRY_ANGLE_INC

#include <math.h>
#include "geometry/Point3.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

namespace Angle
{
    /*!
     * Check if the angle p1-p2-p3 is an obtuse angle
     *
     * NOTE: make sure the three points p1, p2, and p3 are distinct
     */
    template <typename T>
    inline bool is_obtuse(const Point3<T>& p1, const Point3<T>& p2, const Point3<T>& p3)
    {
        return (p3-p2).dotProduct(p1-p2) < 0.;
    }

    /*!
     * Compute the cot value of the angle specified by three points
     */
    template <typename T>
    inline double cot(const Point3<T>& p1, const Point3<T>& p2, const Point3<T>& p3)
    {
        const Vector3<T> u = p1 - p2;
        const Vector3<T> v = p3 - p2;
        const double c = u.dotProduct(v);
        return c / sqrt(u.lengthSqr()*v.lengthSqr() - c*c);
    }
}

#ifdef USE_NAMESPACE
}
#endif
#endif
