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
 *       Filename:  Vector2.hpp
 *
 *        Version:  1.0
 *        Created:  11/29/11 15:07:00
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef VECTOR_2_INC
#   define VECTOR_2_INC

#include "Tuple2.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class Vector2 : public Tuple2<T>
{
    public:
        using Tuple2<T>::x;
        using Tuple2<T>::y;

        static const Vector2<T>  ZERO;

        Vector2() {}
        Vector2(T nx, T ny):Tuple2<T>(nx, ny) {}
        /*! Copy constructor. */
        Vector2(const Tuple2<T>& src):Tuple2<T>(src) {}
        /*! Copy casting constructor. */
        template <typename FromT>
        Vector2(const Tuple2<FromT>& src):Tuple2<T>(src) {}

        Vector2<T>& operator = (const Vector2<T>& rhs)
        {
            x = rhs.x;
            y = rhs.y;
            return *this;
        }

        inline Vector2<T> operator - (const Vector2<T>& rhs) const
        {
            return Vector2<T>(x-rhs.x, y-rhs.y);
        }

        inline T lengthSqr() const
        {   return x*x + y*y; }

        inline T length() const
        {   return sqrt(x*x + y*y); }

        inline void normalize()
        {
            T s = length();
            if ( s > 0 )
            {
                s = (T)1 / s;
                x *= s;
                y *= s;
            }
        }

        /*! Dot product of two vectors. */
        inline T dotProduct(const Vector2<T>& rhs) const 
        {
            return x * rhs.x + y * rhs.y;
        }

};

typedef class Vector2<float>    Vector2f;
typedef class Vector2<double>   Vector2d;
typedef class Vector2<int>      Vector2i;

template <typename T>
const Vector2<T> Vector2<T>::ZERO(0, 0);

#ifdef USE_NAMESPACE
}
#endif

#endif
