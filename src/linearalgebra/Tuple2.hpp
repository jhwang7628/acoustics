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
 *       Filename:  Tuple2.hpp
 *
 *        Version:  1.0
 *        Created:  11/29/11 15:01:12
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef TUPLE_2_INC
#   define TUPLE_2_INC

#include <assert.h>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class Tuple2
{
    public:
        typedef T   element;

        union
        {
            T x;
            T u;
        };

        union
        {
            T y;
            T v;
        };

        // ========== Constructors ===========
        Tuple2():x(0), y(0) { }
        Tuple2(T xv, T yv):x(xv), y(yv) { }
        Tuple2(const Tuple2<T>& src):x(src.x), y(src.y) {}

        template <typename FromT>
        Tuple2(const Tuple2<FromT>& src):x(static_cast<T>(src.x)),
                                         y(static_cast<T>(src.y))
        {}

        template <typename FromT>
        void set(FromT x, FromT y)
        {
            this->x = static_cast<T>(x);
            this->y = static_cast<T>(y);
        }

        Tuple2<T>& operator *= (const T rhs)
        {
            x *= rhs;
            y *= rhs;
            return *this;
        }

        Tuple2<T>& operator -= (const Tuple2<T>& rhs)
        {   
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }

        Tuple2<T>& operator += (const Tuple2<T>& rhs)
        {   
            x += rhs.x;
            y += rhs.y;
            return *this;
        }
};

typedef class Tuple2<int>     Tuple2i;
typedef class Tuple2<float>   Tuple2f;
typedef class Tuple2<double>  Tuple2d;

#ifdef USE_NAMESPACE
}
#endif

#endif
