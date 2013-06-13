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
 *       Filename:  OuterProd.hpp
 *
 *        Version:  1.0
 *        Created:  09/16/11 11:57:20
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef OUTER_PRODUCT_INC
#   define OUTER_PRODUCT_INC

#include "Matrix3.hpp"

/*! out = out + a*b' */
template <typename T>
inline void outer_product_add(const Vector3<T>& a, const Vector3<T>& b, Matrix3<T>& out)
{
    out.cols[0].scaleAdd(b.x, a);
    out.cols[1].scaleAdd(b.y, a);
    out.cols[2].scaleAdd(b.z, a);
}

template <typename T>
inline void outer_product(const Vector3<T>& a, const Vector3<T>& b, Matrix3<T>& out)
{
    out.cols[0] = a*b.x;
    out.cols[1] = a*b.y;
    out.cols[2] = a*b.z;
}

#endif

