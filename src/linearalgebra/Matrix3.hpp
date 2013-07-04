/******************************************************************************
 *  File: Matrix3.hpp
 *  A 3x3 Matrix
 *  Copyright (c) 2007 by Changxi Zheng
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
#ifndef LINEARALGEBRA_MATRIX_3_HPP
#   define LINEARALGEBRA_MATRIX_3_HPP

#include <stdio.h>
#include <iostream>
#include <assert.h>

#include <stdlib.h>

#include "Vector3.hpp"
#include "geometry/Point3.hpp"
#include "utils/math.hpp"
#ifdef DIFF_DEFINE
#include "generic/precision_type.hpp"
#endif /* DIFF_DEFINE */

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class Matrix3 
{
    public:
        Vector3<T>  cols[3];

    public:
        static const Matrix3<T>     I;

        /* ========== constructors ========== */
        Matrix3() { }
        Matrix3(T m11, T m12, T m13,
                T m21, T m22, T m23,
                T m31, T m32, T m33)
        {
            cols[0].set(m11, m21, m31);
            cols[1].set(m12, m22, m32);
            cols[2].set(m13, m23, m33);
        }
        Matrix3(const Tuple3<T>& c1, 
                const Tuple3<T>& c2,
                const Tuple3<T>& c3)
        {
            cols[0] = c1;
            cols[1] = c2;
            cols[2] = c3;
        }

        void zero()
        {
            cols[0].zero();
            cols[1].zero();
            cols[2].zero();
        }

        void set(const Tuple3<T>& c1,
                 const Tuple3<T>& c2,
                 const Tuple3<T>& c3)
        {
            cols[0] = c1;
            cols[1] = c2;
            cols[2] = c3;
        }		

        void set(T m11, T m12, T m13,
                 T m21, T m22, T m23,
                 T m31, T m32, T m33)
        {
            cols[0].set(m11, m21, m31);
            cols[1].set(m12, m22, m32);
            cols[2].set(m13, m23, m33);
        }

        /* ========== copy constructor ========= */
        template <typename FromT>
        Matrix3<T>& operator = (const Matrix3<FromT>& rhs)
        {
            cols[0] = rhs.cols[0];
            cols[1] = rhs.cols[1];
            cols[2] = rhs.cols[2];
            return *this;
        }

        /*
         * The double product (double contraction) of A and B
         *  A : B = tr(A^T * B) = \sum_{i,j=1}^3 A_{ij}B_{ij}
         */
        T colon(const Matrix3<T>& rhs) const
        {
            return cols[0][0] * rhs.cols[0][0] + 
                   cols[0][1] * rhs.cols[0][1] +
                   cols[0][2] * rhs.cols[0][2] +
                   cols[1][0] * rhs.cols[1][0] +
                   cols[1][1] * rhs.cols[1][1] +
                   cols[1][2] * rhs.cols[1][2] +
                   cols[2][0] * rhs.cols[2][0] +
                   cols[2][1] * rhs.cols[2][1] +
                   cols[2][2] * rhs.cols[2][2]; 
        }

        T trace() const
        {
            return cols[0][0] + cols[1][1] + cols[2][2];
        }

        inline T operator () (int i, int j) const 
        { 
            assert(i < 3 && j < 3 && i >= 0 && j >= 0);
            return cols[j][i]; 
        }

        inline T& operator () (int i, int j) 
        { 
            assert(i < 3 && j < 3 && i >= 0 && j >= 0);
            return cols[j][i]; 
        }

        //! right multiply a vector
#if 0
#ifndef DIFF_DEFINE
        Vector3<T> operator * (const Vector3<T>& rhs) const
        {
            return cols[0]*rhs[0] + cols[1]*rhs[1] + cols[2]*rhs[2];
        }

        Point3<T> operator * (const Point3<T>& rhs) const
        {
            return cols[0]*rhs[0] + cols[1]*rhs[1] + cols[2]*rhs[2];
        }
#endif /* DIFF_DEFINE */

        Tuple3<T> operator * (const Tuple3<T>& rhs) const
        {
#ifndef DIFF_DEFINE
            return cols[0]*rhs[0] + cols[1]*rhs[1] + cols[2]*rhs[2];
#else /* DIFF_DEFINE */
            return cols[0]*rhs.x + cols[1]*rhs.y + cols[2]*rhs.z;
#endif /* DIFF_DEFINE */
        }
#endif
        Vector3<T> operator * (const Vector3<T>& rhs) const
        {
            return cols[0]*rhs[0] + cols[1]*rhs[1] + cols[2]*rhs[2];
        }

        Point3<T> operator * (const Point3<T>& rhs) const
        {
            return cols[0]*rhs[0] + cols[1]*rhs[1] + cols[2]*rhs[2];
        }

#ifdef DIFF_DEFINE
        Tuple3<T> multiply(const Tuple3<T>* rhs) const
        {   return cols[0]*rhs->x + cols[1]*rhs->y + cols[2]*rhs->z; }
        Tuple3<T> multiply(T x, T y, T z) const
        {   return cols[0]*x + cols[1]*y + cols[2]*z; }

#endif /* DIFF_DEFINE */
        Matrix3<T> operator * (const Matrix3<T>& rhs) const
        {
            return Matrix3<T>(
                    cols[0]*rhs.cols[0][0] + cols[1]*rhs.cols[0][1] + cols[2]*rhs.cols[0][2],
                    cols[0]*rhs.cols[1][0] + cols[1]*rhs.cols[1][1] + cols[2]*rhs.cols[1][2],
                    cols[0]*rhs.cols[2][0] + cols[1]*rhs.cols[2][1] + cols[2]*rhs.cols[2][2]);
        }

        Matrix3<T> operator * (T rhs) const
        {
            return Matrix3<T>(cols[0]*rhs, cols[1]*rhs, cols[2]*rhs);
        }

        friend Matrix3<T> operator * (T lhs, const Matrix3<T>& rhs)
        {
            return Matrix3<T>(rhs.cols[0] * lhs,
                              rhs.cols[1] * lhs,
                              rhs.cols[2] * lhs);
        }

        Matrix3<T> operator + (const Matrix3<T>& rhs) const
        {
            return Matrix3<T>(cols[0] + rhs.cols[0],
                              cols[1] + rhs.cols[1],
                              cols[2] + rhs.cols[2]);
        }

        Matrix3<T> operator - (const Matrix3<T>& rhs) const
        {
            return Matrix3<T>(cols[0] - rhs.cols[0],
                              cols[1] - rhs.cols[1],
                              cols[2] - rhs.cols[2]);
        }

        Matrix3<T>& operator += (const Matrix3<T>& rhs) 
        {
            cols[0] += rhs.cols[0];
            cols[1] += rhs.cols[1];
            cols[2] += rhs.cols[2];
            return *this;
        }

#ifdef DIFF_DEFINE
        Matrix3<T>& operator -= (const Matrix3<T>& rhs) 
        {
            cols[0] -= rhs.cols[0];
            cols[1] -= rhs.cols[1];
            cols[2] -= rhs.cols[2];
            return *this;
        }

#endif /* DIFF_DEFINE */
        Matrix3<T>& operator *= (const T rhs)
        {
            cols[0] *= rhs;
            cols[1] *= rhs;
            cols[2] *= rhs;
            return *this;
        }

        Matrix3<T>& operator /= (const T rhs)
        {
            T s = (T)1 / rhs;
            cols[0] *= s;
            cols[1] *= s;
            cols[2] *= s;
            return *this;
        }

        inline void scaleAdd(T s, const Matrix3<T>& rhs)
        {
            cols[0].scaleAdd(s, rhs.cols[0]);
            cols[1].scaleAdd(s, rhs.cols[1]);
            cols[2].scaleAdd(s, rhs.cols[2]);
        }

        T det() const
        {
            return cols[1].crossProduct(cols[2]).dotProduct(cols[0]);
        }

        Matrix3<T> inverse() const
        {
            Vector3<T> v0 = cols[1].crossProduct(cols[2]);
            T det = v0.dotProduct(cols[0]);
#ifndef DIFF_DEFINE
            if ( M_ABS(det) < 1E-28 )
#else /* DIFF_DEFINE */
            if ( M_ABS(det) < PrecisionType<T>::MA_EPS )
#endif /* DIFF_DEFINE */
            {
#ifndef DIFF_DEFINE
                fprintf(stderr, "ERROR: Invert matrix(1) with det<1E-28 %.30lf\n", det);
#else /* DIFF_DEFINE */
                fprintf(stderr, "ERROR: Invert matrix(1) with det close to 0 %.20lf\n", det);
                exit(1);
#endif /* DIFF_DEFINE */
                return Matrix3<T>::I;
            }
            Vector3<T> v1 = cols[2].crossProduct(cols[0]);
            Vector3<T> v2 = cols[0].crossProduct(cols[1]);

            Matrix3<T> ret(v0[0], v0[1], v0[2],
                           v1[0], v1[1], v1[2],
                           v2[0], v2[1], v2[2]);
            ret /= det;
            return ret;
        }

        int inverse(Matrix3<T>& out) const
        {
            Vector3<T> v0 = cols[1].crossProduct(cols[2]);
            T det = v0.dotProduct(cols[0]);
#ifndef DIFF_DEFINE
            if ( M_ABS(det) < 1E-20 )
#else /* DIFF_DEFINE */
            if ( M_ABS(det) < PrecisionType<T>::MA_EPS )
#endif /* DIFF_DEFINE */
            {
#ifndef DIFF_DEFINE
                fprintf(stderr, "ERROR: Invert matrix(2) with det<1E-20 %.30lf\n", det);
#endif /* ! DIFF_DEFINE */
                out = Matrix3<T>::I;
                return -1;
            }
            Vector3<T> v1 = cols[2].crossProduct(cols[0]);
            Vector3<T> v2 = cols[0].crossProduct(cols[1]);

            out.set(v0[0], v0[1], v0[2],
                    v1[0], v1[1], v1[2],
                    v2[0], v2[1], v2[2]);
            out /= det;
            return 0;
        }

        Matrix3<T> transpose() const
        {
            return Matrix3(cols[0][0], cols[0][1], cols[0][2],
                           cols[1][0], cols[1][1], cols[1][2],
                           cols[2][0], cols[2][1], cols[2][2]);
        }

        // Read the contents of this matrix from the given file.
        // If rowMajor == true, then we assume that the contents of the
        // file are stored in row major, rather than column major format
        bool read( const char *fileName, bool rowMajor = false )
        {
            FILE* file;
            file = fopen(fileName, "rb");

            if ( !file )
            {
                std::cerr << "Error opening " << fileName << std::endl;
                return false;
            }

            fread( (void *)&cols[ 0 ], sizeof( Vector3<T> ), 1, file );
            fread( (void *)&cols[ 1 ], sizeof( Vector3<T> ), 1, file );
            fread( (void *)&cols[ 2 ], sizeof( Vector3<T> ), 1, file );

            fclose( file );

            // Transpose if necessary
            if ( rowMajor ) {
                *this = this->transpose();
            }

            return true;
        }

        // Write contents of this matrix to the given file.
        // If rowMajor == true, then we will write this data in row major
        // format.
        bool write( const char *fileName, bool rowMajor = false ) const
        {
            if ( rowMajor ) {
                return transpose().write( fileName, false );
            } else {

                FILE* file;
                file = fopen(fileName, "wb");

                if ( !file )
                {
                    std::cerr << "Error opening " << fileName << std::endl;
                    return false;
                }

                fwrite( (void *)&cols[ 0 ], sizeof( Vector3<T> ), 1, file );
                fwrite( (void *)&cols[ 1 ], sizeof( Vector3<T> ), 1, file );
                fwrite( (void *)&cols[ 2 ], sizeof( Vector3<T> ), 1, file );

                fclose( file );

                return true;
            }
        }

        // Writes the contents

        friend std::ostream& operator<<(std::ostream& lhs, const Matrix3<T>& rhs) 
        {
            lhs << rhs.cols[0][0] << ' ' << rhs.cols[1][0] << ' ' << rhs.cols[2][0] << std::endl
                << rhs.cols[0][1] << ' ' << rhs.cols[1][1] << ' ' << rhs.cols[2][1] << std::endl
                << rhs.cols[0][2] << ' ' << rhs.cols[1][2] << ' ' << rhs.cols[2][2];
            return lhs;
        }

        static Matrix3<T> identity(T v)
        {
            return Matrix3<T>(v, 0, 0, 0, v, 0, 0, 0, v);
        }
        static Matrix3<T> diagonal(const Tuple3<T>& v)
        {
            return Matrix3<T>(v.x, 0, 0, 0, v.y, 0, 0, 0, v.z);
        }
        static Matrix3<T> crossProductMatrix(const Tuple3<T>& v)
        {
            return Matrix3<T>( 0.0, -v[2], v[1],
                               v[2], 0.0, -v[0],
                               -v[1], v[0], 0.0 );
        }

};

template <typename T>
const Matrix3<T> Matrix3<T>::I(1, 0, 0, 0, 1, 0, 0, 0, 1);

typedef Matrix3<double> Matrix3d;
typedef Matrix3<float>  Matrix3f;

#ifdef USE_NAMESPACE
}
#endif

#endif
