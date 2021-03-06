/******************************************************************************
 *  File: Vector3.hpp
 *  A vector in 3D space
 *  Copyright (c) 2007 by Changxi Zheng
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
#ifndef LINEARALGEBRA_VECTOR_3_HPP
#   define LINEARALGEBRA_VECTOR_3_HPP

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <assert.h>
#include "Tuple3.hpp"
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif
    
//! Class for three dimensional vector.
template <typename T> 
class Vector3 : public Tuple3<T>
{
    public:
        using Tuple3<T>::x;
        using Tuple3<T>::y;
        using Tuple3<T>::z;

        static const Vector3<T>  ZERO;

        /*! Creates and sets to (0,0,0) */
        Vector3() {}
        /*! Creates and sets to (x,y,z) */
        Vector3(T nx, T ny, T nz):Tuple3<T>(nx, ny, nz) {}
        /*! Copy constructor. */
        Vector3(const Tuple3<T>& src):Tuple3<T>(src) {}
        /*! Copy casting constructor. */
        template <typename FromT>
        Vector3(const Tuple3<FromT>& src):Tuple3<T>(src) {}
       
        Vector3<T>& operator = (const Tuple3<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
        }

        /*! Copy casting operator. */
        template <typename FromT>
        Vector3<T>& operator = (const Tuple3<FromT>& rhs)
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            z = static_cast<T>(rhs.z);
            return *this;
        }

        Vector3<T>& operator = (const T rhs)
        {
            x = y = z = rhs;
            return *this;
        }
        /*
        template <typename FromT>
        Vector3<T>& operator = (const FromT rhs)
        {
            x = y = z = static_cast<T>(rhs);
            return *this;
        }
        */

        Vector3<T> operator + (const Vector3<T>& rhs) const
        {
            return Vector3<T>(x + rhs.x, y + rhs.y, z + rhs.z);
        }

        Vector3<T> operator + (const T &rhs) const
        {
            return Vector3<T>(x+rhs, y+rhs, z+rhs); 
        }

        Vector3<T> operator - (const Vector3<T>& rhs) const 
        {
            return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
        }

        Vector3<T> operator - () const
        {
            return Vector3<T>(-x, -y, -z);
        }

        Vector3<T> operator - (const T &rhs) const
        {
            return Vector3<T>(x-rhs, y-rhs, z-rhs); 
        }

        /*! Multiplication operator */
        Vector3<T> operator * (const T& rhs) const
        {
            return Vector3<T>(x*rhs, y*rhs, z*rhs);
        }

        Vector3<T> operator * (const Vector3<T>& rhs) const 
        {
            return Vector3<T>(x * rhs.x, y * rhs.y, z * rhs.z);
        }

        Vector3<T> operator / (const T& rhs) const
        {
            return Vector3<T>(x/rhs, y/rhs, z/rhs);
        }

        T &operator () (const int &ind)
        {
            if (ind == 0) 
                return x; 
            else if (ind == 1) 
                return y;
            else if (ind == 2) 
                return z; 
            else 
                throw std::runtime_error("**ERROR** vector access out of bounds");
        }

        T operator () (const int &ind) const
        {
            T returnValue;
            if (ind == 0) 
                returnValue = x; 
            else if (ind == 1) 
                returnValue = y;
            else if (ind == 2) 
                returnValue = z; 
            else 
                throw std::runtime_error("**ERROR** vector access out of bounds");

            return returnValue;
        }

        /*! Dot product of two vectors. */
        T dotProduct(const Vector3<T>& rhs) const 
        {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        /*! Cross product opertor */    
        Vector3<T> crossProduct(const Vector3<T>& rhs) const // $$TESTED
        {
            return Vector3<T>(y * rhs.z - rhs.y * z, 
                              z * rhs.x - rhs.z * x, 
                              x * rhs.y - rhs.x * y);
        }
       
        /*! Get lenght of vector.*/
        T length() const 
        {
            return (T)std::sqrt(x * x + y * y + z * z);
        }
       
        /*!
         * Return square of length.
         * @return length ^ 2
         * @note This method is faster then length(). For comparison
         * of length of two vector can be used just this value, instead
         * of computionaly more expensive length() method.
         */
        T lengthSqr() const 
        {
            return x * x + y * y + z * z;
        }

        /*! 
         * Project the vector to the plane that has normal v_n
         */
        void ApplyProjectionInplace(const Vector3<T> &v_n)
        {
            const T vi_cos_t = this->dotProduct(v_n) / v_n.norm(); 
            x = x - v_n.x * vi_cos_t; 
            y = y - v_n.y * vi_cos_t; 
            z = z - v_n.z * vi_cos_t; 
        }
       
        /*! Normalize vector */
        void normalize() 
        {
            T s = length();
            if ( s > 0 )
            {
#ifndef DIFF_DEFINE
                x /= s;
                y /= s;
                z /= s;
#else /* DIFF_DEFINE */
                s = (T)1 / s;
                x *= s;
                y *= s;
                z *= s;
#endif /* DIFF_DEFINE */
            }
        }

        /*! Out-of-place normalize */
        Vector3<T> normalized() const
        {
            Vector3<T> newVector(x,y,z); 
            newVector.normalize(); 
            return newVector; 
        }

        /*!
         * Normalize vector and return its original length
         */
        T normalize2()
        {
            T s = length();
            if ( s > 0 )
            {
#ifndef DIFF_DEFINE
                x /= s;
                y /= s;
                z /= s;
#else /* DIFF_DEFINE */
                s = (T)1 / s;
                x *= s;
                y *= s;
                z *= s;
#endif /* DIFF_DEFINE */
            }
            return s;
        }
       
        /*!
         * Rotate vector around three axis.
         * @param ax Angle (in degrees) to be rotated around X-axis.
         * @param ay Angle (in degrees) to be rotated around Y-axis.
         * @param az Angle (in degrees) to be rotated around Z-axis.
         */
        void rotate(T ax, T ay, T az) 
        {
            T a = (T)cos(DEG2RAD(ax));
            T b = (T)sin(DEG2RAD(ax));
            T c = (T)cos(DEG2RAD(ay));
            T d = (T)sin(DEG2RAD(ay));
            T e = (T)cos(DEG2RAD(az));
            T f = (T)sin(DEG2RAD(az));
            T nx = c * e * x - c * f * y + d * z;
            T ny = (a * f + b * d * e) * x + (a * e - b * d * f) * y - b * c * z;
            T nz = (b * f - a * d * e) * x + (a * d * f + b * e) * y + a * c * z;
            x = nx; y = ny; z = nz;
        }

        void set(T x, T y, T z) 
        {
            this->x = x; 
            this->y = y; 
            this->z = z; 
        }
        
}; // end of Vector3

typedef class Vector3<float>    Vector3f;
typedef class Vector3<double>   Vector3d;
typedef class Vector3<int>      Vector3i;

template <typename T>
const Vector3<T> Vector3<T>::ZERO(0, 0, 0);

inline Vector3d str2vec3d( const std::string& str )
{
	Vector3d v;
	std::istringstream ss( str );
	ss >> v[0] >> v[1] >> v[2];
	return v;
}

inline Vector3f str2vec3f( const std::string& str )
{
	Vector3f v;
	std::istringstream ss( str );
	ss >> v[0] >> v[1] >> v[2];
	return v;
}

#ifdef USE_NAMESPACE
}
#endif

#endif
