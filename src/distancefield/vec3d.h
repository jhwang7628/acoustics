#ifndef _MINIVEC3D_H_
#define _MINIVEC3D_H_

#include <math.h>
#include <ostream>
#include "vec_types.h"
                                                                                                                                                             
class Vec3d {

public:

	inline Vec3d() {}
	inline Vec3d(double x_g, double y_g, double z_g) {elt[0]=x_g; elt[1]=y_g; elt[2]=z_g;}
	inline Vec3d(ZeroOrOne k);

  inline Vec3d & operator=(const Vec3d & source);
  inline bool operator==(const Vec3d & vec2);
	inline Vec3d & operator=(ZeroOrOne k);

	inline Vec3d operator+ (const Vec3d & vec2);
	inline Vec3d & operator+= (const Vec3d & vec2);

	inline Vec3d operator- (const Vec3d & vec2);
	inline Vec3d & operator-= (const Vec3d & vec2);

	inline Vec3d operator* (double scalar) const;
	inline Vec3d & operator*= (double scalar);

	inline Vec3d operator/ (double scalar) const;
	inline Vec3d & operator/= (double scalar);

	friend inline Vec3d operator* (double scalar, const Vec3d & vec2);
	friend inline Vec3d operator/ (double scalar, const Vec3d & vec2);

	friend inline double dot(const Vec3d & vec1, const Vec3d & vec2);

	friend inline Vec3d cross(const Vec3d & vec1, const Vec3d & vec2);

  friend inline Vec3d norm(const Vec3d & vec1);

	friend class Mat3d;

	inline double & operator[] (int index);
	inline const double & operator[] (int index) const;

protected:
  double elt[3];

};


inline Vec3d::Vec3d(ZeroOrOne k)
{
  elt[0] = k;
  elt[1] = k;
  elt[2] = k;
}

inline Vec3d & Vec3d::operator=(const Vec3d & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];
  elt[2] = source.elt[2];

  return *this;
}

inline bool Vec3d::operator==(const Vec3d & vec2)
{
  return ((elt[0] == vec2[0]) &&
          (elt[1] == vec2[1]) &&
          (elt[2] == vec2[2]));
}

inline Vec3d & Vec3d::operator=(ZeroOrOne k)
{
  elt[0] = k;
  elt[1] = k;
  elt[2] = k;

  return *this;
}


inline Vec3d operator* (double scalar, const Vec3d & vec2)
{
  Vec3d result = vec2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;
  result.elt[2] *= scalar;

  return result;
}

inline Vec3d operator/ (double scalar, const Vec3d & vec2)
{
  Vec3d result = vec2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;
  result.elt[2] /= scalar;

  return result;
}


inline Vec3d Vec3d::operator+ (const Vec3d & vec2)
{
  Vec3d sum = *this;
  sum.elt[0] += vec2.elt[0];
  sum.elt[1] += vec2.elt[1];
  sum.elt[2] += vec2.elt[2];

  return sum;
}

inline Vec3d & Vec3d::operator+= (const Vec3d & vec2)
{
  elt[0] += vec2.elt[0];
  elt[1] += vec2.elt[1];
  elt[2] += vec2.elt[2];

  return *this;
}

inline Vec3d Vec3d::operator- (const Vec3d & vec2)
{
  Vec3d sum = *this;
  sum.elt[0] -= vec2.elt[0];
  sum.elt[1] -= vec2.elt[1];
  sum.elt[2] -= vec2.elt[2];

  return sum;
}

inline Vec3d & Vec3d::operator-= (const Vec3d & vec2)
{
  elt[0] -= vec2.elt[0];
  elt[1] -= vec2.elt[1];
  elt[2] -= vec2.elt[2];

  return *this;
}

inline double & Vec3d::operator[] (int index)
{
  return elt[index];
}

inline const double & Vec3d::operator[] (int index) const
{
  return elt[index];
}

inline double dot(const Vec3d & vec1, const Vec3d & vec2)
{
  return (vec1.elt[0] * vec2.elt[0] + vec1.elt[1] * vec2.elt[1] + vec1.elt[2] * vec2.elt[2]);
}

inline Vec3d cross(const Vec3d & vec1, const Vec3d & vec2)
{
  Vec3d result(vec1.elt[1] * vec2.elt[2] - vec2.elt[1] * vec1.elt[2],
	          -vec1.elt[0] * vec2.elt[2] + vec2.elt[0] * vec1.elt[2],
			   vec1.elt[0] * vec2.elt[1] - vec2.elt[0] * vec1.elt[1]);

  return result;
}

inline Vec3d norm(const Vec3d & vec1)
{
  double norm2 = dot(vec1,vec1);
  Vec3d result = vec1;
  result *= 1.0 / sqrt(norm2);
  
  return result;
}

inline Vec3d & Vec3d::operator*= (double scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  elt[2] *= scalar;
  return *this;
}

inline Vec3d Vec3d::operator* (double scalar) const
{
  return (Vec3d(elt[0]*scalar,elt[1]*scalar,elt[2]*scalar));
}

inline Vec3d Vec3d::operator/ (double scalar) const
{
  return (Vec3d(elt[0]/scalar,elt[1]/scalar,elt[2]/scalar));
}

inline Vec3d & Vec3d::operator/= (double scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  elt[2] /= scalar;
  return *this;
}


inline double len(const Vec3d & vec1)
{
  return(sqrt(dot(vec1,vec1)));
}

inline std::ostream &operator << (std::ostream &s, const Vec3d &v)
{
  double a = v[0];
  double b = v[1];
  double c = v[2];
  
  return(s << '[' << a << ' ' << b << ' ' << c << ']');
}

#endif
