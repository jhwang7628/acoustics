#ifndef _MINIVEC2D_H_
#define _MINIVEC2D_H_

#include <math.h>
#include <ostream>
#include "vec_types.h"
                                                                                                                                                             
class Vec2d {

public:

	inline Vec2d() {}
	inline Vec2d(double x_g, double y_g) {elt[0]=x_g; elt[1]=y_g;}
	inline Vec2d(ZeroOrOne k);

  inline Vec2d & operator=(const Vec2d & source);
  inline bool operator==(const Vec2d & vec2);
	inline Vec2d & operator=(ZeroOrOne k);

	inline Vec2d operator+ (const Vec2d & vec2);
	inline Vec2d & operator+= (const Vec2d & vec2);

	inline Vec2d operator- (const Vec2d & vec2);
	inline Vec2d & operator-= (const Vec2d & vec2);

	inline Vec2d operator* (double scalar) const;
	inline Vec2d & operator*= (double scalar);

	inline Vec2d operator/ (double scalar) const;
	inline Vec2d & operator/= (double scalar);

	friend inline Vec2d operator* (double scalar, const Vec2d & vec2);
	friend inline Vec2d operator/ (double scalar, const Vec2d & vec2);

	friend inline double dot(const Vec2d & vec1, const Vec2d & vec2);

  friend inline Vec2d norm(const Vec2d & vec1);

	friend class Mat3d;

	inline double & operator[] (int index);
	inline const double & operator[] (int index) const;

protected:
  double elt[2];

};


inline Vec2d::Vec2d(ZeroOrOne k)
{
  elt[0] = k;
  elt[1] = k;
}

inline Vec2d & Vec2d::operator=(const Vec2d & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];

  return *this;
}

inline bool Vec2d::operator==(const Vec2d & vec2)
{
  return ((elt[0] == vec2[0]) &&
          (elt[1] == vec2[1]));
}

inline Vec2d & Vec2d::operator=(ZeroOrOne k)
{
  elt[0] = k;
  elt[1] = k;

  return *this;
}


inline Vec2d operator* (double scalar, const Vec2d & vec2)
{
  Vec2d result = vec2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;

  return result;
}

inline Vec2d operator/ (double scalar, const Vec2d & vec2)
{
  Vec2d result = vec2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;

  return result;
}


inline Vec2d Vec2d::operator+ (const Vec2d & vec2)
{
  Vec2d sum = *this;
  sum.elt[0] += vec2.elt[0];
  sum.elt[1] += vec2.elt[1];

  return sum;
}

inline Vec2d & Vec2d::operator+= (const Vec2d & vec2)
{
  elt[0] += vec2.elt[0];
  elt[1] += vec2.elt[1];

  return *this;
}

inline Vec2d Vec2d::operator- (const Vec2d & vec2)
{
  Vec2d sum = *this;
  sum.elt[0] -= vec2.elt[0];
  sum.elt[1] -= vec2.elt[1];

  return sum;
}

inline Vec2d & Vec2d::operator-= (const Vec2d & vec2)
{
  elt[0] -= vec2.elt[0];
  elt[1] -= vec2.elt[1];

  return *this;
}

inline double & Vec2d::operator[] (int index)
{
  return elt[index];
}

inline const double & Vec2d::operator[] (int index) const
{
  return elt[index];
}

inline double dot(const Vec2d & vec1, const Vec2d & vec2)
{
  return (vec1.elt[0] * vec2.elt[0] + vec1.elt[1] * vec2.elt[1]);
}

inline Vec2d norm(const Vec2d & vec1)
{
  double norm2 = dot(vec1,vec1);
  Vec2d result = vec1;
  result *= 1.0 / sqrt(norm2);
  
  return result;
}

inline Vec2d & Vec2d::operator*= (double scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  return *this;
}

inline Vec2d Vec2d::operator* (double scalar) const
{
  return (Vec2d(elt[0]*scalar,elt[1]*scalar));
}

inline Vec2d Vec2d::operator/ (double scalar) const
{
  return (Vec2d(elt[0]/scalar,elt[1]/scalar));
}

inline Vec2d & Vec2d::operator/= (double scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  return *this;
}


inline double len(const Vec2d & vec1)
{
  return(sqrt(dot(vec1,vec1)));
}

inline std::ostream &operator << (std::ostream &s, const Vec2d &v)
{
  double a = v[0];
  double b = v[1];
  
  return(s << '[' << a << ' ' << b << ']');
}

#endif
