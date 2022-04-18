#ifndef _MINIMAT3D_H_
#define _MINIMAT3D_H_

#include "vec3d.h"

class Mat3d {

public:
  
	inline Mat3d() {}

	inline Mat3d(double x0_g, double x1_g, double x2_g,
		           double x3_g, double x4_g, double x5_g,
				       double x6_g, double x7_g, double x8_g);

  inline Mat3d(double mat[9]);

	inline Mat3d(ZeroOrOne k);

	inline Mat3d & operator=(const Mat3d & source);
	inline Mat3d & operator=(ZeroOrOne k);

	inline Mat3d operator+ (const Mat3d & );
	inline Mat3d & operator+= (const Mat3d & );

	inline Mat3d operator- (const Mat3d & );
	inline Mat3d & operator-= (const Mat3d & );

	inline Mat3d & operator*= (double scalar);
	inline Mat3d & operator/= (double scalar);

  friend inline Mat3d operator* (double scalar, const Mat3d & mat2);
  friend inline Mat3d operator/ (double scalar, const Mat3d & mat2);

	friend inline Mat3d TensorProduct(const Vec3d & vec1, const Vec3d & vec2);

  friend inline Mat3d inv(const Mat3d &);
  friend inline double det(Mat3d & mat);

  friend inline Mat3d trans(Mat3d & mat);

	inline Vec3d & operator[] (int index);

  inline void convertToArray(double * array); // in row-major order

	// matrix-vector multiply
	inline Vec3d operator* (Vec3d & vec);
	// matrix-matrix multiply
	inline Mat3d operator* (Mat3d & mat2);

  // computes eigenvalues and eigenvectors of the 3x3 matrix
  // assumes symmetric matrix; contents of matrix a are overwritten (destroyed)
  // returns true if iteration succeeded in making the sum of abs values of non-diagonal values below epsilon, and false otherwise
  // default epsilon is machine precision (which always converged with our matrices)
  // eigenvalues are sorted in decreasing abs order
  // returned eigenvectors are unit length
  bool eigen_sym(Mat3d & a, Vec3d & eig_val, Vec3d eig_vec[3],
                        int maxIterations=50, double epsilon=0.0);

protected:
  Vec3d elt[3];

};


inline Mat3d::Mat3d(double x0_g, double x1_g, double x2_g,
					double x3_g, double x4_g, double x5_g,
					double x6_g, double x7_g, double x8_g)
{
  elt[0] = Vec3d(x0_g,x1_g,x2_g);
  elt[1] = Vec3d(x3_g,x4_g,x5_g);
  elt[2] = Vec3d(x6_g,x7_g,x8_g);
}

inline Mat3d::Mat3d(double mat[9])
{
  elt[0] = Vec3d(mat[0],mat[1],mat[2]);
  elt[1] = Vec3d(mat[3],mat[4],mat[5]);
  elt[2] = Vec3d(mat[6],mat[7],mat[8]);
}

inline Mat3d::Mat3d(ZeroOrOne k)
{
  elt[0] = Vec3d(k,0,0);
  elt[1] = Vec3d(0,k,0);
  elt[2] = Vec3d(0,0,k);
}


inline Mat3d & Mat3d::operator=(const Mat3d & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];
  elt[2] = source.elt[2];

  return *this;
}

inline Mat3d & Mat3d::operator=(ZeroOrOne k)
{

  elt[0] = Vec3d(k,0,0);
  elt[1] = Vec3d(0,k,0);
  elt[2] = Vec3d(0,0,k);

  return *this;
}


inline Mat3d Mat3d::operator+ (const Mat3d & mat2)
{
  Mat3d sum = *this;
  sum.elt[0] += mat2.elt[0];
  sum.elt[1] += mat2.elt[1];
  sum.elt[2] += mat2.elt[2];

  return sum;
}

inline Mat3d & Mat3d::operator+= (const Mat3d & mat2)
{
  elt[0] += mat2.elt[0];
  elt[1] += mat2.elt[1];
  elt[2] += mat2.elt[2];
  return *this;
}


inline Mat3d Mat3d::operator- (const Mat3d & mat2)
{
  Mat3d sum = *this;
  sum.elt[0] -= mat2.elt[0];
  sum.elt[1] -= mat2.elt[1];
  sum.elt[2] -= mat2.elt[2];

  return sum;
}

inline Mat3d & Mat3d::operator-= (const Mat3d & mat2)
{
  elt[0] -= mat2.elt[0];
  elt[1] -= mat2.elt[1];
  elt[2] -= mat2.elt[2];

  return *this;
}



inline Vec3d & Mat3d::operator[] (int index)
{
  return elt[index];
}

inline Mat3d & Mat3d::operator*= (double scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  elt[2] *= scalar;

  return *this;

}

inline Mat3d operator* (double scalar, const Mat3d & mat2)
{
  Mat3d result = mat2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;
  result.elt[2] *= scalar;

  return result;
}

inline Mat3d operator/ (double scalar, const Mat3d & mat2)
{
  Mat3d result = mat2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;
  result.elt[2] /= scalar;

  return result;
}
 
inline Mat3d TensorProduct(Vec3d & vecA, Vec3d & vecB)
{
  Mat3d result(vecA[0]*vecB[0],vecA[0]*vecB[1],vecA[0]*vecB[2],
	       vecA[1]*vecB[0],vecA[1]*vecB[1],vecA[1]*vecB[2],
	       vecA[2]*vecB[0],vecA[2]*vecB[1],vecA[2]*vecB[2]);

  return result;
}


inline Mat3d & Mat3d::operator/= (double scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  elt[2] /= scalar;

  return *this;

}


inline Vec3d Mat3d::operator* (Vec3d & vec)
{
  return(Vec3d(
	  dot(elt[0],vec),
	  dot(elt[1],vec),
	  dot(elt[2],vec)));
}

inline Mat3d Mat3d::operator* (Mat3d & mat2)
{
  return(Mat3d(
    dot(elt[0],Vec3d(mat2.elt[0][0],mat2.elt[1][0],mat2.elt[2][0])),
	  dot(elt[0],Vec3d(mat2.elt[0][1],mat2.elt[1][1],mat2.elt[2][1])),
	  dot(elt[0],Vec3d(mat2.elt[0][2],mat2.elt[1][2],mat2.elt[2][2])),

    dot(elt[1],Vec3d(mat2.elt[0][0],mat2.elt[1][0],mat2.elt[2][0])),
	  dot(elt[1],Vec3d(mat2.elt[0][1],mat2.elt[1][1],mat2.elt[2][1])),
	  dot(elt[1],Vec3d(mat2.elt[0][2],mat2.elt[1][2],mat2.elt[2][2])),

    dot(elt[2],Vec3d(mat2.elt[0][0],mat2.elt[1][0],mat2.elt[2][0])),
	  dot(elt[2],Vec3d(mat2.elt[0][1],mat2.elt[1][1],mat2.elt[2][1])),
	  dot(elt[2],Vec3d(mat2.elt[0][2],mat2.elt[1][2],mat2.elt[2][2])) ));
}

inline Mat3d inv(Mat3d & mat)
{
  double invDeterminant = 1.0 / 
          (-mat[0][2] * mat[1][1] * mat[2][0] + 
            mat[0][1] * mat[1][2] * mat[2][0] + 
    	      mat[0][2] * mat[1][0] * mat[2][1] - 
	          mat[0][0] * mat[1][2] * mat[2][1] - 
	          mat[0][1] * mat[1][0] * mat[2][2] + 
	          mat[0][0] * mat[1][1] * mat[2][2] );


  return Mat3d(
    invDeterminant * (-mat[1][2] * mat[2][1] + mat[1][1] * mat[2][2]),
    invDeterminant * (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]),
    invDeterminant * (-mat[0][2] * mat[1][1] + mat[0][1] * mat[1][2]),
    invDeterminant * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]),
    invDeterminant * (-mat[0][2] * mat[2][0] + mat[0][0] * mat[2][2]),
    invDeterminant * (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]),
    invDeterminant * (-mat[1][1] * mat[2][0] + mat[1][0] * mat[2][1]),
    invDeterminant * (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]),
    invDeterminant * (-mat[0][1] * mat[1][0] + mat[0][0] * mat[1][1])
  );
}

inline double det(Mat3d & mat)
{
  return
          (-mat[0][2] * mat[1][1] * mat[2][0] + 
            mat[0][1] * mat[1][2] * mat[2][0] + 
    	      mat[0][2] * mat[1][0] * mat[2][1] - 
	          mat[0][0] * mat[1][2] * mat[2][1] - 
	          mat[0][1] * mat[1][0] * mat[2][2] + 
	          mat[0][0] * mat[1][1] * mat[2][2] );

}

inline Mat3d trans(Mat3d & mat)
{
  return
       Mat3d( mat[0][0], mat[1][0], mat[2][0],
              mat[0][1], mat[1][1], mat[2][1],
              mat[0][2], mat[1][2], mat[2][2] );
}

inline void Mat3d::convertToArray(double * array) // in row-major order
{
  array[0] = elt[0][0]; array[1] = elt[0][1]; array[2] = elt[0][2];
  array[3] = elt[1][0]; array[4] = elt[1][1]; array[5] = elt[1][2];
  array[6] = elt[2][0]; array[7] = elt[2][1]; array[8] = elt[2][2];
}

#endif
