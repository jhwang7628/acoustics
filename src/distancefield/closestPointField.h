#ifndef _CLOSEST_POINT_FIELD_H_
#define _CLOSEST_POINT_FIELD_H_

#include "distanceField.h"

#if 0
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/VEC3.h>
#endif

#include <config.h>

// same as distance field, except it also stores the nearest point, not just the distance
class ClosestPointField : public DistanceField
{
public:
	// A feature class used to track the closest feature (point on
	// a triangle represented by barycentric coordinates) to a point.
	class Feature {
		public:
			Feature()
			{
			}

			Feature( int index1, int index2, int index3,
							 REAL alpha, REAL beta, REAL gamma )
				: index1( index1 ), index2( index2 ), index3( index3 ),
					alpha( alpha ), beta( beta ), gamma( gamma )
			{
			}

			int index1, index2, index3;
			REAL alpha, beta, gamma;
	};

	ClosestPointField();
	~ClosestPointField();

	int computeUnsignedField(const std::string& filename,
													 int resolution, int maxTriCount=15,
													 int maxDepth=10);
	int computeSignedField(const std::string& filename,
												 int resolution, int maxTriCount=15,
												 int maxDepth=10);

	int computeUnsignedField(ObjFile * objFile, int resolution,
													 int maxTriCount=15, int maxDepth=10);
	int computeSignedField(ObjFile * objFile, int resolution,
												 int maxTriCount=15, int maxDepth=10);

#if 0
	int computeUnsignedField(TriMesh *mesh, int resolution,
													 int maxTriCount=15, int maxDepth=10);
#endif

	int load(const std::string& filename);
	int save(const std::string& filename, bool doublePrecision=true);
	void set(int resolution, Vec3d bmin_, Vec3d bmax_,
					 float * distanceData, float * closestPointData,
					 Feature *closestFeatureData);

	void closestPoint(float pos[3], float result[3]);

	inline void closestPoint(int i, int j, int k, float result[3])
	{
		int pos = 3*((k * (resolution+1) + j ) * (resolution+1) + i);
		result[0] = closestPointData[pos];
		result[1] = closestPointData[pos+1];
		result[2] = closestPointData[pos+2];
	}

	bool sanityCheck();

	void getClosestPointData(float ** floatBuffer);

	inline void setClosestPoint(int i, int j, int k, float value[3])
	{
		int pos = 3*((k * (resolution+1) + j ) * (resolution+1) + i);
		closestPointData[pos] = value[0];
		closestPointData[pos+1] = value[1];
		closestPointData[pos+2] = value[2];
	}

	void closestFeature(int i, int j, int k, Feature &result)
	{
		int pos = (k * (resolution+1) + j ) * (resolution+1) + i;
		result = closestFeatureData[pos];
	}

#if 0
	// Interpolate 3-vector and scalar quantities using
	// closest point data.  For scalar quantities, the
	// input vector should have size equal to the number
	// of vertices in the input mesh.  For vector quantities,
	// the input vector should be three times this large.
	void interpolateScalar( const VEC3F &pos, VECTOR &data, REAL &result );
	void interpolateVector( const VEC3F &pos, VECTOR &data, VEC3F &result );
#endif

protected:
	float * closestPointData;

	// Store closest feature data, in addition to closest points
	Feature * closestFeatureData;
};

#endif

