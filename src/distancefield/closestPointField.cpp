#include "closestPointField.h"
#include "trilinearInterpolation.h"
#include <string.h>

ClosestPointField::ClosestPointField()
	: DistanceField()
{
	closestPointData = NULL;
	closestFeatureData = NULL;
}

ClosestPointField::~ClosestPointField()
{
	free(closestPointData);
	delete[] closestFeatureData;
}

// the routines for signed/unsigned closest point field computation
#define COMPUTE_CLOSEST_POINT
	#define COMPUTE_SIGNED_FIELD
		#include "computeField.cpp"
	#undef COMPUTE_SIGNED_FIELD
		#include "computeField.cpp"
#undef COMPUTE_CLOSEST_POINT

void ClosestPointField::set(int resolution, Vec3d bmin_, Vec3d bmax_,
														float * distanceData, float * closestPointData,
														Feature *closestFeatureData)
{
	DistanceField::set(resolution, bmin_, bmax_, distanceData);

	free(this->closestPointData);
	this->closestPointData = (float*) malloc
		(sizeof(float)*3*(resolution+1)*(resolution+1)*(resolution+1));
	memcpy(this->closestPointData, closestPointData,
		sizeof(float)*3*(resolution+1)*(resolution+1)*(resolution+1));

	if(this->closestFeatureData != NULL){
		delete this->closestFeatureData;
		this->closestFeatureData = NULL;
	}
	int size = (resolution+1) * (resolution+1) * (resolution+1);
	this->closestFeatureData = new Feature[size];
	for (int i = 0; i < size; i++)
	{
		this->closestFeatureData[i] = closestFeatureData[i];
	}
}

int ClosestPointField::load(const std::string& filename)
{
	ifstream fin(filename.c_str(),ios::binary);
	if (!fin)
		return 1;

	fin.read((char*)&resolution,4);

	// the type of data (single-precision or double-precision) is encoded	//	 as the sign of the x-resolution
	bool floatData = (resolution < 0);
	
	fin.read((char*)&resolution,4); // note: all three resolutions (x,y,z) are assumed to be equal
	if (resolution >=0)
	{
		return 1; // by convention the second resolution must be negative, so that we can differentiate closest point files from distance files
	}
	fin.read((char*)&resolution,4);

	fin.read((char*)&(bmin_[0]),8);
	fin.read((char*)&(bmin_[1]),8);
	fin.read((char*)&(bmin_[2]),8);
	
	fin.read((char*)&(bmax_[0]),8);
	fin.read((char*)&(bmax_[1]),8);
	fin.read((char*)&(bmax_[2]),8);

	side = bmax_ - bmin_;
	gridX = side[0] / resolution;
	gridY = side[1] / resolution;
	gridZ = side[2] / resolution;
	setInvGridXYZ();

	double * buffer = (double*) malloc (sizeof(double) * 3 * (resolution+1) * (resolution+1)); // place for one slice

	// load distances
	distanceData = (float*) malloc (sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));

	float * distanceDataPos = distanceData;

	for(int k=0; k <= resolution; k++)
	{
		if (floatData)
			fin.read((char*)buffer,sizeof(float)*(resolution+1)*(resolution+1));
		else
			fin.read((char*)buffer,sizeof(double)*(resolution+1)*(resolution+1));

		char * bufferPos = (char*)buffer;
		int bufferPosInc = floatData ? sizeof(float) : sizeof(double);

		// copy data to internal distance field buffer
		for(int j = 0; j <= resolution; j++)
			for(int i = 0; i <= resolution; i++)
			{
				if (floatData)
					*distanceDataPos = *(float*)bufferPos;
				else
					*distanceDataPos = *(double*)bufferPos;
				distanceDataPos++;
				bufferPos += bufferPosInc;
			}
	}

	// load closest points
	closestPointData = (float*) malloc (sizeof(float)*3*(resolution+1)*(resolution+1)*(resolution+1));

	float * closestPosDataPos = closestPointData;

	for(int k=0; k <= resolution; k++)
	{
		if (floatData)
			fin.read((char*)buffer,sizeof(float)*3*(resolution+1)*(resolution+1));
		else
			fin.read((char*)buffer,sizeof(double)*3*(resolution+1)*(resolution+1));

		char * bufferPos = (char*)buffer;
		int bufferPosInc = floatData ? sizeof(float) : sizeof(double);

		// copy data to internal closest point field buffer
		for(int j = 0; j <= resolution; j++)
			for(int i = 0; i <= resolution; i++)
				for(int coord=0; coord<3; coord++)
				{
					if (floatData)
						*closestPosDataPos = *(float*)bufferPos;
					else
						*closestPosDataPos = *(double*)bufferPos;
					closestPosDataPos++;
					bufferPos += bufferPosInc;
				}
	}

	// Load closest features
	int size = (resolution+1) * (resolution+1) * (resolution+1);
	int sliceSize = (resolution+1) * (resolution+1);
	closestFeatureData = new Feature[size];

	char *featureBuffer = (char*)closestFeatureData;

	for (int k = 0; k <= resolution; k++)
	{
		fin.read(featureBuffer, sizeof(Feature) * sliceSize);

		featureBuffer += sliceSize * sizeof(Feature);
	}

	// FIXME - remove this
	Feature &f1 = closestFeatureData[0];
	Feature &f2 = closestFeatureData[size-1];

	cout << "Loading distance field..." << endl;
	cout << "First feature = " << f1.index1;
	cout << ", " << f1.index2;
	cout << ", " << f1.index3;
	cout << ", " << f1.alpha;
	cout << ", " << f1.beta;
	cout << ", " << f1.gamma << endl;;
	cout << "Last feature = " << f2.index1;
	cout << ", " << f2.index2;
	cout << ", " << f2.index3;
	cout << ", " << f2.alpha;
	cout << ", " << f2.beta;
	cout << ", " << f2.gamma << endl;;

	fin.close();

	return 0;
}

int ClosestPointField::save(const std::string& filename, bool doublePrecision)
{
	if (doublePrecision)
		printf("Error: double precision output not supported. Using float instead.\n");
	doublePrecision = false;

	ofstream fout(filename.c_str(),ios::binary);

	int data = resolution;
	if (!doublePrecision)
		data = -data;
	fout.write((char*)&data,4);
	data = -resolution;
	fout.write((char*)&data,4);
	fout.write((char*)&resolution,4);

	fout.write((char*)&(bmin_[0]),8);
	fout.write((char*)&(bmin_[1]),8);
	fout.write((char*)&(bmin_[2]),8);

	fout.write((char*)&(bmax_[0]),8);
	fout.write((char*)&(bmax_[1]),8);
	fout.write((char*)&(bmax_[2]),8);

	int size = (resolution+1)*(resolution+1)*(resolution+1);

/*
	float * buffer = NULL;
	if (!doublePrecision)
		convertToFloat(&buffer);
*/

	// write out distances
	// float precision
	fout.write((char*)distanceData, sizeof(float) * size);
	
/*
	if (doublePrecision)
		fout.write((char*)distanceData,sizeof(double)*(resolution+1)*(resolution+1)*(resolution+1));
	else
	{
		fout.write((char*)buffer,sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
	}
*/

	// write out closest points
	// float precision
	fout.write((char*)closestPointData, sizeof(float) * 3 * size);

/*
	if (doublePrecision)
		fout.write((char*)closestPointData,sizeof(double)*3*(resolution+1)*(resolution+1)*(resolution+1));
	else
	{
		fout.write((char*)(buffer+(resolution+1)*(resolution+1)*(resolution+1)),
			sizeof(float)*3*(resolution+1)*(resolution+1)*(resolution+1));
		free(buffer);
	}
*/

	// Write out the closest feature data
	fout.write((char*)closestFeatureData,
		sizeof(ClosestPointField::Feature) * size);

	fout.close();

	// FIXME - remove this
	Feature f1 = closestFeatureData[0];
	Feature f2 = closestFeatureData[size-1];

	cout << "Saving distance field..." << endl;
	cout << "Testing my patience..." << endl;
	cout << "First feature = " << f1.index1;
	cout << ", " << f1.index2;
	cout << ", " << f1.index3;
	cout << ", " << f1.alpha;
	cout << ", " << f1.beta;
	cout << ", " << f1.gamma << endl;;
	cout << "Last feature = " << f2.index1;
	cout << ", " << f2.index2;
	cout << ", " << f2.index3;
	cout << ", " << f2.alpha;
	cout << ", " << f2.beta;
	cout << ", " << f2.gamma << endl;;
	cout << "Testing my patience AGAIN..." << endl;

	return 0;
}

void ClosestPointField::getClosestPointData(float ** floatBuffer)
{
	unsigned int size = 3 * (resolution+1) * (resolution+1) * (resolution+1);
	*floatBuffer = (float*) malloc (sizeof(float) * size);
	memcpy(*floatBuffer, closestPointData, sizeof(float) * size);
}

void ClosestPointField::closestPoint(const Vector3d &pos, Vector3d &result)
{
    float result_float[3]; 
    float position_float[3] = {(float)pos[0], (float)pos[1], (float)pos[2]}; 

    closestPoint(position_float, result_float); 

    result[0] = (double) result_float[0]; 
    result[1] = (double) result_float[1]; 
    result[2] = (double) result_float[2]; 

}

void ClosestPointField::closestPoint(float pos[3], float result[3])
{
	int i,j,k;

	// get the indices
	i = (int)((pos[0] - bmin_[0]) * invGridX);
	j = (int)((pos[1] - bmin_[1]) * invGridY);
	k = (int)((pos[2] - bmin_[2]) * invGridZ);

	if ((i<0) || (i>=resolution) || (j<0) || (j>=resolution) || (k<0) || (k>=resolution))
	{
		printf("Warning: querying the closest point field outside of the bounding box.\n");
		result[0] = result[1] = result[2] = 0;
	}

	double wx,wy,wz;
	wx = ((pos[0]-bmin_[0]) / gridX) - i;
	wy = ((pos[1]-bmin_[1]) / gridY) - j;
	wz = ((pos[2]-bmin_[2]) / gridZ) - k;

	float c0[3], c1[3], c2[3], c3[3], 
				c4[3], c5[3], c6[3], c7[3];
	closestPoint(i,j,k,c0);
	closestPoint(i+1,j,k,c1);
	closestPoint(i+1,j+1,k,c2);
	closestPoint(i,j+1,k,c3);
	closestPoint(i,j,k+1,c4);
	closestPoint(i+1,j,k+1,c5);
	closestPoint(i+1,j+1,k+1,c6);
	closestPoint(i,j+1,k+1,c7);

	result[0] = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 c0[0], c1[0], c2[0], c3[0], c4[0], c5[0], c6[0], c7[0]);
	result[1] = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 c0[1], c1[1], c2[1], c3[1], c4[1], c5[1], c6[1], c7[1]);
	result[2] = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 c0[2], c1[2], c2[2], c3[2], c4[2], c5[2], c6[2], c7[2]);
}

void ClosestPointField::interpolateScalar( const Vector3d &pos,
                                           VECTOR &data, REAL &result )
{
	int i,j,k;

	// get the indices
	i = (int)((pos[0] - bmin_[0]) * invGridX);
	j = (int)((pos[1] - bmin_[1]) * invGridY);
	k = (int)((pos[2] - bmin_[2]) * invGridZ);

	if ((i<0) || (i>=resolution) || (j<0) || (j>=resolution) || (k<0) || (k>=resolution))
	{
		printf("Warning: querying the closest point field outside of the bounding box.\n");
		result = 0;
	}

	double wx,wy,wz;
	wx = ((pos[0]-bmin_[0]) / gridX) - i;
	wy = ((pos[1]-bmin_[1]) / gridY) - j;
	wz = ((pos[2]-bmin_[2]) / gridZ) - k;

	ClosestPointField::Feature f[8];
	closestFeature(i,j,k,f[0]);
	closestFeature(i+1,j,k,f[1]);
	closestFeature(i+1,j+1,k,f[2]);
	closestFeature(i,j+1,k,f[3]);
	closestFeature(i,j,k+1,f[4]);
	closestFeature(i+1,j,k+1,f[5]);
	closestFeature(i+1,j+1,k+1,f[6]);
	closestFeature(i,j+1,k+1,f[7]);

	// Find a scalar value pertaining to each corner
	REAL s[8];
	REAL s_i1, s_i2, s_i3;
	for (int i = 0; i < 8; i++)
	{
		// Get the values at triangle corners
		s_i1 = data( f[i].index1 );
		s_i2 = data( f[i].index2 );
		s_i3 = data( f[i].index3 );

		s[i] = f[i].alpha * s_i1 + f[i].beta * s_i2 + f[i].gamma * s_i3;
	}

	result = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
}

void ClosestPointField::interpolateVector( const Vector3d &pos,
                                           VECTOR &data, Vector3d &result )
{
	int i,j,k;

	// get the indices
	i = (int)((pos[0] - bmin_[0]) * invGridX);
	j = (int)((pos[1] - bmin_[1]) * invGridY);
	k = (int)((pos[2] - bmin_[2]) * invGridZ);

	if ((i<0) || (i>=resolution) || (j<0) || (j>=resolution) || (k<0) || (k>=resolution))
	{
		printf("Warning: querying the closest point field outside of the bounding box.\n");
		result[0] = result[1] = result[2] = 0;
	}

	double wx,wy,wz;
	wx = ((pos[0]-bmin_[0]) / gridX) - i;
	wy = ((pos[1]-bmin_[1]) / gridY) - j;
	wz = ((pos[2]-bmin_[2]) / gridZ) - k;

	ClosestPointField::Feature f[8];
	closestFeature(i,j,k,f[0]);
	closestFeature(i+1,j,k,f[1]);
	closestFeature(i+1,j+1,k,f[2]);
	closestFeature(i,j+1,k,f[3]);
	closestFeature(i,j,k+1,f[4]);
	closestFeature(i+1,j,k+1,f[5]);
	closestFeature(i+1,j+1,k+1,f[6]);
	closestFeature(i,j+1,k+1,f[7]);

	// Find a scalar value pertaining to each corner
	Vector3d v[8];
	Vector3d v_i1, v_i2, v_i3;
	for (int i = 0; i < 8; i++)
	{
		// Get the values at triangle corners
		v_i1[0] = data( 3 * f[i].index1 );
		v_i1[1] = data( 3 * f[i].index1 + 1 );
		v_i1[2] = data( 3 * f[i].index1 + 2 );

		v_i2[0] = data( 3 * f[i].index2 );
		v_i2[1] = data( 3 * f[i].index2 + 1 );
		v_i2[2] = data( 3 * f[i].index2 + 2 );

		v_i3[0] = data( 3 * f[i].index3 );
		v_i3[1] = data( 3 * f[i].index3 + 1 );
		v_i3[2] = data( 3 * f[i].index3 + 2 );

		v[i] = f[i].alpha * v_i1 + f[i].beta * v_i2 + f[i].gamma * v_i3;
	}

	result[0] = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 v[0][0], v[1][0], v[2][0], v[3][0], v[4][0], v[5][0], v[6][0], v[7][0]);
	result[1] = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 v[0][1], v[1][1], v[2][1], v[3][1], v[4][1], v[5][1], v[6][1], v[7][1]);
	result[2] = TRILINEAR_INTERPOLATION(wx,wy,wz,
		 v[0][2], v[1][2], v[2][2], v[3][2], v[4][2], v[5][2], v[6][2], v[7][2]);
}

bool ClosestPointField::sanityCheck()
{
	bool exitCode = DistanceField::sanityCheck();

	bool myExitCode = true;
	for (int k=0; k <= resolution; k++)
		for (int j=0; j <= resolution; j++)
			for (int i=0; i <= resolution; i++)
			{
				float distanceScalar = distance(i,j,k);
				float closestPoint_[3];
				closestPoint(i,j,k, closestPoint_);
				Vec3d gridPos = position(i,j,k);
				float distanceNorm = 
					sqrt((closestPoint_[0] - gridPos[0])*(closestPoint_[0] - gridPos[0]) +
							 (closestPoint_[1] - gridPos[1])*(closestPoint_[1] - gridPos[1]) +
							 (closestPoint_[2] - gridPos[2])*(closestPoint_[2] - gridPos[2]));

				if (distanceScalar - distanceNorm > 1E-6)
				{
					printf("(i,j,k)=(%d,%d,%d): distance=%G | norm of closest vector=%G\n", i,j,k, distanceScalar, distanceNorm);
					myExitCode = false;
				}
			}

	return exitCode && myExitCode;
}
