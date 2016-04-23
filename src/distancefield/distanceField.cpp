#include "distanceField.h"
#include "trilinearInterpolation.h"
#include <string.h>
#include <float.h>

//#define GENERATE_DEBUG_DATA

DistanceField::DistanceField()
{ 
    distanceData = NULL;
    gradientData = NULL;
    pseudoData = NULL;
    setAutomaticBoundingBox(); 
}

DistanceField::~DistanceField()
{
    free(distanceData); 
    free(gradientData); 
    free(pseudoData); 
}

void DistanceField::setAutomaticBoundingBox(bool allBoxSidesEqual, double expansionRatio)
{
    useAutomaticBox = true;
    this->expansionRatio = expansionRatio;
    this->allBoxSidesEqual = allBoxSidesEqual;
}


// the routines for signed and unsigned distance field computation
#undef COMPUTE_CLOSEST_POINT
#define COMPUTE_SIGNED_FIELD
#include "computeField.cpp"
#undef COMPUTE_SIGNED_FIELD
#include "computeField.cpp"

int DistanceField::load(const std::string& filename)
{
    ifstream fin(filename.c_str(),ios::binary);
    if (!fin)
        return 1;

    fin.read((char*)&resolution,sizeof(int));

    // the type of data (single-precision or double-precision) is encoded
    //   as the sign of the x-resolution
    bool floatData = (resolution < 0);

    // note: all three resolutions (x,y,z) are assumed to be equal
    fin.read((char*)&resolution,sizeof(int));

    if (resolution < 0) // negative second resolution implies closest point data
        return 1;

    fin.read((char*)&resolution,sizeof(int));

    fin.read((char*)&(bmin_[0]),sizeof(double));
    fin.read((char*)&(bmin_[1]),sizeof(double));
    fin.read((char*)&(bmin_[2]),sizeof(double));

    fin.read((char*)&(bmax_[0]),sizeof(double));
    fin.read((char*)&(bmax_[1]),sizeof(double));
    fin.read((char*)&(bmax_[2]),sizeof(double));

    side = bmax_ - bmin_;
    gridX = side[0] / resolution;
    gridY = side[1] / resolution;
    gridZ = side[2] / resolution;
    setInvGridXYZ();

    distanceData = (float*) malloc
        (sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
    // place for one slice
    double * buffer = (double*) malloc
        (sizeof(double) * (resolution+1) * (resolution+1));

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

    /*
       if (!floatData) 
       fin.read((char*)distanceData,sizeof(double)*(resolution+1)*(resolution+1)*(resolution+1));
       else
       {
       fin.read((char*)distanceData,sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
    // position the data into place in a backwards pass
    for(int i=(resolution+1)*(resolution+1)*(resolution+1)-1; i>=0; i--)
    distanceData[i] = ((float*)distanceData)[i];       
     */

    fin.close();

    return 0;
}


int DistanceField::openStreamDistanceField(const std::string& filename,
        Vec3d * bmin, Vec3d * bmax, int * resolution,
        bool * floatData, ifstream & fin)
{
    fin.open(filename.c_str(),ios::binary);
    if (!fin)
        return 1;

    fin.read((char*)resolution,4);

    // the type of data (single-precision or double-precision) is encoded
    //   as the sign of the x-resolution
    *floatData = (*resolution < 0);

    // note: all three resolutions (x,y,z) are assumed to be equal
    fin.read((char*)resolution,4);

    if (*resolution < 0) // negative second resolution implies closest point data
        return 1;

    fin.read((char*)resolution,4);

    double bminx, bminy, bminz;
    fin.read((char*)&bminx,8);
    fin.read((char*)&bminy,8);
    fin.read((char*)&bminz,8);

    double bmaxx, bmaxy, bmaxz;
    fin.read((char*)&bmaxx,8);
    fin.read((char*)&bmaxy,8);
    fin.read((char*)&bmaxz,8);

    (*bmin)[0] = bminx;
    (*bmin)[1] = bminy;
    (*bmin)[2] = bminz;

    (*bmax)[0] = bmaxx;
    (*bmax)[1] = bmaxy;
    (*bmax)[2] = bmaxz;

    return 0;
}

void DistanceField::retrieveZSlice(ifstream & fin, bool floatData,
        int resolution, float * slice)
{
    // place for one slice
    double * buffer = (double*) malloc
        (sizeof(double) * (resolution+1) * (resolution+1));

    float * distanceDataPos = slice;

    if (floatData)
        fin.read((char*)buffer,sizeof(float)*(resolution+1)*(resolution+1));
    else
        fin.read((char*)buffer,sizeof(double)*(resolution+1)*(resolution+1));

    char * bufferPos = (char*)buffer;
    int bufferPosInc = floatData ? sizeof(float) : sizeof(double);

    // copy data to the slice
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

    /*
       if (!floatData) 
       fin.read((char*)slice,sizeof(double)*(resolution+1)*(resolution+1));
       else
       {
       fin.read((char*)slice,sizeof(float)*(resolution+1)*(resolution+1));
    // position the data into place in a backwards pass
    for(int i=(resolution+1)*(resolution+1)-1; i>=0; i--)
    slice[i] = ((float*)slice)[i];       
    }
     */
}

void DistanceField::set(int resolution, Vec3d bmin_, Vec3d bmax_,
        float * distanceData)
{
    this->resolution = resolution;
    this->bmin_ = bmin_; 
    this->bmax_ = bmax_; 

    side = bmax_ - bmin_;
    gridX = side[0] / resolution;
    gridY = side[1] / resolution;
    gridZ = side[2] / resolution;
    setInvGridXYZ();

    free(this->distanceData);
    this->distanceData = (float*) malloc
        (sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
    memcpy(this->distanceData, distanceData,
            sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
}

int DistanceField::isDoublePrecision(const std::string & filename,
        bool & doublePrecision)
{
    ifstream fin(filename.c_str(),ios::binary);
    if (!fin)
        return 1;

    int res;
    fin.read((char*)&res,4);

    // the type of data (single-precision or double-precision) is encoded
    //   as the sign of the x-resolution
    doublePrecision = (res >= 0);

    fin.close();

    return 0;
}

int DistanceField::saveToText(const std::string& filename)
{
    FILE * fout = fopen((char*)filename.c_str(), "wa");
    if (!fout)
        return 1;

    fprintf(fout, "%d\n", resolution);
    fprintf(fout, "%d\n", resolution);
    fprintf(fout, "%d\n", resolution);

    for(int i=0; i<=resolution; i++)
        for(int j=0; j<=resolution; j++)
            for(int k=0; k<=resolution; k++)
            {
                fprintf(fout,"%G\n", distance(i,j,k));
            }

    fclose(fout);

    return 0;
}

int DistanceField::save(const std::string& filename, bool doublePrecision)
{
    if (doublePrecision)
        printf("Error: double precision output not supported. Using float instead.\n");
    doublePrecision = false;

    ofstream fout(filename.c_str(),ios::binary);

    int data = resolution;
    if (!doublePrecision)
        data = -data;

    fout.write((char*)&data,sizeof(int));
    fout.write((char*)&resolution,sizeof(int));
    fout.write((char*)&resolution,sizeof(int));

    fout.write((char*)&(bmin_[0]),sizeof(double));
    fout.write((char*)&(bmin_[1]),sizeof(double));
    fout.write((char*)&(bmin_[2]),sizeof(double));

    fout.write((char*)&(bmax_[0]),sizeof(double));
    fout.write((char*)&(bmax_[1]),sizeof(double));
    fout.write((char*)&(bmax_[2]),sizeof(double));

    // float precision
    fout.write((char*)distanceData,
            sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));

    /*
       if (doublePrecision)
       fout.write((char*)distanceData,sizeof(double)*(resolution+1)*(resolution+1)*(resolution+1));
       else
       {
       float * buffer;
       convertToFloat(&buffer);
       fout.write((char*)buffer,sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
       free(buffer);
       }
     */
    fout.close();

#ifdef GENERATE_DEBUG_DATA
    ofstream fout1("debugPseudo",ios::binary);
    fout1.write((char*)pseudoData,
            6*sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
    fout1.close();
#endif

    return 0;
}

typedef struct
{
    int i,j,k;
    int di,dj,dk;
    double fieldDist, gridDist, relError; 
} errorData;

struct more_errorData : public std::binary_function< errorData, errorData, bool > {
    bool operator()(const errorData& x, const errorData& y) {
        return((x.relError) > (y.relError)); 
    }
};

bool DistanceField::sanityCheck()
{
    //double diagonal = sqrt(gridX*gridX + gridY*gridY + gridZ*gridZ);

    bool result = true;
    int count = 0;

    const int numErrorsPrinted = 3;

    vector<errorData> relErrors;
    errorData emptyEntry = {0,0,0,0,0,0,0.0,0.0,-1.0};
    for(int i=0; i<numErrorsPrinted; i++)
        relErrors.push_back(emptyEntry);

    for (int k=0; k <= resolution; k++)
    {
        for (int j=0; j <= resolution; j++)
            for (int i=0; i <= resolution; i++)
            {
                float d = distance(i,j,k);
                float d1;
                float sanity;
                float relError;
                float fieldDist;
                float gridDist;

#define PROCESS(di,dj,dk)\
                if ((i+(di) <= resolution) && (i+(di) >= 0) &&\
                        (j+(dj) <= resolution) && (j+(dj) >= 0) &&\
                        (k+(dk) <= resolution) && (k+(dk) >= 0))\
                {\
                    d1 = distance(i+(di),j+(dj),k+(dk));\
                    gridDist = len(Vec3d((di)*gridX,(dj)*gridY,(dk)*gridZ));\
                    fieldDist = fabs(d-d1);\
                    sanity = fieldDist - gridDist;\
                    if (sanity > 1E-6)\
                    {\
                        relError = sanity/gridDist;\
                        if (relError > relErrors[numErrorsPrinted-1].relError)\
                        {\
                            errorData errorEntry = {i,j,k,di,dj,dk,fieldDist,gridDist,relError};\
                            relErrors[numErrorsPrinted-1] = errorEntry;\
                            sort(relErrors.begin(),relErrors.end(),more_errorData());\
                        }\
                        result = false;\
                        count++;\
                    }\
                }

                PROCESS(1,0,0);
                PROCESS(1,1,0);
                PROCESS(0,1,0);
                PROCESS(-1,1,0);
                PROCESS(-1,0,0);
                PROCESS(-1,-1,0);
                PROCESS(0,-1,0);
                PROCESS(1,-1,0);

                PROCESS(0,0,1);
                PROCESS(1,0,1);
                PROCESS(1,1,1);
                PROCESS(0,1,1);
                PROCESS(-1,1,1);
                PROCESS(-1,0,1);
                PROCESS(-1,-1,1);
                PROCESS(0,-1,1);
                PROCESS(1,-1,1);

                PROCESS(0,0,-1);
                PROCESS(1,0,-1);
                PROCESS(1,1,-1);
                PROCESS(0,1,-1);
                PROCESS(-1,1,-1);
                PROCESS(-1,0,-1);
                PROCESS(-1,-1,-1);
                PROCESS(0,-1,-1);
                PROCESS(1,-1,-1);

            }
        cout << "." << flush;
    }



    cout << endl;

    if (count == 0)
        cout << "Success: sanity check passed." << endl;
    else
    {
        cout << "Encountered " << count << " possible errors. Largest top " << numErrorsPrinted << " errors (or all errors if fewer):" << endl;
        for(int i=0; i< (count < numErrorsPrinted ? count : numErrorsPrinted); i++)
        {
            errorData * errorEntry = &relErrors[i]; 
            float d1 = distance(errorEntry->i,errorEntry->j,errorEntry->k);
            float d2 = distance(errorEntry->i + errorEntry->di ,errorEntry->j + errorEntry->dj, errorEntry->k + errorEntry->dk);
            cout << "Distance field change too large. [" << errorEntry->i << "," << errorEntry->j << "," << errorEntry->k << "] to [" << errorEntry->i + (errorEntry->di) << "," << errorEntry->j + (errorEntry->dj) << "," << errorEntry->k + (errorEntry->dk) << "]" << " Dist 1: " << d1 << " Dist 2: " << d2 << " Reported change in distance field: " << errorEntry->fieldDist << " Grid distance: " << errorEntry->gridDist << " Relative error: " << errorEntry->relError << endl;
        }
    }
    return result;
}


void DistanceField::setBoundingBox(Vec3d bmin_g, Vec3d bmax_g)
{
    useAutomaticBox = false;
    bmin_ = bmin_g;
    bmax_ = bmax_g;
}

REAL DistanceField::distance(const Vector3d &pos) const
{
    Vec3d v_pos( pos[0], pos[1], pos[2] );

    return (REAL)( distance( v_pos ) );
}

float DistanceField::distance(const double &x, const double &y, const double &z, int constrainToBox) const
{
    int i,j,k;

    // get the index coordinate of the lower-right-bottom corner of the voxel containing 'pos'
    i = (int)((x - bmin_[0]) * invGridX);
    j = (int)((y - bmin_[1]) * invGridY);
    k = (int)((z - bmin_[2]) * invGridZ);

    if (((i<0) || (i>=resolution) || (j<0) || (j>=resolution) || (k<0) || (k>=resolution)) && (!constrainToBox))
    {
#if 0
        printf("Warning: querying the distance field outside of the bounding box: (i,j,k)=(%d,%d,%d), resolution=%d\n",i,j,k,resolution);
#endif
        return FLT_MAX;
    }

    if (constrainToBox)
    {
        if (i >= resolution)
            i = resolution - 1;
        if (i < 0)
            i = 0;
        if (j >= resolution)
            j = resolution - 1;
        if (j < 0)
            j = 0;
        if (k >= resolution)
            k = resolution - 1;
        if (k < 0)
            k = 0;
    }

    double wx,wy,wz;
    wx = ((x-bmin_[0]) / gridX) - i;
    wy = ((y-bmin_[1]) / gridY) - j;
    wz = ((z-bmin_[2]) / gridZ) - k;

    return TRILINEAR_INTERPOLATION(wx,wy,wz,distance(i,j,k),distance(i+1,j,k),distance(i+1,j+1,k),distance(i,j+1,k),
            distance(i,j,k+1),distance(i+1,j,k+1),distance(i+1,j+1,k+1),distance(i,j+1,k+1));
}

float DistanceField::distance(Vec3d pos, int constrainToBox) const
{
    distance(pos[0],pos[1],pos[2],constrainToBox); 
}

float DistanceField::maxValue()
{
    float maxValue=-FLT_MAX;
    for (int i=0; i <= resolution; i++)
        for (int j=0; j <= resolution; j++)
            for (int k=0; k <= resolution; k++)
            {
                float dist = distance(i,j,k);
                if (dist > maxValue)
                    maxValue = dist;
            }
    return maxValue;
}

float DistanceField::minValue()
{
    float minValue=FLT_MAX;
    for (int i=0; i <= resolution; i++)
        for (int j=0; j <= resolution; j++)
            for (int k=0; k <= resolution; k++)
            {
                float dist = distance(i,j,k);
                if (dist < minValue)
                    minValue = dist;
            }
    return minValue;
}

float DistanceField::maxAbsValue()
{
    float maxValue=0;
    for (int i=0; i <= resolution; i++)
        for (int j=0; j <= resolution; j++)
            for (int k=0; k <= resolution; k++)
            {
                float dist = fabs(distance(i,j,k));
                if (dist > maxValue)
                    maxValue = dist;
            }
    return maxValue;
}

float DistanceField::maxNonInftyAbsValue()
{
    float maxValue=0;
    for (int i=0; i <= resolution; i++)
        for (int j=0; j <= resolution; j++)
            for (int k=0; k <= resolution; k++)
            {
                float dist = fabs(distance(i,j,k));
                if ((dist > maxValue) && (dist != FLT_MAX))
                    maxValue = dist;
            }
    return maxValue;
}

float DistanceField::maxAbsValue(float threshold)
{
    float maxValue=0;
    for (int i=0; i <= resolution; i++)
        for (int j=0; j <= resolution; j++)
            for (int k=0; k <= resolution; k++)
            {
                float dist = fabs(distance(i,j,k));
                if ((dist > maxValue) && (dist < threshold))
                    maxValue = dist;
            }
    return maxValue;
}

Vec3d DistanceField::gradient(Vec3d pos) const
{
    int i,j,k;

    // get the indices
    i = (int)((pos[0] - bmin_[0]) * invGridX);
    j = (int)((pos[1] - bmin_[1]) * invGridY);
    k = (int)((pos[2] - bmin_[2]) * invGridZ);

    if ((i<=0) || (i>=resolution) || (j<=0) || (j>=resolution) || (k<=0) || (k>=resolution))
    {
        return Vec3d(0,0,0);
    }

    double wx,wy,wz;
    wx = ((pos[0]-bmin_[0]) / gridX) - i;
    wy = ((pos[1]-bmin_[1]) / gridY) - j;
    wz = ((pos[2]-bmin_[2]) / gridZ) - k;

    // gradient with respect to trilinear interpolation

    float v000 = distance(i,j,k);
    float v100 = distance(i+1,j,k);
    float v110 = distance(i+1,j+1,k);
    float v010 = distance(i,j+1,k);
    float v001 = distance(i,j,k+1);
    float v101 = distance(i+1,j,k+1);
    float v111 = distance(i+1,j+1,k+1);
    float v011 = distance(i,j+1,k+1);

    return Vec3d(
            GRADIENT_COMPONENT_X(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011),
            GRADIENT_COMPONENT_Y(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011),
            GRADIENT_COMPONENT_Z(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011) );

}

Vector3d DistanceField::gradient( const Vector3d &pos ) const
{
    Vec3d v_pos( pos[0], pos[1], pos[2] );
    Vec3d v_gradient = gradient( v_pos );

    return Vector3d( v_gradient[ 0 ], v_gradient[ 1 ], v_gradient[ 2 ] );
}

// Checks to see if a voxel is either fully or partially inside an object
bool DistanceField::voxelInside( int i, int j, int k, bool full,
        REAL tolerance )
{
    bool                       inside = full ? true : false;

    if ( full )
    {
        for ( int voxel_idx_x = i; voxel_idx_x <= i + 1; voxel_idx_x++ )
            for ( int voxel_idx_y = j; voxel_idx_y <= j + 1; voxel_idx_y++ )
                for ( int voxel_idx_z = k; voxel_idx_z <= k + 1; voxel_idx_z++ )
                {
                    inside = inside
                        && distance( voxel_idx_x, voxel_idx_y, voxel_idx_z ) < tolerance;
                }
    }
    else
    {
        for ( int voxel_idx_x = i; voxel_idx_x <= i + 1; voxel_idx_x++ )
            for ( int voxel_idx_y = j; voxel_idx_y <= j + 1; voxel_idx_y++ )
                for ( int voxel_idx_z = k; voxel_idx_z <= k + 1; voxel_idx_z++ )
                {
                    inside = inside
                        || distance( voxel_idx_x, voxel_idx_y, voxel_idx_z ) < tolerance;
                }
    }

    return inside;
}

void DistanceField::getDistanceData(float ** floatBuffer)
{
    unsigned int size = (resolution+1) * (resolution+1) * (resolution+1);
    *floatBuffer = (float*) malloc (sizeof(float) * size);
    memcpy(*floatBuffer, distanceData, sizeof(float) * size);
    //for(unsigned int i=0; i< size; i++)
    //(*floatBuffer)[i] = (float) (distanceData[i]);
}


/*
   void DistanceField::convertToFloat(float ** floatBuffer)
   {
   unsigned int size = (resolution+1) * (resolution+1) * (resolution+1);
 *floatBuffer = (float*) malloc (sizeof(float) * size);
 for(unsigned int i=0; i< size; i++)
 (*floatBuffer)[i] = (float) (distanceData[i]);
 }
 */

/*
   void DistanceField::convertToDouble(double ** doubleBuffer)
   {
   unsigned int size = (resolution+1) * (resolution+1) * (resolution+1);
 *doubleBuffer = (double*) malloc (sizeof(double) * size);
 memcpy(*doubleBuffer, distanceData, sizeof(double) * size);
 }
 */

bool DistanceField::isSurfaceVoxel(int i, int j, int k)
{
    float v[8];
    v[0] = distance(i,j,k);
    v[1] = distance(i+1,j,k);
    v[2] = distance(i+1,j+1,k);
    v[3] = distance(i,j+1,k);
    v[4] = distance(i,j,k+1);
    v[5] = distance(i+1,j,k+1);
    v[6] = distance(i+1,j+1,k+1);
    v[7] = distance(i,j+1,k+1);

    bool allPositive = true;
    for(int l=0; l<8; l++)
        if (v[l] < 0)
        {
            allPositive = false;
            break;
        }

    bool allNegative = true;
    for(int l=0; l<8; l++)
        if (v[l] >= 0)
        {
            allNegative = false;
            break;
        }

    return (!allNegative && !allPositive);
}

bool DistanceField::isSurfaceVoxel(int customResolution, int i, int j, int k, float levelSetValue)
{
    //printf("i: %d j:%d k:%d\n",i,j,k);
    Vec3d customGrid = 1.0 * side / customResolution;

    Vec3d basePos = bmin_ + Vec3d( 1.0 * i * customGrid[0], 1.0 * j * customGrid[1], 1.0 * k * customGrid[2] );

    if (((i < 0) || (i >= customResolution)) && ((j < 0) || (j >= customResolution)) && ((k < 0) || (k >= customResolution)))
    {
        printf("Warning (inside isSurfaceVoxel): voxel insides specified outside of the bounding box: (i,j,k)=(%d,%d,%d), customResolution=%d\n",i,j,k,customResolution);
    }

    // pass parameter 1 to constrain to box
    float v[8];
    v[0] = distance(basePos, 1) - levelSetValue;
    v[1] = distance(basePos + Vec3d(customGrid[0], 0, 0), 1) - levelSetValue;
    v[2] = distance(basePos + Vec3d(customGrid[0], customGrid[1], 0), 1) - levelSetValue;
    v[3] = distance(basePos + Vec3d(0, customGrid[1], 0), 1) - levelSetValue;
    v[4] = distance(basePos + Vec3d(0, 0, customGrid[2]), 1) - levelSetValue;
    v[5] = distance(basePos + Vec3d(customGrid[0], 0, customGrid[2]), 1) - levelSetValue;
    v[6] = distance(basePos + Vec3d(customGrid[0], customGrid[1], customGrid[2]), 1) - levelSetValue;
    v[7] = distance(basePos + Vec3d(0, customGrid[1], customGrid[2]), 1) - levelSetValue;

    /*
       double v[8];
       v[0] = distance(i,j,k);
       v[1] = distance(i+1,j,k);
       v[2] = distance(i+1,j+1,k);
       v[3] = distance(i,j+1,k);
       v[4] = distance(i,j,k+1);
       v[5] = distance(i+1,j,k+1);
       v[6] = distance(i+1,j+1,k+1);
       v[7] = distance(i,j+1,k+1);
     */

    bool allPositive = true;
    for(int l=0; l<8; l++)
        if (v[l] < 0)
        {
            allPositive = false;
            break;
        }

    bool allNegative = true;
    for(int l=0; l<8; l++)
        if (v[l] >= 0)
        {
            allNegative = false;
            break;
        }

    return (!allNegative && !allPositive);
}

int DistanceField::numSurfaceVoxels(float levelSetValue)
{
    return numSurfaceVoxels(resolution, levelSetValue);
}

int DistanceField::numSurfaceVoxels(int customResolution, float levelSetValue)
{
    printf("num surface voxels... custom res: %d\n", customResolution);
    int count = 0;

    for(int i=0; i<customResolution; i++)
        for(int j=0; j<customResolution; j++)
            for(int k=0; k<customResolution; k++)
                if (isSurfaceVoxel(customResolution,i,j,k, levelSetValue))
                    count++;

    return count;
}

float DistanceField::minBoundaryValue()
{
    float minDistance = FLT_MAX;
        
        for(int i=0; i <= resolution; i += resolution)
            for(int j=0; j <= resolution; j++)
                for(int k=0; k <= resolution; k++)
                {
                    float value = distance(i,j,k);
                        if (value < minDistance)
                            minDistance = value;
                }
    
        for(int i=0; i <= resolution; i++)
            for(int j=0; j <= resolution; j += resolution)
                for(int k=0; k <= resolution; k++)
                {
                    float value = distance(i,j,k);
                        if (value < minDistance)
                            minDistance = value;
                }
    
        for(int i=0; i <= resolution; i++)
            for(int j=0; j <= resolution; j++)
                for(int k=0; k <= resolution; k += resolution)
                {
                    float value = distance(i,j,k);
                        if (value < minDistance)
                            minDistance = value;
                }
    
        return minDistance;
}

