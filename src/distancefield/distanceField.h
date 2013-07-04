#ifndef _DISTANCEFIELD_H_
#define _DISTANCEFIELD_H_

#include "obj.h"

#if 0
class TriMesh;
#endif

#if 0
#include <linearalgebra/VEC3.h>
#endif

#if 0
#include <SETTINGS.h>
#endif
#include <config.h>

/*
   Jernej Barbic
   jernej.barbic@gmail.com
   CMU, 2003-2007
   (Research performed together with Prof Doug James)

   DO NOT DISTRIBUTE (meant for internal use at CMU and Cornell (Doug James's group))
   A version for public release is under preparation as of December 2007,
   but this is not yet that version!!

   This code generates a distance filed, either signed or unsigned, to the given triangle mesh.
   The field can be loaded/saved to a file.
   You can also lookup the field (once computed or loaded) at arbitrary locations inside the field (using trilinear interpolation).

   Input mesh need not have triangle faces (e.g. can have quads; they will be triangulated internally).
   For signed field generation, the input mesh must be a closed manifold mesh.
   Input mesh must be given in the .obj format:
http://www.royriggs.com/obj.html

By default, the bounding box will be a cube, obtained by fitting the smallest cube to the geometry, and then expanded (scaled) from its center by a factor of 1.5. You can provide your own bounding boxes. However, note: (1) the provided bounding box must cover all the geometry, and (2) bounding boxes that are not cubes were not (yet) tested

Currently, only equal distance field resolutions along the three axes are supported. E.g., you can do 256 x 256 x 256, but not, for example, 256 x 128 x 128.

The bounding box will be divided into resolution number of smaller cubes along each axis. The distance field will be computed at the vertices of these cubes. So, if resolution is 256, the bounding box will get divided into 256 x 256 x 256 voxels, and the distance field will be computed on the resulting 257 x 257 x 257 grid of voxel vertices. Note that when indexing voxels, the indices (i,j,k) will run from 0 to 255 inclusive, whereas when indexing voxel vertices (also called "grid vertices"), they will run from 0 to 256 inclusive.

Distance field data is stored at voxel vertices. In memory, distance field value at voxel vertex (i,j,k) is stored at location k * (resolution+1)*(resolution+1) + j * (resolution+1) + i .

Internally, the code builds an octree on top of the triangle mesh. There are two parameters that control this process (you can keep them at default values, which worked well in practice for us) :
the max depth of the octree (auxiliary internal structure) is "maxDepth"
the max number of triangles intersecting an octree cell is "maxTriCount"
Note: once max depth level is reached, the maxTriCount bound is not imposed any more
 */

class DistanceField
{
    public:

        DistanceField();
        virtual ~DistanceField();

        // computes unsigned distance field
        // filename is the obj filename.  Alternately, this can be used
        // with a TriMesh object.
        virtual int computeUnsignedField(const std::string& filename,
                int resolution,
                int maxTriCount=15,
                int maxDepth=10);
        virtual int computeUnsignedField(ObjFile * objFile, int resolution,
                int maxTriCount=15, int maxDepth=10);
#if 0
        virtual int computeUnsignedField(TriMesh *mesh, int resolution,
                int maxTriCount=15, int maxDepth=10);
#endif

        // computes signed distance field
        // filename is the obj filename
        virtual int computeSignedField(const std::string& filename, int resolution,
                int maxTriCount=15, int maxDepth=10);
        virtual int computeSignedField(ObjFile * objFile, int resolution,
                int maxTriCount=15, int maxDepth=10);

        // sets the bounding box within which the distance field will be computed
        // set it before calling computeSigned/UnsignedField
        void setBoundingBox(Vec3d bmin, Vec3d bmax);
        // set a tight-fitting bounding box around the model, and expand it by the
        // expansion ratio given
        // note: if you don't set any bounding boxes, you get an automaticBoundingBox
        // with 1.5 expansion ratio by default
        void setAutomaticBoundingBox(bool allBoxSidesEqual=true,
                double expansionRatio=1.5);

        // loads a previously computed distance field from a disk file
        virtual int load(const std::string& filename); // returns 0 on success

        // opens the distance field for stream processing
        // this is advanced routine, don't use in normal circumstances;
        // use "load" instead
        int openStreamDistanceField(const std::string& filename,
                Vec3d * bmin, Vec3d * bmax,
                int * resolution, bool * floatData,
                ifstream & fin);
        void retrieveZSlice(ifstream & fin, bool floatData, int resolution,
                float * slice);

        // saves the current distance field to a disk file (e.g. after computing
        // it once, for later fast reloading) 
        virtual int save(const std::string& filename, bool doublePrecision);

        // exports the distance field to a file, in text format
        int saveToText(const std::string& filename);

        // sets the distance field to the given external data
        void set(int resolution, Vec3d bmin_, Vec3d bmax_, float * distanceData);

        // Is data from the given file in single or double precision?
        // note: this class uses single precision everywhere, but some older code
        // used double precision, so routines were necessary to load binary distance
        // field files computed using old versions
        int isDoublePrecision(const std::string & filename, bool & doublePrecision);

        // return distance field value at grid vertex (i,j,k)
        // each of i,j,k must be an integer from {0, ..., resolution}
        inline float distance(int i, int j, int k) const
        {
            return distanceData[(k * (resolution+1) + j ) * (resolution+1) + i];
        }
        // returns a 32-bit unique voxel index
        inline int packedVoxelIndex(int i, int j, int k);

        // returns indices of voxel containing pos
        inline void voxelIndices(Vec3d pos, int * i, int * j, int * k);

        // tells whether the pos is inside box or not
        inline bool insideBox(Vec3d pos);

        // returns a 32-bit unique voxel index, for voxel constaining pos
        inline int packedVoxelIndex(Vec3d pos);

        // alters the distance at a particular grid vertex (i,j,k)
        inline void setDistance(int i, int j, int k, float value)
        {
            distanceData[(k * (resolution+1) + j ) * (resolution+1) + i] = value;
        }

#if 0
        REAL distance(const VEC3F &pos) const;
#endif
        // computes distance and gradient at arbitrary position
        float distance(Vec3d pos, int constrainToBox=0) const;
        // gradient is computed with respect to trilinear interpolation
        // note: gradient is discontinuous at the cell boundaries
        Vec3d gradient(Vec3d pos) const;
#if 0
        VEC3F gradient(const VEC3F &pos) const;
#endif

        // Checks to see if a voxel is either fully or partially inside an object
        bool voxelInside( int i, int j, int k, bool full,
                          REAL tolerance = 0.0 );

        // returns the world-coordinate position of the grid vertex with indices (i,j,k)
        // must have: 0 <= i,j,k <= resolution   (i,j,k of course not necessarily sorted)
        inline Vec3d position(int i, int j, int k);

#if 0
        inline VEC3F gridPosition( int i, int j, int k )
        {
            Vec3d pos = position( i, j, k );

            return VEC3F( pos[ 0 ], pos[ 1 ], pos[ 2 ] );
        }
#endif

        // get distance field resolution
        inline int getResolution() { return resolution;}

        // get the diagonal of the bounding box
        inline double diagonal() { return len(bmax_-bmin_);}

        // get the lower-left-front corner of bounding box
        inline Vec3d bmin() const { return bmin_; }
        // get the upper-right-back corner of bounding box
        inline Vec3d bmax() const { return bmax_; }
        // alternative interface to bmin, bmax functions above
        inline void bmin(float * bmin)
        {
            bmin[0] = bmin_[0]; bmin[1] = bmin_[1]; bmin[2] = bmin_[2];
        }
        inline void bmax(float * bmax)
        {
            bmax[0] = bmax_[0]; bmax[1] = bmax_[1]; bmax[2] = bmax_[2];
        }
        inline void bmin(double * bmin)
        {
            bmin[0] = bmin_[0]; bmin[1] = bmin_[1]; bmin[2] = bmin_[2];
        }
        inline void bmax(double * bmax)
        {
            bmax[0] = bmax_[0]; bmax[1] = bmax_[1]; bmax[2] = bmax_[2];
        }

        // checks if distance for any two adjacent voxels is less than voxel grid
        // spacing apart (which it must be by triangle inequality, for both signed
        // and unsigned fields)
        virtual bool sanityCheck();

        float maxValue();
        float minValue();
        float maxAbsValue();
        float maxAbsValue(float threshold); // only abs values up to threshold
        float maxNonInftyAbsValue();

        // minimum distance value on the surface of the bounding box
        float minBoundaryValue();

        // returns the entire distance data, by allocating the buffer and copying
        // the distance data into it 
        virtual void getDistanceData(float ** floatBuffer);

        // returns the spatial dimensions of the voxels
        // i.e., this is the spatial distance between consecutive grid vertices
        // along the three dimensions
        inline void getGridSpacing(double * gridX, double * gridY, double * gridZ);

        // Is voxel with indices (i,j,k) intersecting the zero-isocontour?
        // (only applies to sign distance fields)
        bool isSurfaceVoxel(int i, int j, int k); 
        // total number of such voxels
        int numSurfaceVoxels(float levelSetValue = 0.0); 
        // If distance field were resampled to customResolution, is voxel with
        // indices (i,j,k) intersecting the isocontour with value levelSetValue ?
        bool isSurfaceVoxel(int customResolution, int i, int j, int k,
                float levelSetValue);
        // total number of such voxels
        int numSurfaceVoxels(int customResolution, float levelSetValue); 

    protected:
        int maxTriCount;
        int resolution;
        int maxDepth;
        double expansionRatio;

        bool useAutomaticBox, allBoxSidesEqual;

        float * distanceData; // the raw distance data
        float * gradientData; // the raw gradient data
        float * pseudoData; // for debug

        Vec3d bmin_, bmax_;
        Vec3d side;
        double gridX, gridY, gridZ;
        double invGridX, invGridY, invGridZ;

        ObjFileOctree<TriangleWithCollisionInfo> * objFileOctree;
        ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals> * objFileOrientedOctree;

        inline void setInvGridXYZ();
};

inline void DistanceField::getGridSpacing(double * gridX, double * gridY,
        double * gridZ)
{
    *gridX = this->gridX;
    *gridY = this->gridY;
    *gridZ = this->gridZ;
}

inline Vec3d DistanceField::position(int i, int j, int k)
{
    return(
            Vec3d (bmin_[0] + 1.0 * i / resolution * side[0],
                bmin_[1] + 1.0 * j / resolution * side[1],
                bmin_[2] + 1.0 * k / resolution * side[2])
          );
}

inline void DistanceField::setInvGridXYZ()
{
    invGridX = 1.0 / gridX;
    invGridY = 1.0 / gridY;
    invGridZ = 1.0 / gridZ;
}

// returns indices of voxel containing pos
inline void DistanceField::voxelIndices(Vec3d pos, int * i, int * j, int * k)
{
    *i = (int)((pos[0] - bmin_[0]) * invGridX);
    *j = (int)((pos[1] - bmin_[1]) * invGridY);
    *k = (int)((pos[2] - bmin_[2]) * invGridZ);
}

inline bool DistanceField::insideBox(Vec3d pos)
{
    int i = (int)((pos[0] - bmin_[0]) * invGridX);
    int j = (int)((pos[1] - bmin_[1]) * invGridY);
    int k = (int)((pos[2] - bmin_[2]) * invGridZ);

    return ((i >= 0) && (i < resolution) && (j >= 0) && (j < resolution)
            && (k >= 0) && (k < resolution));
}

inline int DistanceField::packedVoxelIndex(Vec3d pos)
{
    // get the indices
    int i = (int)((pos[0] - bmin_[0]) * invGridX);
    int j = (int)((pos[1] - bmin_[1]) * invGridY);
    int k = (int)((pos[2] - bmin_[2]) * invGridZ);

    return (k * resolution + j ) * resolution + i;
}

inline int DistanceField::packedVoxelIndex(int i,int j,int k)
{
    return (k * resolution + j ) * resolution + i;
}

#endif
