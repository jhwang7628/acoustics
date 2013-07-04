#ifndef __OBJFILE_H__
#define __OBJFILE_H__

#ifdef WIN32
#include <windows.h>
#endif

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include "minivector.h"

#ifndef PI
#define PI 3.141592653589793238462643
#endif

// Jernej Barbic, CMU, 2003-2006
//
// Obj loader by Chris Twigg.
//
// This is a wrapper for an .obj file.  We can construct one by simply
// passing in a filename; assuming the obj file does not have any errors
// we will then be able to access any of the vertices, etc. of the obj file.

// The .obj file format is documented fully here:
//    http://www.royriggs.com/obj.html
// A quick summary:
//   1.  vertices, normals, and texture coordinates are all specified in 
//       a global 1-based namespace.  
//   2.  faces are divided into groups
//   3.  each face consists of a listing of vertices, like this:
//          f 1/1/1 2/2/2 3/3/3
//       The numbers are references to the vertices, normals, and texture
//       coordinates, all of which were specified (as mentioned above) in
//       a global 1-based namespace.  
//
// So, with that said, to get an object out of the ObjFile object below
// once it has been constructed, do the following:
//   1.  Get the list of groups out of the .obj file using the "groupNames"
//       function.
//   2.  Select the group you want, and retrieve it using the "group" function.
//   3.  Iterate through the faces in the group using the "face" member function
//       of the Group class.
//   4.  Iterate through the vertices using the vertex member function of 
//       the Face class.
//   5.  Retrieve the various indices for the position(), textureCoordinate(), and 
//       normal() attributes of the vertex.
//   6.  Look these up in the .obj file global namspace using the vertexPosition,
//       textureCoordinate, and normal member functions of the ObjFile class.


// thrown by the ObjFile constructor
class ObjFileException
{
  public:
    ObjFileException( const std::string& reason );
    ObjFileException( const std::string& reason, const std::string& filename, unsigned int line);

    std::string reason() const { return reason_; }

  private:
    std::string reason_;
};

class ObjFile
{
  public:

    class Vertex
    {
      public:
        Vertex( const unsigned int& pos )
          :	position_(pos),
          texture_( std::make_pair( false, 0 ) ),
          normal_( std::make_pair( false, 0 ) ) {}

        Vertex( const unsigned int& pos, const unsigned int& texture )
          :	position_(pos),
          texture_( std::make_pair( true, texture ) ),
          normal_( std::make_pair( false, 0 ) ) {}

        Vertex( const unsigned int& pos, const unsigned int& texture, const unsigned int& normal )
          :	position_(pos),
          texture_( std::make_pair( true, texture ) ),
          normal_( std::make_pair( true, normal ) ) {}

        Vertex( const unsigned int pos, std::pair<bool, unsigned int> texture, std::pair<bool, unsigned int> normal )
          :	position_(pos),
          texture_(texture),
          normal_(normal) {}

        unsigned int position() const			{ return position_; }
        unsigned int normal() const				{ assert( normal_.first ); return normal_.second; }
        unsigned int textureCoordinate() const	{ assert( texture_.first ); return texture_.second;}
        std::pair< bool, unsigned int > texturePair() const { return texture_; }
        std::pair< bool, unsigned int > normalPair() const { return normal_; }


        // Normals and texture coordinates are not considered "required" in the
        // obj file format standard.  Check these before retrieving them.
        bool hasNormal() const	{ return normal_.first; }
        bool hasTextureCoordinate() const { return texture_.first; }

        void setPosition(unsigned int position_g) { position_ = position_g; }
        void setNormal(unsigned int normal_g) { normal_ = std::pair<bool,unsigned int>(true,normal_g); }
        void setTextureCoordinate(unsigned int textureCoordinate_g) { texture_ = std::pair<bool,unsigned int>(true,textureCoordinate_g); }

      protected:
        unsigned int position_; 
        std::pair< bool, unsigned int > texture_;
        std::pair< bool, unsigned int > normal_;	
    };

    class Material
    {
      public:

        explicit Material(const std::string name, Vec3d Ka, Vec3d Kd, Vec3d Ks, double shininess=0): 
          Ka_(Ka), Kd_(Kd), Ks_(Ks), shininess_(shininess), name_(name) {}

        explicit Material(): 
          Ka_(Vec3d(1,1,1)), Kd_(Vec3d(1,1,1)), Ks_(Vec3d(1,1,1)), name_(std::string("default")) {}

        std::string name() const { return name_; }
        inline Vec3d Ka() { return Ka_; }
        inline Vec3d Kd() { return Kd_; }
        inline Vec3d Ks() { return Ks_; }
        inline double shininess() { return shininess_; }

        bool operator==(Material & material2);

      protected:
        Vec3d Ka_, Kd_, Ks_;
        double shininess_;
        std::string name_;
    };

    class Face
    {
      public:
        Face() { vertices_.reserve( 3 ); }
        Face( const Vertex& v1, const Vertex& v2, const Vertex& v3 )
        {
          vertices_.reserve(3);
          vertices_.push_back(v1);
          vertices_.push_back(v2);
          vertices_.push_back(v3);
        }

        void addVertex( const Vertex& v )
        {
          vertices_.push_back( v );
        }

        size_t vertexCount() const { return vertices_.size(); }
        Vertex vertex(unsigned int vert) { return vertices_[vert]; }
        Vertex * vertexHandle(unsigned int vert) { return &(vertices_[vert]); }

        void reverseVertices() { reverse(vertices_.begin(),vertices_.end()); }
        void displayVertices() { for(unsigned int i=0; i<vertices_.size(); i++) std::cout << vertices_[i].position() << " ";}

      protected:
        std::vector< Vertex > vertices_;
    };

    class Group
    {
      public:
        explicit Group( std::string name, unsigned int materialIndex=0 )
          : name_(name), materialIndex_(materialIndex) {}

        void addFace( const Face& face ) { faces_.push_back( face ); }

        size_t faceCount() const { return faces_.size(); }
        Face face(unsigned int face) const { return faces_[face]; }
        Face * faceHandle(unsigned int face) { return &(faces_[face]); }
        std::string name() const { return name_; }
        unsigned int materialIndex() const { return materialIndex_; }
        void setMaterialIndex(unsigned int index) { materialIndex_ = index;}

        void reverseFace(unsigned int i) { faces_[i].reverseVertices(); }

        void removeFace(unsigned int i);

      protected:
        std::string name_;
        unsigned int materialIndex_;
        std::vector< Face > faces_;
    };

    // Constructs the OBJ file and reads it in.  Throws an ObjFileException if
    // it fails for any reason (e.g., file not there, etc.).
    ObjFile( const std::string& filename, int verbose=1 );

    // makes an empty structure
    ObjFile() {}

    // Retrieve a group given a group name.
    Group group( std::string name ) const;

    // Retrieve a list of all the group names in the obj file.
    std::vector<std::string> groupNames() const;

    // Prints info on the obj model
    void PrintInfo() const; 

    bool isTriangularMesh();
    void triangulate();

    bool isQuadrilateralMesh();

    // expansionRatio of 1 just encloses the model
    // the box is always a cube
    void boundingBoxRendering(Vec3d & bmin, Vec3d & bmax); // assumes default expansion ratio (2.0)
    void boundingBox(double expansionRatio,Vec3d & bmin, Vec3d & bmax);
    void boundingBoxCube(double expansionRatio,Vec3d & bmin, Vec3d & bmax);
    double diameter(); // for a box with expansion ratio of 1

    // all locations are 0-indexed
    inline Vec3d vertexPosition( int iPos ) const { return vertexPositions_[iPos]; }
    inline Vec3d textureCoordinate( int iPos ) const { return textureCoords_[iPos]; }
    inline Vec3d normal( int iPos ) const { return normals_[iPos]; }
    inline Material material(unsigned int index) { return materials_[index]; }
    inline Material * materialHandle(unsigned int index) { return &materials_[index]; }

    inline size_t numVertices() { return vertexPositions_.size(); }
    unsigned int numFaces(); // total number of faces in all groups
    inline size_t numNormals() { return normals_.size(); }
    inline size_t numTextureCoordinates() { return textureCoords_.size(); }
    inline size_t numGroups() { return groups_.size(); }
    inline size_t numMaterials() { return materials_.size(); }
    int removeDuplicatedMaterials();

    inline void addMaterial( Material material ) { materials_.push_back(material); }
    inline void addMaterial( std::string name, Vec3d Ka, Vec3d Kd, Vec3d Ks, double shininess ) { materials_.push_back(Material(name,Ka,Kd,Ks,shininess));}
    inline void addGroup( Group & group ) { groups_.push_back(group);}
    inline void addGroup( std::string name ) { groups_.push_back(Group(name));}
    inline void addVertexPosition (Vec3d pos) { vertexPositions_.push_back(pos); }
    inline void addVertexNormal (Vec3d normal) { normals_.push_back(normal); }
    inline void addTextureCoordinate (Vec3d uv) { textureCoords_.push_back(uv); }
    inline void addFaceToGroup(Face face, unsigned int group) { groups_[group].addFace(face); }

    inline void setVertexPosition(int iPos, Vec3d pos) { vertexPositions_[iPos] = pos; }

    // returns a vector, specifying the surface area belonging to each vertex
    void computeSurfaceAreaPerVertex(); 
    inline double surfaceAreaPerVertex(unsigned int i) { return surfaceAreaPerVertex_[i]; }

    // these routines assume that the faces are oriented coherently
    void computePseudoNormals();
    inline Vec3d pseudoNormal(unsigned int i) { return pseudoNormals_[i]; }
    void setNormalsToPseudoNormals();

    // these routine assume that the faces are oriented coherently
    void computeEdgePseudoNormals();
    int edgePseudoNormal(unsigned int i, unsigned int j, Vec3d * pseudoNormal);

    // erases zero area faces from all groups
    int removeZeroAreaFaces();

    // scales the model uniformly, with center being the center of the scaling
    void uniformScale(Vec3d center, double factor);

    double computeMass_Triangles(std::vector<double> & groupDensities); // second argument gives surface density per each group
    Vec3d computeCenterOfMass_Vertices(); // of the vertices
    Vec3d computeCenterOfMass_Triangles(); // of the triangular surface
    Vec3d computeCenterOfMass_Triangles(std::vector<double> & groupDensities); // second argument gives surface density per each group

    void computeInertiaTensor_Triangles(double mass, double IT[6]); // of the triangular surface, with respect to the center of mass, assumes uniform density
    void computeInertiaTensor_Triangles(std::vector<double> & groupDensities, double IT[6]);

    inline Group * groupHandle(unsigned int group) { return &(groups_[group]); }

    // returns the integer position of a specified vertex
    inline int vertexIndex(unsigned int group, unsigned int faceIndex, unsigned int vert) { return groups_[group].face(faceIndex).vertex(vert).position(); }

    // writes the obj to a disk file (including materials to filename.mtl if so requested)
    void writeToFile(const std::string & filename, int outputMaterials = 0);

    // only writes geometry (not materials)
    void writeToSTLFile(const std::string & filename);

    // Abaqus format
    // only writes geometry (not materials)
    void writeToAbqFile(const std::string & filename);

    Vec3d faceCentroid(Face & face);
    double faceSurfaceArea(Face & face); // of a single face
    void faceSurfaceAreas(std::vector<double> & surfaceAreas);
    // warning: this only cares about the first triangle in a face (assumes planar face):
    Vec3d faceNormal(Face & face);

    double computeSurfaceArea(); // of the whole mesh
    void computeSurfaceAreaPerGroup(std::vector<double> & surfaceAreas);
    // computes masses belonging to each vertex (different groups can have different surface mass densities)
    void computeMassPerVertex(std::vector<double> & groupSurfaceMassDensity, std::vector<double> & masses);

    // computes the 3D volume enclosed by the orientable surface
    // assumes a triangle mesh
    double computeVolume();

    // generates the normals by averaging normals for adjacent faces
    // any pre-specified normals are overwritten by these new normals
    // does not assume a triangular mesh
    void setNormalsToAverageFaceNormals(); 
    // sets normals face normals
    void setNormalsToFaceNormals();       

    double minEdgeLength(); // computes minimum edge length in the mesh
    double medianEdgeLength(); // computes median edge length in the mesh
    double maxEdgeLength(); // computes maximum edge length in the mesh

    double minEdgeLength(int * vtxa, int * vtxb); // also returns 0-indexed vertices achieving min
    double maxEdgeLength(int * vtxa, int * vtxb); // also returns 0-indexed vertices achieving max

    unsigned int maxFaceDegree();

    // finds nearest mesh vertex to query position pos (using exhaustive search)
    unsigned int nearestVertex(Vec3d queryPos, double & dist);

    void getCentroids(std::vector<Vec3d> & centroids);
    void interpolateToCentroids(std::vector<double> & nodalData, std::vector<double> & centroidData);
    void interpolateToCentroids(std::vector<Vec3d> & nodalData, std::vector<Vec3d> & centroidData);

    void initSurfaceSampling();
    Vec3d sampleSurfacePosition(double sample); // sample should be between 0 and 1
    void dumpGeometry(int * numVertices, double ** vertices, 
        int * numTriangles = NULL, int ** triangles = NULL);

    void computeBoundingBox(); // sets the following:

  protected:
    std::vector< Material > materials_;
    std::vector< Group > groups_;
    std::vector< Vec3d > vertexPositions_;
    std::vector< Vec3d > textureCoords_;
    std::vector< Vec3d > normals_;
    //std::string filename_;

    // computes internal bounding box (rectangular)
    double diameter_;
    Vec3d bmin_,bmax_;
    Vec3d center_, cubeHalf_;

    std::vector<double> surfaceAreaPerVertex_;
    std::vector<Vec3d> pseudoNormals_;

    // index assumes that the first int is smaller than the second
    std::map< std::pair<unsigned int,unsigned int>, Vec3d > edgePseudoNormals_;

    // inertia tensor around the origin, assuming a triangle mass of 1
    void computeSpecificInertiaTensor(Vec3d & v0, Vec3d & v1, Vec3d & v2, double t[6]);
    void dirname(char * path, char * result);

    void stripBlanks(char * s);

    void parseMaterials(char * objFilename, char * materialFilename);

    std::vector<std::pair<double,Face*> > surfaceSamplingAreas;

    friend class MeshSmoothing;

};


class ObjFileOrientable //: public ObjFile
{
  public:
    // calls ObjFile
    // if generateHalfEdges flag is on, it also generates the half edges (otherwise class not fully initialized)
    // if numOrientationFlips is not NULL, returns the number of edges that were flipped to orient the surface coherently
    ObjFileOrientable( const std::string& filename, int generateHalfEdges=1, int * numOrientationFlips = NULL);

    ObjFileOrientable( ObjFile * objFile, int generateHalfEdges=1, int * numOrientationFlips = NULL);

    ~ObjFileOrientable();

    // makes an empty structure
    ObjFileOrientable() {}

    class HalfEdge
    {
      public:
        explicit HalfEdge( const unsigned int & position_g,
            const unsigned int& startVertex_g, const unsigned int& endVertex_g, 
            const unsigned int& startV_g, const unsigned int& endV_g, 
            const unsigned int groupID_g, const unsigned int& face_g,
            const int opposite_g, // value of -1 denotes boundary edge
            const unsigned int next_g)
          :	position_(position_g),
          startVertex_(startVertex_g), endVertex_(endVertex_g), 
          startV_(startV_g), endV_(endV_g),
          groupID_(groupID_g),face_(face_g),
          opposite_(opposite_g), next_(next_g) {}

        // accessors for getting global edge position
        unsigned int position() { return position_;}

        // accessors for starting and ending vertices of the edge (global indexing)
        unsigned int startVertex() { return startVertex_; }
        unsigned int endVertex() { return endVertex_; }

        // accessors for starting and ending vertices of the edge (local indexing on the local face)
        unsigned int startV() { return startV_; }
        unsigned int endV() { return endV_; }

        // accessors for the face on the left of the edge
        unsigned int groupID() { return groupID_; }
        unsigned int face() { return face_; }

        // accessors for opposite and next edges
        int opposite() { return opposite_; }
        unsigned int next() { return next_; }

        // mutator for opposite and next edges
        void setOpposite(int opposite_g) { opposite_ = opposite_g; }
        void setNext(unsigned int next_g) { next_ = next_g; }

        // is this edge a boundary edge
        bool isBoundary() { return (opposite_ == -1);}

        void flipOrientation(); // flips orientation of the edge (careful, structure not coherent any more now)

        bool operator== (HalfEdge & halfEdge2) { return (position_ == halfEdge2.position_);} 
        bool operator!= (HalfEdge & halfEdge2) { return (position_ != halfEdge2.position_);} 


      protected:
        unsigned int position_; // the global position of the half-edge in the data structure
        unsigned int startVertex_, endVertex_; // global vertex indices
        unsigned int startV_, endV_; // local vertex indices, on the face
        unsigned int groupID_;
        unsigned int face_;
        int opposite_;
        unsigned int next_;


    };



  protected:
    ObjFile * objFile;
    void Init(int generateHalfEdges, int * numOrientationFlips_ );

    std::vector< HalfEdge > halfEdges_;
    std::vector< int > boundaryEdges_;
    unsigned int connectedComponents; // the total number of connected components of the mesh

    bool hasBoundary_; // does the surface have boundary

    std::vector<int> edgesAtVertices_; // for every vertex, contains one half-edge emanating out of it
    std::vector<std::vector<int> > edgesAtFaces_; // for every face, contains one half-edge on this face

    void determineIfSurfaceHasBoundary();

  public:
    size_t numHalfEdges() { return halfEdges_.size(); }
    HalfEdge & halfEdge(unsigned int i) { return halfEdges_[i]; }

    // this function is mostly called internally, but can sometimes also be called from the outside
    int GenerateHalfEdgeDataStructure(); // generates the whole datastructure, assuming the base objFile class has been initialized
    // returns the number of edges that were flipped to orient the surface coherently
    // this will be zero if the input mesh already is oriented coherently

    void CopyHalfEdgeTopologyFrom(ObjFileOrientable * source); // makes the half-edge topological info equal to that of source

    // returns the opposite halfedge to the given half-edge
    // this function will fail for boundary edges, should always check first with isBoundary()
    HalfEdge & edgeOpposite ( HalfEdge & halfedge ) { return halfEdges_[halfedge.opposite()]; } 

    // returns the next halfedge to the given half-edge
    HalfEdge & edgeNext ( HalfEdge & halfedge ) { return halfEdges_[halfedge.next()]; } 

    // returns the previous halfedge to the given half-edge
    // does so by looping around the face (pointers to previous edges are not explicitly stored), so this is slower than edgeNext
    HalfEdge & edgePrevious ( HalfEdge & halfedge );

    // loops around the vertex (vertex is defined as the ending position of the halfedge)
    // consists of taking the next edge, then taking the opposite edge
    // if boundary edge encountered, can't take the opposite edge; it this case flag=1 is returned and the edge returned is the boundary edge pointing away from the vertex
    // if taking the opposite edge is possible, the returned edge points into the vertex and flag is set to 0
    // this is effectively looping in the clockwise (negative) orientation
    HalfEdge & loopVertex(HalfEdge & halfedge, int & flag); 

    // returns the the group that contains the given half-edge
    ObjFile::Group edgeGroup ( HalfEdge & halfedge ) { return *(objFile->groupHandle(halfedge.groupID())); } 

    // returns the face to the left of the given half-edge
    ObjFile::Face edgeFace ( HalfEdge & halfedge ) 
    { ObjFile::Group * groupHandle = objFile->groupHandle(halfedge.groupID());
      return groupHandle->face(halfedge.face()); }

      // returns the starting vertex of the given half-edge
      ObjFile::Vertex edgeStartVertex ( HalfEdge & halfedge ) 
      { ObjFile::Group * groupHandle = objFile->groupHandle(halfedge.groupID());
        return (groupHandle->face(halfedge.face())).vertex(halfedge.startV()); }

        // returns the ending vertex of the given half-edge
        ObjFile::Vertex edgeEndVertex ( HalfEdge & halfedge ) 
        { ObjFile::Group * groupHandle = objFile->groupHandle(halfedge.groupID());
          return (groupHandle->face(halfedge.face())).vertex(halfedge.endV()); }

          unsigned int nConnectedComponents() { return connectedComponents;}

          bool isIsolatedVertex( unsigned int vertex ) { return (edgesAtVertices_[vertex] == -1); } // returns true if vertex is isolated

          bool isBoundaryVertex( unsigned int vertex ) { return halfEdges_[edgesAtVertices_[vertex]].isBoundary(); } // returns true if vertex is a mesh boundary vertex

          void PrintHalfEdges(); // prints the half-edges out

          // returns some halfedge emanating out of a given vertex (returns always the same edge)
          // in case vertex is a boundary vertex, it will return the edge such that there is no clockwise edge to the given edge around the given vertex
          // this function will fail for isolated vertices; should always check first with isIsolatedVertex()
          HalfEdge & edgeAtVertex( unsigned int vertex ) { return halfEdges_[edgesAtVertices_[vertex]]; }

          // returns some halfedge on the given face (returns always the same edge)
          HalfEdge & edgeAtFace( unsigned int groupID, unsigned int faceID ) { return halfEdges_[edgesAtFaces_[groupID][faceID]]; }

          // returns true if surface has boundary and false if it closed
          bool hasBoundary() { return hasBoundary_;}

          size_t numBoundaryEdges() { return boundaryEdges_.size(); }
          int boundaryEdge(int i) { return boundaryEdges_[i]; }

          int internalAllocation;
};

inline bool my_isnan(double x)
{
  return (x != x);
}


#endif
