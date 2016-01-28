#include <vector>
#include <queue>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string.h>
#include <assert.h>
using namespace std;

#include "float.h"
#include "objfile.h"

// Half-edge data structure derived class implementation by Jernej Barbic <barbic@cs.cmu.edu>

using std::ifstream;
using std::string;



// removes all whitespace characters from string s
void ObjFile::stripBlanks(char * s)
{

  char * w = s;
  while (*w != '\0')
  {
    while (*w == ' ') // erase blank
    {
      char * u = w;
      while (*u != '\0') // shift everything left one char
      {
        *u = *(u+1);
        u++;
      }
    }
    w++;
  }
}


// shrinks all whitespace to a single space character; also removes any whitespace before end of string
void convertWhitespaceToSingleBlanks(char * s)
{
  char * w = s;
  while (*w != '\0')
  {
    while ((*w == ' ') && ((*(w+1) == '\0') || (*(w+1) == ' ') )) // erase consecutive blanks or end-of-string blanks
    {
      char * u = w;
      while (*u != '\0') // shift everything left one char
      {
        *u = *(u+1);
        u++;
      }
    }
    w++;
  }
}

ObjFile::ObjFile( const std::string& filename, int verbose )
  //: filename_(filename)
{
  std::string filename_ = filename;

  unsigned int numFaces = 0;

  const int maxline = 1000;
  std::ifstream ifs(filename.c_str());
  char line[maxline];

  unsigned int currentGroup=0;
  unsigned int ignoreCounter=0;

  //unsigned int  currentMaterialIndex = 0;
  addMaterial("default",Vec3d(1,1,1),Vec3d(1,1,1),Vec3d(1,1,1),0);

  if (verbose)
    std::cout << "Parsing .obj file '" << filename << "'." << std::endl;

  if( !ifs )
  {
    std::string message = "Couldn't open .obj file '";
    message.append( filename );
    message.append( "'" );
    throw ObjFileException( message );
  }

  int lineNum = 0;

  while( ifs )
  {
    ++lineNum;

    ifs.getline( line, maxline );

    convertWhitespaceToSingleBlanks(line);

    char command = line[0];

    if( strncmp(line,"v ",2) == 0 ) // vertex
    {
      //std::cout << "v " ;
      Vec3d pos;
      double x,y,z;
      if (sscanf(line,"v %lf %lf %lf\n",&x,&y,&z) < 3)
      {
        throw ObjFileException("Invalid vertex", filename, lineNum);
      }
      pos = Vec3d(x,y,z);
      vertexPositions_.push_back( pos );
      //cout << vertexPositions_.size() << ":" << vertexPositions_.back() << endl;
    }
    else if( strncmp(line,"vn ",3) == 0)
    {
      //std::cout << "vn " ;
      Vec3d normal;
      double x,y,z;
      if (sscanf(line,"vn %lf %lf %lf\n",&x,&y,&z) < 3)
      {
        throw ObjFileException("Invalid normal", filename, lineNum);
      }
      normal = Vec3d(x,y,z);
      normals_.push_back( normal );
    }
    else if( strncmp(line,"vt ",3) == 0 )
    {
      //std::cout << "vt " ;
      Vec3d tex = vl_0;	// w defaults to 0
      double x,y;
      if (sscanf(line,"vt %lf %lf\n",&x,&y) < 2)
      {
        throw ObjFileException("Invalid texture coordinate", filename, lineNum);
      }
      tex = Vec3d(x,y,0);
      textureCoords_.push_back( tex );
    }
#if 0
    else if( strncmp(line,"g ",2) == 0 )
    {
      char s[4096];
      if (sscanf(line,"g %s\n",s) < 1)
        if (verbose)
          cout << "Warning:  Empty group name encountered: " << filename << " " << lineNum << endl;
      std::string name = s;

      // check if this group already exists
      bool groupFound = false;
      unsigned int counter = 0;
      for( std::vector< Group >::const_iterator itr = groups_.begin(); itr != groups_.end(); ++itr )
      {
        if( itr->name() == name )
        {
          currentGroup = counter;
          groupFound = true;
          break;
        }
        counter++;
      }
      if (!groupFound)
      {
        groups_.push_back( Group(name,currentMaterialIndex) );
        currentGroup = groups_.size()-1;
      }
    }
#endif
    else if ( strncmp(line,"g",1) == 0 )
    {
      // Ignore groups
      continue;
    }
    else if(( strncmp(line,"f ",2) == 0 ) || ( strncmp(line,"fo ",3) == 0 ))
    {
      char * faceLine = &line[2];
      if  ( strncmp(line,"fo",2) == 0 )
        faceLine = &line[3];

      //std::cout << "f " ;
      if( groups_.empty() )
      {
        groups_.push_back( Group("default") );
        currentGroup = 0;
      }

      Face face;

      // the faceLine string now looks like the following:
      //   vertex1 vertex2 ... vertexn
      // where vertexi is v/t/n, v//n, v/t, or v

      char * curPos = faceLine;
      while( *curPos != '\0' )
      {
        // seek for next whitespace or eof
        char * tokenEnd = curPos;
        while ((*tokenEnd != ' ') && (*tokenEnd != '\0'))
          tokenEnd++;

        bool whiteSpace = false;
        if (*tokenEnd == ' ')
        {
          *tokenEnd = '\0';
          whiteSpace = true;
        }

        unsigned int pos;
        unsigned int nor;
        unsigned int tex;
        std::pair< bool, unsigned int > texPos;
        std::pair< bool, unsigned int > normal;

        // now, parse curPos
        if (strstr(curPos,"//") != NULL) 
        {
          // v//n
          if (sscanf(curPos,"%u//%u",&pos,&nor) < 2)
          {
            throw ObjFileException( "Invalid face", filename, lineNum);
          }
          texPos = make_pair(false,0);
          normal = make_pair(true,nor);
        }
        else 
        {
          if (sscanf(curPos,"%u/%u/%u",&pos,&tex,&nor) != 3)
          {
            if (strstr(curPos,"/") != NULL)
            {
              // v/t
              if (sscanf(curPos,"%u/%u",&pos,&tex) == 2)
              {
                texPos = make_pair(true,tex);
                normal = make_pair(false,0);
              }
              else
              {
                throw ObjFileException( "Invalid face", filename, lineNum);
              }
            }
            else
            { // v
              if (sscanf(curPos,"%u",&pos) == 1)
              {
                texPos = make_pair(false,0);
                normal = make_pair(false,0);
              }
              else
              {
                throw ObjFileException( "Invalid face", filename, lineNum);
              }
            }
          }
          else
          { // v/t/n
            texPos = make_pair(true,tex);
            normal = make_pair(true,nor);
          }
        }

        // decrease indices to make them 0-indexed
        pos--;
        if (texPos.first)
          texPos.second--;
        if (normal.first)
          normal.second--;

        face.addVertex( Vertex( pos, texPos, normal ) );

        if(whiteSpace)
        {
          *tokenEnd = ' ';
          curPos = tokenEnd + 1; 
        }
        else
          curPos = tokenEnd;
      }

      ++numFaces;
      groups_[currentGroup].addFace( face );
    }
    else if  ((strncmp(line,"#",1) == 0 ) ||  (strncmp(line,"\0",1) == 0 ))
    { // ignore comment lines and empty lines
    }
    else if (strncmp(line,"usemtl",6) == 0)
    {
      // Ignore material lines
      continue;
#if 0
      // switch to a new material
      bool materialFound = false;
      unsigned int counter = 0;
      char * materialName = &line[7];
      for( std::vector< Material >::const_iterator itr = materials_.begin(); itr != materials_.end(); ++itr )
      {
        if( itr->name() == string(materialName))
        {
          currentMaterialIndex = counter;

          // update current group
          if( groups_.empty() )
          {
            groups_.push_back( Group("default") );
            currentGroup = 0;
          }

          groups_[currentGroup].setMaterialIndex(currentMaterialIndex);
          materialFound = true;
          break;
        }
        counter++;
      }
      if (!materialFound)
      {
        char msg[4096];
        sprintf(msg,"Material %s does not exist.\n",materialName);
        throw ObjFileException(msg);
      }
#endif
    }
    else if( strncmp(line,"mtllib",6) == 0 )
    {
      char mtlFilename[4096];
      strcpy(mtlFilename,filename.c_str());
      parseMaterials(mtlFilename, &line[7]);
    }
    else if (strncmp(line,"s ",2) == 0 )
    {
      //std::cout << command << " ";
      // ignore for now
      if (ignoreCounter < 5)
      {
        if (verbose)
          std::cout << "Warning: ignoring '" << command << "' line" << std::endl;
        ignoreCounter++;
      }
      if (ignoreCounter == 5)
      {
        if (verbose)
          std::cout << "(suppressing further output of ignored lines)" << std::endl;
        ignoreCounter++;
      }
    }
    else
    {
      //std::cout << "invalid ";
      std::ostringstream msg;
      msg << "Invalid line in .obj file '" << filename << "': " << line;
      throw ObjFileException(msg.str(), filename, lineNum);
    }
  }

  computeBoundingBox();

  // statistics
  if (verbose)
  {
    std::cout << "Parsed obj file '" << filename << "'; statistics:" << std::endl;
    std::cout << "   " << groups_.size() << " groups," << std::endl;
    std::cout << "   " << numFaces << " faces," << std::endl;
    std::cout << "   " << vertexPositions_.size() << " vertices," << std::endl;
    std::cout << "   " << normals_.size() << " normals, " << std::endl;
    std::cout << "   " << textureCoords_.size() << " texture coordinates, " << std::endl;
  }
}

std::vector<std::string> ObjFile::groupNames() const
{
  std::vector<std::string> result;
  result.reserve( groups_.size() );
  for( std::vector<Group>::const_iterator groupItr = groups_.begin();
      groupItr != groups_.end();
      ++groupItr )
  {
    result.push_back( groupItr->name() );
  }

  return result;
}

ObjFile::Group ObjFile::group( std::string name ) const
{
  for( std::vector< Group >::const_iterator itr = groups_.begin();
      itr != groups_.end();
      ++itr )
  {
    if( itr->name() == name )
      return *itr;
  }

  std::ostringstream oss;
  oss << "Invalid group name: '" << name << "'.";
  throw ObjFileException( oss.str() );
}

void ObjFile::PrintInfo() const
{
  typedef std::vector<std::string> SVec;
  SVec groupNames1 = groupNames();
  for( SVec::const_iterator nameItr = groupNames1.begin();
      nameItr != groupNames1.end();
      ++nameItr )
  {
    std::cout << "Found obj group '" << *nameItr << std::endl;
    ObjFile::Group group1 = group( *nameItr ); // retrieve group named *nameItr, and store it into "group"
    std::cout << "Iterating through group faces..." << std::endl;
    for( unsigned int iFace = 0; iFace < group1.faceCount(); ++iFace )
    {
      ObjFile::Face face = group1.face(iFace); // get face number iFace
      if( face.vertexCount() == 3 )
        std::cout << "  found triangle ";
      else if( face.vertexCount() == 4 )
        std::cout << "  found quadrilateral ";
      else
        std::cout << "  found " << face.vertexCount() << "-gon ";

      // Since the vertex positions are unique within the files, we can
      // use these to cross-index the polygons.
      for( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {
        if( iVertex != 0 )
          std::cout << " -> ";
        std::cout << face.vertex( iVertex ).position(); // print out integer indices of the vertices
      }
      std::cout << std::endl;

      // Now we will retrieve positions, normals, and texture coordinates of the
      // files by indexing into the global vertex namespace.
      for( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {
        ObjFile::Vertex vertex = face.vertex( iVertex );
        std::cout << "    vertex " << iVertex << "; " << std::endl;
        std::cout << "      position = "
          << vertexPosition( vertex.position() ) << ";" << std::endl;
        if( vertex.hasNormal() )
          std::cout << "      normal = " << normal( vertex.normal() ) << ";" << std::endl;
        if( vertex.hasTextureCoordinate() )
          std::cout << "      texture coordinate = " <<
            textureCoordinate( vertex.textureCoordinate() ) << ";" << std::endl;
      }
    }
  }

}

bool ObjFile::isTriangularMesh()
{
  for(unsigned int i=0; i < groups_.size(); i++) // over all groups
    for (unsigned int j=0; j < groups_[i].faceCount(); j++) // over all faces
    {
      if (groups_[i].face(j).vertexCount() != 3)
        return false;
    }
  return true;
}

bool ObjFile::isQuadrilateralMesh()
{
  for(unsigned int i=0; i < groups_.size(); i++) // over all groups
    for (unsigned int j=0; j < groups_[i].faceCount(); j++) // over all faces
    {
      if (groups_[i].face(j).vertexCount() != 4)
        return false;
    }
  return true;
}


unsigned int ObjFile::maxFaceDegree()
{
  unsigned int maxDegree = 0;
  for(unsigned int i=0; i < groups_.size(); i++) // over all groups
    for (unsigned int j=0; j < groups_[i].faceCount(); j++) // over all faces
    {
      Face * face = groups_[i].faceHandle(j);
      if (face->vertexCount() > maxDegree)
        maxDegree = face->vertexCount();
    }

  return maxDegree;
}

void ObjFile::triangulate()
{
  for(unsigned int i=0; i < groups_.size(); i++) // over all groups
    for (unsigned int j=0; j < groups_[i].faceCount(); j++) // over all faces
    {
      Face * face = groups_[i].faceHandle(j);
      if (face->vertexCount() < 3)
      {
        printf("Warning: encountered a face with fewer than 3 vertices.\n");
      }

      unsigned int faceDegree = face->vertexCount();

      if (faceDegree > 3)
      {
        // triangulate the face

        // get the vertices:
        vector<Vertex> vertices;
        for(unsigned int k=0; k<face->vertexCount(); k++)
          vertices.push_back(face->vertex(k));

        Face newFace;
        newFace.addVertex(vertices[0]);
        newFace.addVertex(vertices[1]);
        newFace.addVertex(vertices[2]);

        // overwrite old face
        *face = newFace;

        for(unsigned int k=2; k<faceDegree-1; k++)
        {
          // tesselate the remainder of the old face
          Face newFace;
          newFace.addVertex(vertices[0]);
          newFace.addVertex(vertices[k]);
          newFace.addVertex(vertices[k+1]);
          groups_[i].addFace(newFace);
        }

      }
    }
}


void ObjFile::computeBoundingBox()
{
  bmin_ =  vertexPositions_[0]; 
  bmax_ =  vertexPositions_[0]; 

  for(unsigned int i=1; i < vertexPositions_.size(); i++) // over all vertices
  {
    Vec3d p = vertexPositions_[i]; 

    if (p[0] < bmin_[0])
      bmin_[0] = p[0];
    if (p[0] > bmax_[0])
      bmax_[0] = p[0];

    if (p[1] < bmin_[1])
      bmin_[1] = p[1];
    if (p[1] > bmax_[1])
      bmax_[1] = p[1];

    if (p[2] < bmin_[2])
      bmin_[2] = p[2];
    if (p[2] > bmax_[2])
      bmax_[2] = p[2];
  }

  center_ = 0.5 * (bmin_ + bmax_);

  Vec3d halfside_ = 0.5* (bmax_ - bmin_);

  /*
     double maxHalf = halfside_[0];
     if (halfside_[1] > maxHalf)
     maxHalf = halfside_[1];
     if (halfside_[2] > maxHalf)
     maxHalf = halfside_[2];
     cubeHalf_ = Vec3d(maxHalf,maxHalf,maxHalf);
   */

  cubeHalf_ = halfside_;
  diameter_ = len(bmin_-bmax_);
}

void ObjFile::boundingBoxRendering(Vec3d & bmin, Vec3d & bmax)
{
  boundingBox(2.0, bmin, bmax);
}

void ObjFile::boundingBox(double expansionRatio, Vec3d & bmin, Vec3d & bmax)
{
  bmin = center_ - expansionRatio * cubeHalf_;
  bmax = center_ + expansionRatio * cubeHalf_;
}

void ObjFile::boundingBoxCube(double expansionRatio, Vec3d & bmin, Vec3d & bmax)
{
  double maxHalf = cubeHalf_[0];

  if (cubeHalf_[1] > maxHalf)
    maxHalf = cubeHalf_[1];

  if (cubeHalf_[2] > maxHalf)
    maxHalf = cubeHalf_[2];

  Vec3d cubeHalfCube_ = Vec3d(maxHalf,maxHalf,maxHalf);

  bmin = center_ - expansionRatio * cubeHalfCube_;
  bmax = center_ + expansionRatio * cubeHalfCube_;
}

double ObjFile::diameter()
{
  return diameter_;
}

void ObjFileOrientable::HalfEdge::flipOrientation()
{
  int buffer = startVertex_;
  startVertex_ = endVertex_;
  endVertex_ = buffer;
}

ObjFileOrientable::~ObjFileOrientable()
{
  if (internalAllocation)
    delete(objFile);
}

ObjFileOrientable::ObjFileOrientable( const std::string& filename, 
    int generateHalfEdges, int * numOrientationFlips_ ) //: ObjFile(filename)
{
  internalAllocation = 1;
  objFile = new ObjFile(filename);

  Init(generateHalfEdges, numOrientationFlips_);
}

ObjFileOrientable::ObjFileOrientable( ObjFile * objFile, 
    int generateHalfEdges, int * numOrientationFlips_ ) //: ObjFile(filename)
{
  internalAllocation = 0;
  this->objFile = objFile;

  Init(generateHalfEdges, numOrientationFlips_);
}

void ObjFileOrientable::Init(int generateHalfEdges, int * numOrientationFlips_ )
{
  if (generateHalfEdges)
  {
    int numOrientationFlips = GenerateHalfEdgeDataStructure();

    if (numOrientationFlips_ != NULL)
      *numOrientationFlips_ = numOrientationFlips;
  }
}


void ObjFileOrientable::PrintHalfEdges()
{
  for (unsigned int i=0; i<halfEdges_.size(); i++)
  {
    cout << "Half edge "<< i << " :" << endl;
    cout << "  Opposite edge: " << halfEdges_[i].opposite() << endl;
    cout << "  Next edge: " << halfEdges_[i].next() << endl;
    cout << "  Group: " << halfEdges_[i].groupID() << endl;
    cout << "  Face: " << halfEdges_[i].face() << endl;
    cout << "  Start vertex: " << halfEdges_[i].startVertex() << endl;
    cout << "  End vertex: " << halfEdges_[i].endVertex() << endl;
    cout << "  Start vertex (local): " << halfEdges_[i].startV() << endl;
    cout << "  End vertex (local): " << halfEdges_[i].endV() << endl;
    cout << "  Is boundary: " << halfEdges_[i].isBoundary() << endl;
  }
}

// returns the number of edges flipped
int ObjFileOrientable::GenerateHalfEdgeDataStructure()
{

  std::cout << "Building the half edge data structure..." << std::endl;

  // Step 1: iterate over all faces
  // for each face, add all the edges onto the list of half-edges

  std::cout << "Step 1: Generating the list of half edges..." << std::endl;

  //typedef std::vector<ObjFile::Group> SGroup;

  int coutCounter = 0;

  for(unsigned int i = 0; i < objFile->numGroups(); i++ )
  {
    ObjFile::Group * currentGroup = objFile->groupHandle(i);

    std::cout << "  Processing obj group '" << currentGroup->name() << std::endl;
    std::cout << "  Iterating through group faces..." << std::endl;

    for( unsigned int iFace = 0; iFace < currentGroup->faceCount(); ++iFace )
    {
      ObjFile::Face face = currentGroup->face(iFace); // get face whose number is iFace

      if (coutCounter < 100)
      {
        std::cout << face.vertexCount() ;
        coutCounter++;
      }
      if (coutCounter == 100)
      {
        cout << "...[and more]";
        coutCounter++;
      }

      unsigned int edgesSoFar = halfEdges_.size();

      for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {
        // create a half edge for each edge, store -1 for half-edge adjacent edge for now
        // index vertices starting from 0
        int nextEdge = edgesSoFar + ((iVertex + 1) % face.vertexCount());
        HalfEdge halfEdge(edgesSoFar + iVertex, face.vertex(iVertex).position(), face.vertex((iVertex + 1) % face.vertexCount()).position(),
            iVertex, (iVertex + 1) % face.vertexCount(),	
            i, iFace, -1, nextEdge); 

        halfEdges_.push_back(halfEdge);
      }
    }
    std::cout << std::endl;
  }

  /*  
      for (unsigned int i=0; i<halfEdges_.size(); i++)
      {
      cout << "Half edge "<< i << " :" << endl;
      cout << "  Opposite edge: " << halfEdges_[i].opposite() << endl; 
      cout << "  Next edge: " << halfEdges_[i].next() << endl; 
      cout << "  Group: " << halfEdges_[i].groupID() << endl; 
      cout << "  Face: " << halfEdges_[i].face() << endl; 
      cout << "  Start vertex: " << halfEdges_[i].startVertex() << endl; 
      cout << "  End vertex: " << halfEdges_[i].endVertex() << endl; 
      cout << "  Start vertex (local): " << halfEdges_[i].startV() << endl; 
      cout << "  End vertex (local): " << halfEdges_[i].endV() << endl; 
      cout << "  Is boundary: " << halfEdges_[i].isBoundary() << endl; 
      }
   */

  // Step 2: build correspondence among half-dges
  // for each half-edge, search for the opposite half-edge, if it exists

  std::cout << "Step 2: Building correspondence among the half-edges..." << std::endl;
  std::cout << "Boundary edges: ";

  // insert all edges into a binary tree

  typedef std::multimap<std::pair< unsigned int, unsigned int > , unsigned int> BinaryTree;
  BinaryTree edges;

  for (unsigned int i=0; i < halfEdges_.size(); i++)
  {

    int vertex1 = halfEdges_[i].startVertex();
    int vertex2 = halfEdges_[i].endVertex();

    if (vertex1 == vertex2)
    {
      std::cout << "Error: encountered a degenerated edge with equal starting and ending vertex." << std::endl;
      std::cout << "  Group:" << halfEdges_[i].groupID() << "  Face #: " << halfEdges_[i].face() << "Vertex ID: " << vertex1 << std::endl;
      exit(1);
    }

    if (vertex1 > vertex2) // swap
    {
      int buffer = vertex1;
      vertex1 = vertex2;
      vertex2 = buffer;
    }

    std::pair<unsigned int, unsigned int> vertices(vertex1,vertex2);
    edges.insert(std::make_pair(vertices,i)); 
  }  

  // retrieve one by one and build correspondence
  for (unsigned int i=0; i < halfEdges_.size(); i++)
  {
    int vertex1 = halfEdges_[i].startVertex();
    int vertex2 = halfEdges_[i].endVertex();

    if (vertex1 > vertex2) // swap
    {
      int buffer = vertex1;
      vertex1 = vertex2;
      vertex2 = buffer;
    }

    std::pair<unsigned int, unsigned int> vertices(vertex1,vertex2);

    // search for the edge

    int hits = 0;
    int candidates = 0;
    BinaryTree::iterator pos;
    for (pos = edges.lower_bound(vertices); 
        pos != edges.upper_bound(vertices);
        ++pos)
    {
      candidates++;
      // check if we found ourselves
      if (pos->second != i) 
      { // not ourselves
        halfEdges_[i].setOpposite(pos->second);
        hits++;
      }
    }

    if (candidates >= 3) 
    {
      std::cout << "Error: encountered an edge that appears in more than two triangles. Geometry is non-manifold. Exiting." << std::endl;
      int faceNum = halfEdges_[i].face();
      std::cout << "  Group:" << halfEdges_[i].groupID() << std::endl << "  Face #: " << faceNum << std::endl;
      std::cout << "  Edge: " << vertex1 << " " << vertex2 << std::endl;
      std::cout << "  Vertices: " << objFile->vertexPosition(vertex1) << " " << objFile->vertexPosition(vertex2) << std::endl;
      exit(1);
    }

    if (hits == 0) // boundary edge
    {  
      std::cout << "B"; 
      boundaryEdges_.push_back(i);
    }

  }  

  std::cout << " " << boundaryEdges_.size() << std::endl;
  /*  
      for (unsigned int i=0; i<halfEdges_.size(); i++)
      {
      cout << "Half edge "<< i << " :" << endl;
      cout << "  Opposite edge: " << halfEdges_[i].opposite() << endl; 
      cout << "  Next edge: " << halfEdges_[i].next() << endl; 
      cout << "  Group: " << halfEdges_[i].groupID() << endl; 
      cout << "  Face: " << halfEdges_[i].face() << endl; 
      cout << "  Start vertex: " << halfEdges_[i].startVertex() << endl; 
      cout << "  End vertex: " << halfEdges_[i].endVertex() << endl; 
      cout << "  Start vertex (local): " << halfEdges_[i].startV() << endl; 
      cout << "  End vertex (local): " << halfEdges_[i].endV() << endl; 
      cout << "  Is boundary: " << halfEdges_[i].isBoundary() << endl; 
      }
   */
  // now, each half-edge knows its mirror edge, but orientations of faces might be inconsistent

  // orient all half-edges correctly
  std::cout << "Step 3: Attempting to orient the faces coherently..." << std::endl;

  // build marks for all edges
  std::vector<int> marks;
  for (unsigned int i=0; i < halfEdges_.size(); i++)
    marks.push_back(0);

  // initialize queue
  std::priority_queue<int> queue;

  connectedComponents = 0;
  int numOrientationFlips = 0;

  while(1) // breakable
  {

    // find the first un-marked edge and queue it
    unsigned int unmarkedEdge = 0;
    for (unmarkedEdge = 0; unmarkedEdge < halfEdges_.size(); unmarkedEdge++)
    {
      if (marks[unmarkedEdge] == 0)
        break; // found an unmarked edge
    }

    if (unmarkedEdge == halfEdges_.size()) // no unmarked edge was found
    {
      break; // out of while; we are done
    }
    else
    {
      cout << "Found new connected component. Seed half-edge is: " << unmarkedEdge << endl;
      connectedComponents++;
      queue.push(unmarkedEdge);

      while(!queue.empty())
      {
        int edge = queue.top();
        queue.pop();

        //std::cout << "Retrieved edge from queue: " << edge << std::endl;
        //cout << "The edge is boundary: " << halfEdges_[edge].isBoundary() << endl;

        //std::cout << "Marking all the edges on this face: ";
        // first, mark all the edges on this face as visited
        int loop = edge;
        do 
        {
          marks[loop] = 1;
          //std::cout << loop << " ";
          loop = halfEdges_[loop].next();
        }
        while (loop != edge);
        //std::cout << std::endl; 

        // check if edge is consistent with the opposite edge orientation
        // careful: edge might be on the boundary
        // find a non-boundary edge on the same face (if it exists)
        //std::cout << "Seeking for a non-boundary edge on this face...";
        loop = edge;
        int exitFlag = 0;

        while ((halfEdges_[loop].isBoundary()) && (exitFlag == 0))
        {
          //cout << loop << " ";
          loop = halfEdges_[loop].next();
          if (loop == edge) // all edges are boundary
            exitFlag = 1;
        }


        if (exitFlag == 1) // all edges are boundary; this is an isolated face
        {
          //cout << "none found." << endl;
          continue; // no need to queue anything or flip anything, this was an isolated face
          // also, this case can only happen during the first iteration of the while loop, which will also be the last one
        }

        edge = loop; // now, edge is an halfedge not on a boundary

        //std::cout << "found non-boundary edge: " << edge << std::endl;
        //std::cout << "opposite edge is: " << halfEdges_[edge].opposite() << std::endl;


        bool orientationFlipNecessary = (halfEdges_[edge].startVertex() == (edgeOpposite(halfEdges_[edge])).startVertex());

        //std::cout << "Orientation flip necessary for this face: " << orientationFlipNecessary << std::endl;

        if (orientationFlipNecessary) 
        { // flip all edges along this face
          //cout << "Orientation flip" << endl;
          numOrientationFlips++;
          loop = edge;
          int cache = 0;
          do
          {
            int nextO = halfEdges_[loop].next();
            halfEdges_[loop].setNext(cache);
            cache = loop;
            halfEdges_[loop].flipOrientation(); // flip orientation
            loop = nextO;
          }
          while (loop != edge);
          halfEdges_[loop].setNext(cache);

          int groupID = halfEdges_[loop].groupID();
          int faceID = halfEdges_[loop].face();

          ObjFile::Group * currentGroup = objFile->groupHandle(groupID);
          currentGroup->face(faceID).displayVertices();
          currentGroup->reverseFace(faceID);
          currentGroup->face(faceID).displayVertices();
        }

        // check if new orientation is consistent eveywhere along the face
        // if not, surface is not orientable
        // at the same time, queue the opposite edges if they are not marked already
        loop = edge;
        do
        {
          if (!halfEdges_[loop].isBoundary()) // skip boundary edges
          {

            // if opposite unmarked, queue the opposite edge 
            if (marks[halfEdges_[loop].opposite()] == 0)
            {
              queue.push(halfEdges_[loop].opposite());
              marks[halfEdges_[loop].opposite()] = 1;
              //std::cout << "visiting edge: " << loop << " pushing opposite: " << halfEdges_[loop].opposite() << std::endl;
            }
            else
            { // opposite edge is marked as already visited
              // if orientation consistent, do nothing 
              // if orientation not consistent, surface is not orientable

              bool orientationConsistent = (halfEdges_[loop].startVertex() == (edgeOpposite(halfEdges_[loop])).endVertex()); 

              //std::cout << "visiting edge: " << loop << " opposite marked " << std::endl;

              if (!orientationConsistent)
              {
                std::cout << "Error: surface is non-orientable." << std::endl;
                exit(1);
              }
            }
          }
          loop = halfEdges_[loop].next();
        }
        while (loop != edge);

      }
    }
  } // end of while  

  printf("Consistent orientation generated. Performed %d orientation flips.\n", numOrientationFlips);


  //PrintHalfEdges();

  // step 4: for every vertex, find a half-edge emanating out of it
  std::cout << "Step 4: For every vertex, caching a half-edge emanating out of it..." << std::endl;
  std::cout << "        For every face, caching a half-edge on it..." << std::endl;

  for (unsigned int i=0; i< objFile->numVertices(); i++)
    edgesAtVertices_.push_back(-1); // value of -1 corresponds to no edge (i.e. isolated vertex)

  for (unsigned int i=0; i < halfEdges_.size(); i++)
  {
    //cout << i << " " << halfEdges_[i].startVertex() << " " << halfEdges_[i].endVertex() << endl;
    edgesAtVertices_[halfEdges_[i].endVertex()] = i;
  }

  // if vertex is on the boundary, rotate the edge until it is an incoming boundary edge
  // rotate edge until it is either on the boundary, or we come around to the same edge
  int numIsolatedVertices = 0;
  for (unsigned int i=0; i < objFile->numVertices(); i++)
  {
    if (isIsolatedVertex(i))
    {
      numIsolatedVertices++;
      continue;
    }
    HalfEdge * loop = &edgeAtVertex(i);
    HalfEdge * start = loop; 
    do
    {
      if (loop->isBoundary())
      {
        // set the edge to the current edge
        edgesAtVertices_[i] = loop->position();
        break;
      }
      loop = &edgePrevious(edgeOpposite(*loop));
    }
    while (*loop != *start);
    // if we came around, no need to change edgeAtVertices[i]
  } 

  if (numIsolatedVertices > 0)
    printf("Warning: mesh has %d isolated vertices.\n", numIsolatedVertices);

  // build the cache for faces, first reset to -1
  for (unsigned int i=0; i < objFile->numGroups(); i++)
  {
    ObjFile::Group * currentGroup = objFile->groupHandle(i);
    std::vector<int> dataForThisGroup;
    dataForThisGroup.clear();
    for (unsigned int j=0; j < currentGroup->faceCount(); j++)
    { 
      dataForThisGroup.push_back(-1);
    }
    edgesAtFaces_.push_back(dataForThisGroup);
  }   
  for (unsigned int i=0; i < halfEdges_.size(); i++)
    edgesAtFaces_[halfEdges_[i].groupID()][halfEdges_[i].face()] = i;

  // sanity check: none of the face entries should be -1
  for (unsigned int i=0; i < objFile->numGroups(); i++)
  {
    ObjFile::Group * currentGroup = objFile->groupHandle(i);
    for (unsigned int j=0; j < currentGroup->faceCount(); j++)
      if (edgesAtFaces_[i][j] == -1)
        cout << "Warning: face on group " << i << "(" << currentGroup->name() << "), position " << j << " has no edges." << endl;
  }

  determineIfSurfaceHasBoundary();

  // testing: previous edge capability
  /*
     cout << "Testing previous edge capability..." << endl;
     for (unsigned int i=0; i < halfEdges_.size(); i++)
     {
     cout << i << ": " << edgePrevious(halfEdges_[i]).position() << endl;
     }

  // testing: print out associated edges for every vertex
  for (unsigned int i=0; i < vertexPositions_.size(); i++)
  {
  cout << "Halfedge into vertex " << i << ": " << edgeAtVertex(i).position() << endl;
  }

  // testing: print out associated edges for every face
  for (unsigned int i=0; i < groups_.size(); i++)
  for (unsigned int j=0; j < groups_[i].faceCount(); j++)
  {
  cout << "Halfedge on face " << i << " " << j << ": " << edgeAtFace(i,j).position() << endl;
  }

  // testing: loop around every vertex
  for (unsigned int i=0; i < vertexPositions_.size(); i++)
  {
  cout << "Looping around vertex " << i << ":"; 
  int flag = 0;
  HalfEdge * start = &edgeAtVertex(i);
  HalfEdge * loop = start;
  do 
  {
  cout << loop->position() << " ";

  if (flag != 0) // boundary edge, exit the loop
  {
  cout << " B";
  break;
  }
  loop = &loopVertex(*loop,flag);
  }
  while (loop->position() != start->position());

  cout << endl;
  }
   */

  std::cout << "Half-edge datastructure constructed successfully." << std::endl;

  std::cout << "Statistics: " << std::endl;
  std::cout << "  Half-edges: " << halfEdges_.size() << std::endl;
  std::cout << "  Boundary half-edges: " << boundaryEdges_.size() << std::endl;
  std::cout << "  Connected components: " << connectedComponents << std::endl;

  return numOrientationFlips;
} 

// returns the previous halfedge to the given half-edge
// does so by looping around the face (pointers to previous edges are not explicitly stored), so this is slower than edgeNext
ObjFileOrientable::HalfEdge & ObjFileOrientable::edgePrevious ( HalfEdge & halfedge )
{
  HalfEdge * loop = &halfedge;
  while (edgeNext(*loop) != halfedge)
    loop = &(edgeNext(*loop));

  HalfEdge & prevEdge = *loop;

  return prevEdge;
}

// loops around the vertex (vertex is defined as the ending position of the halfedge)
// consists of taking the next edge, then taking the opposite edge
// if boundary edge encountered, can't take the opposite edge; it this case flag=1 is returned 
//     and the edge returned is the boundary edge pointing away from the vertex
// if taking the opposite edge is possible, the returned edge points into the vertex and flag is set to 0
ObjFileOrientable::HalfEdge & ObjFileOrientable::loopVertex(HalfEdge & halfedge, int & flag)
{
  HalfEdge * loop = &halfedge;
  loop = &(edgeNext(*loop));

  if (loop->isBoundary())
  {
    flag = 1;
    // return boundary edge pointing away from the vertex (there is no corresponding edge pointing into the vertex)
    HalfEdge & result = *loop;
    return result;
  }
  else
  {
    flag = 0;
    loop = &(edgeOpposite(*loop));
    // return edge pointing into the vertex
    HalfEdge & result = *loop;
    return result;
  }
}



void ObjFile::writeToFile(const string & filename, int outputMaterials)
{
  string materialFilename;
  string materialFilenameLocal;
  if (outputMaterials)
  {
    materialFilename = filename + ".mtl";
    // remove directory part from materialFilename
    char * materialFilenameTempC = (char*)materialFilename.c_str();
    char * beginString = materialFilenameTempC;
    // seek for last '/'
    for(unsigned int i=0; i< strlen(materialFilenameTempC); i++)
      if ((materialFilenameTempC[i] == '/') || (materialFilenameTempC[i] == '\\'))
        beginString = &materialFilenameTempC[i+1];

    materialFilenameLocal = string(beginString);
  }

  cout << "Writing obj to file " << filename << " ." << endl;
  if (outputMaterials)
    cout << "Outputting materials to " << materialFilename << " ." << endl;

  // open file
  ofstream fout(filename.c_str());

  if (!fout)
  {
    cout << "Error: could not write to file " << filename << endl;
    return;
  }

  // count total number of triangles
  int numTriangles = 0;
  for(unsigned int i = 0; i < groups_.size(); i++ )
    numTriangles += groups_[i].faceCount();

  fout << "# Generated automatically by the objfile class" << endl;
  fout << "# Number of vertices: " << vertexPositions_.size() << endl;
  fout << "# Number of texture coordinates: " << textureCoords_.size() << endl;
  fout << "# Number of normals: " << normals_.size() << endl;
  fout << "# Number of faces: " << numTriangles << endl;
  fout << "# Number of groups: " << groups_.size() << endl;

  if (outputMaterials)
  {
    fout << endl << "mtllib " << materialFilenameLocal << endl << endl;
  }

  // vertices...
  for (unsigned int i=0; i < vertexPositions_.size(); i++)
  {
    Vec3d pos = vertexPosition(i);

    fout << "v " << pos[0] << " " << pos[1] << " " << pos[2] << endl; 
  }

  // texture coordinates
  for (unsigned int i=0; i < textureCoords_.size(); i++)
  {
    Vec3d pos = textureCoordinate(i);

    fout << "vt " << pos[0] << " " << pos[1] << endl; 
  }

  // normals...
  for (unsigned int i=0; i < normals_.size(); i++)
  {
    Vec3d pos = normal(i);

    fout << "vn " << pos[0] << " " << pos[1] << " " << pos[2] << endl; 
  }

  // groups and faces...
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    fout << "g " << groups_[i].name() << endl;
    if (outputMaterials)
    {
      fout << "usemtl " << materials_[groups_[i].materialIndex()].name() << endl;
    }
    for( unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      Face face = groups_[i].face(iFace); // get face whose number is iFace

      fout << "f";   

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {

        Vertex vertex = face.vertex(iVertex);

        fout << " " << vertex.position()+1;

        if (vertex.hasTextureCoordinate() || vertex.hasNormal())
        {
          fout << "/";

          if (vertex.hasTextureCoordinate())
          {
            fout << vertex.textureCoordinate()+1;
          }

          if (vertex.hasNormal())
          {
            fout << "/";

            if (vertex.hasNormal())
            {
              fout << vertex.normal()+1;
            }
          }
        }


      }

      fout << endl;
    }
  }

  fout.close(); 

  if (outputMaterials)
  {
    ofstream fout(materialFilename.c_str());

    if (!fout)
    {
      cout << "Error: could not write to file " << materialFilename << endl;
      return;
    }

    for(unsigned int i=0; i< numMaterials(); i++)
    {
      fout << "newmtl " << materials_[i].name() << endl;
      fout << "illum 4" << endl;

      Vec3d Ka = materials_[i].Ka();
      Vec3d Kd = materials_[i].Kd();
      Vec3d Ks = materials_[i].Ks();
      double shininess = materials_[i].shininess() * 1000.0 / 128.0;

      fout << "Ka " << Ka[0] << " " << Ka[1] << " " << Ka[2] << endl;
      fout << "Kd " << Kd[0] << " " << Kd[1] << " " << Kd[2] << endl;
      fout << "Ks " << Ks[0] << " " << Ks[1] << " " << Ks[2] << endl;
      fout << "Ns " << shininess << endl;
      fout << endl;
    }

    fout.close();
  }
}

void ObjFile::writeToAbqFile(const string & filename)
{
  cout << "Writing obj to abq file " << filename << " ." << endl;

  if (maxFaceDegree() > 4)
  {
    cout << "Error: mesh has faces with more than 4 vertices." << endl;
    return;
  } 

  vector<double> surfaceAreas ;
  computeSurfaceAreaPerGroup(surfaceAreas);
  for(unsigned int i=0; i<surfaceAreas.size(); i++)
  {
    printf("Surface area of group %d: %G\n",i,surfaceAreas[i]);
  }

  // open file
  FILE * fout = fopen(filename.c_str(),"wa");

  if (!fout)
  {
    cout << "Error: could not write to file " << filename << endl;
    return;
  }

  // vertices...
  fprintf(fout, "*NODE\n");
  for (unsigned int i=0; i < vertexPositions_.size(); i++)
  {
    Vec3d pos = vertexPosition(i);
    fprintf(fout,"   %d,   %.15f,   %.15f,   %.15f\n",i+1,pos[0],pos[1],pos[2]);
  }

  // groups and faces...
  int faceCount=0;
  vector<int> startIndex; // for generation of element sets
  vector<int> endIndex;
  vector<std::string> groupNames;
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    printf("Num faces in group %d: %d\n",i+1,(int)(groups_[i].faceCount()));

    if (groups_[i].faceCount() == 0)
      continue;

    startIndex.push_back(faceCount+1);
    groupNames.push_back(groups_[i].name());

    // two passes: triangles and quads
    for(unsigned int numFaceVertices=3; numFaceVertices<=4; numFaceVertices++)
    {
      bool firstElement = true;
      for( unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
      {
        Face face = groups_[i].face(iFace); // get face whose number is iFace

        if (face.vertexCount() != numFaceVertices)
          continue;

        if (firstElement)
        {
          fprintf(fout,"*ELEMENT, TYPE=S%d\n",numFaceVertices);
          firstElement = false;
        }

        faceCount++;
        fprintf(fout,"%d   ",faceCount);

        for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
        {
          Vertex vertex = face.vertex(iVertex);
          fprintf(fout,",%d",vertex.position()+1);
        }

        fprintf(fout,"\n");
      }
    }

    endIndex.push_back(faceCount);
  }

  for(unsigned int i=0; i<startIndex.size(); i++)
  {
    fprintf(fout,"*ELSET,ELSET=%s,GENERATE\n",groupNames[i].c_str());
    fprintf(fout,"  %d,%d\n",startIndex[i],endIndex[i]);
  }

  fprintf(fout,"*ELSET,ELSET=EALL,GENERATE\n");
  fprintf(fout,"  1,%d\n",faceCount);

  fclose(fout);
}

void ObjFile::writeToSTLFile(const string & filename)
{
  cout << "Writing obj to STL file " << filename << " ." << endl;

  // open file
  ofstream fout(filename.c_str());

  if (!fout)
  {
    cout << "Error: could not write to file " << filename << endl;
    return;
  }

  // check if mesh is triangular
  if (!isTriangularMesh())
  {
    cout << "Error: input mesh is not triangular. " << endl;
    return;
  }

  // count total number of triangles
  int numTriangles = 0;
  for(unsigned int i = 0; i < groups_.size(); i++ )
    numTriangles += groups_[i].faceCount();

  fout << "# Generated automatically by the objfile class" << endl;
  fout << "# Number of vertices: " << vertexPositions_.size() << endl;
  fout << "# Number of faces: " << numTriangles << endl;


  fout << "solid" << endl;

  // groups and faces...
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for( unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      fout << "  facet normal ";

      // get the face data
      Vertex v0 = face.vertex(0);
      Vertex v1 = face.vertex(1);
      Vertex v2 = face.vertex(2);

      Vec3d p0 = vertexPosition(v0.position());
      Vec3d p1 = vertexPosition(v1.position());
      Vec3d p2 = vertexPosition(v2.position());

      // compute the face normal
      Vec3d normal = norm(cross(p1-p0,p2-p0));

      fout << normal[0] << " " << normal[1] << " " << normal[2] << endl;

      fout << "    outer loop" << endl;

      fout << "      vertex " << p0[0] << " " << p0[1] << " " << p0[2] << endl;
      fout << "      vertex " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
      fout << "      vertex " << p2[0] << " " << p2[1] << " " << p2[2] << endl;

      fout << "    endloop" << endl;
      fout << "  endfacet" << endl;

    }
  }

  fout << "endsolid" << endl;

  fout.close();

}



double triangleSurfaceArea(Vec3d p0, Vec3d p1, Vec3d p2)
{
  return 0.5 * (len(cross(p1-p0,p2-p0)));
}

void ObjFile::getCentroids(std::vector<Vec3d> & centroids)
{
  interpolateToCentroids(vertexPositions_,centroids);
}

void ObjFile::interpolateToCentroids(std::vector<double> & nodalData, std::vector<double> & centroidData)
{
  int faceIndex = 0;
  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      double data = 0;
      for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {
        unsigned int index = face.vertex(iVertex).position();
        data += nodalData[index];
      }

      data /= face.vertexCount();
      centroidData[faceIndex] = data;

      faceIndex++;
    }
  }
}

void ObjFile::interpolateToCentroids(std::vector<Vec3d> & nodalData, std::vector<Vec3d> & centroidData)
{
  int faceIndex = 0;
  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      Vec3d data(0,0,0);
      for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {
        unsigned int index = face.vertex(iVertex).position();
        data += nodalData[index];
      }

      data /= face.vertexCount();
      centroidData[faceIndex] = data;

      faceIndex++;
    }
  }
}

void ObjFile::computeSurfaceAreaPerVertex()
{
  for (unsigned int i=0; i < numVertices(); i++)
    surfaceAreaPerVertex_.push_back(0);


  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      double faceSurfaceArea_ = faceSurfaceArea(face);

      for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
      {
        unsigned int index = face.vertex(iVertex).position();
        surfaceAreaPerVertex_[index] += faceSurfaceArea_ / face.vertexCount(); // each vertex owns an equal share of the face
      }
    }
  }
}


void ObjFile::uniformScale(Vec3d center, double factor)
{
  for (unsigned int i=0; i < vertexPositions_.size(); i++) // over all vertices
    vertexPositions_[i] = center + factor * (vertexPositions_[i] - center);
}



double ObjFile::computeVolume()
{

  double volume = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() != 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") which is not a triangle." << endl;

      // base vertex
      Vec3d v0 = vertexPositions_[face.vertex(0).position()];
      Vec3d v1 = vertexPositions_[face.vertex(1).position()];
      Vec3d v2 = vertexPositions_[face.vertex(2).position()];

      Vec3d normal = cross(v1-v0,v2-v0);
      Vec3d center = 1.0 / 3 * (v0 + v1 + v2);

      volume += dot(normal,center);

    }
  }

  volume /= 6.0;

  return volume;
}

Vec3d ObjFile::computeCenterOfMass_Vertices()
{
  Vec3d center(0,0,0);

  for (unsigned int i=0; i < vertexPositions_.size(); i++) // over all vertices
    center += vertexPositions_[i];

  center /= vertexPositions_.size();

  return center;
}

Vec3d ObjFile::computeCenterOfMass_Triangles(vector<double> & groupDensities)
{
  Vec3d centerOfMass = vl_zero;

  double totalMass=0.0;
  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    double density = groupDensities[i];
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      double area = faceSurfaceArea(face);
      double mass = density*area;
      totalMass += mass;
      Vec3d centroid = faceCentroid(face);       
      centerOfMass += mass * centroid;
    }
  }

  centerOfMass /= totalMass;

  return centerOfMass;
}

void ObjFile::computeInertiaTensor_Triangles(double mass, double IT[6])
{
  double surface = computeSurfaceArea();
  double surfaceMassDensity = mass / surface;
  vector<double> groupDensities;
  for(unsigned int i=0; i<groups_.size(); i++)
    groupDensities.push_back(surfaceMassDensity);
  computeInertiaTensor_Triangles(groupDensities, IT);
}

double ObjFile::computeMass_Triangles(vector<double> & groupDensities)
{
  double totalMass=0.0;
  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    double density = groupDensities[i];
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
      {
        printf("Warning: encountered a face with fewer than three vertices.\n");
        continue;
      }

      double area = faceSurfaceArea(face);
      double mass = density*area;
      totalMass += mass;
    }
  }

  return totalMass;
}

void ObjFile::computeInertiaTensor_Triangles(vector<double> & groupDensities, double IT[6])
{
  Vec3d centerOfMass = vl_zero;
  memset(IT,0,sizeof(double)*6);

  double totalMass=0.0;
  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    double density = groupDensities[i];
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
      {
        printf("Warning: encountered a face with fewer than three vertices.\n");
        continue;
      }

      double area = faceSurfaceArea(face);
      double mass = density*area;
      totalMass += mass;
      Vec3d centroid = faceCentroid(face);       
      centerOfMass += mass * centroid;

      unsigned int index0 = face.vertex(0).position();
      Vec3d v0 = vertexPositions_[index0];
      for (unsigned int iVertex = 1; iVertex < face.vertexCount()-1; ++iVertex )
      {
        unsigned int index1 = face.vertex(iVertex).position();
        Vec3d v1 = vertexPositions_[index1];
        unsigned int index2 = face.vertex(iVertex+1).position();
        Vec3d v2 = vertexPositions_[index2];
        double ITTriangle[6];
        computeSpecificInertiaTensor(v0, v1, v2, ITTriangle);

        double triangleArea = 0.5 * len(cross(v1-v0,v2-v0));
        for(int j=0; j<6; j++)
          IT[j] += triangleArea * density * ITTriangle[j];
      }
    }
  }

  centerOfMass /= totalMass;

  // IT is now the center around the origin  
  // transfer tensor to the center of mass
  double a = centerOfMass[0];
  double b = centerOfMass[1];
  double c = centerOfMass[2];

  double correction[6] = 
  { b*b + c*c, -a*b, -a*c,
    a*a + c*c, -b*c,
    a*a + b*b };

  for(int i=0; i<6; i++)
    IT[i] -= totalMass * correction[i];

}

Vec3d ObjFile::computeCenterOfMass_Triangles()
{
  Vec3d centerOfMass = vl_zero;

  double totalArea=0.0;
  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      double area = faceSurfaceArea(face);
      totalArea += area;
      Vec3d centroid = faceCentroid(face);       
      centerOfMass += area * centroid;
    }
  }

  centerOfMass /= totalArea;

  return centerOfMass;
}

void ObjFile::faceSurfaceAreas(vector<double> & surfaceAreas)
{
  int faceIndex = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      surfaceAreas[faceIndex] = faceSurfaceArea(face);

      faceIndex++;
    }
  }

}

double ObjFile::computeSurfaceArea()
{
  double area = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      /*
         if (face.vertexCount() != 3)
         cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") which is not a triangle." << endl;
       */

      area += faceSurfaceArea(face);

    }
  }

  return area;
}

void ObjFile::computeSurfaceAreaPerGroup(vector<double> & surfaceAreas)
{
  surfaceAreas.clear();

  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    double area = 0;
    // over all faces
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      /*
         if (face.vertexCount() != 3)
         cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") which is not a triangle." << endl;
       */

      area += faceSurfaceArea(face);

    }

    surfaceAreas.push_back(area);
  }

}

void ObjFile::computeMassPerVertex(vector<double> & groupSurfaceMassDensity, vector<double> & masses)
{
  masses.clear();
  for(unsigned int i=0; i<numVertices(); i++)
    masses.push_back(0.0);

  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    // over all faces
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      /*
         if (face.vertexCount() != 3)
         cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") which is not a triangle." << endl;
       */

      double faceSurfaceArea_ = faceSurfaceArea(face);
      for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
        masses[face.vertex(iVertex).position()] += groupSurfaceMassDensity[i] * faceSurfaceArea_ / face.vertexCount();
    }
  }
}

void ObjFileOrientable::determineIfSurfaceHasBoundary()
{
  for (unsigned int i=0; i < halfEdges_.size(); i++)
  {
    if (halfEdges_[i].isBoundary())
    {
      hasBoundary_ = true;
      return;
    }
  }

  hasBoundary_ = false;
}

void ObjFileOrientable::CopyHalfEdgeTopologyFrom(ObjFileOrientable * source) // makes the half-edge topological info equal to that of source
{

  halfEdges_ = source->halfEdges_;
  boundaryEdges_ = source->boundaryEdges_;
  connectedComponents = source->connectedComponents;
  edgesAtVertices_ = source->edgesAtVertices_;
  edgesAtFaces_ = source->edgesAtFaces_;
  hasBoundary_ = source->hasBoundary_;

}

// warning: this only cares about the first triangle in a face (assumes planar face)
Vec3d ObjFile::faceNormal(Face & face)
{
  // the three vertices
  unsigned int index0 = face.vertex(0).position();
  unsigned int index1 = face.vertex(1).position();
  unsigned int index2 = face.vertex(2).position();

  Vec3d pos0 = vertexPositions_[index0];
  Vec3d pos1 = vertexPositions_[index1];
  Vec3d pos2 = vertexPositions_[index2];

  Vec3d normal = norm(cross(pos1-pos0,pos2-pos0));
  return normal;
}

void ObjFile::setNormalsToFaceNormals()
{
  // over all faces
  normals_.clear();
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face * faceHandle = groups_[i].faceHandle(iFace); // get face whose number is iFace

      if (faceHandle->vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      Vec3d normal = faceNormal(*faceHandle);
      addVertexNormal(normal);

      // over all vertices of the face
      for (unsigned k=0; k<faceHandle->vertexCount(); k++)
      {
        Vertex * vertex = faceHandle->vertexHandle(k);
        vertex->setNormal(numNormals()-1);
      }
    }
  }
}


void ObjFile::setNormalsToAverageFaceNormals()
{
  vector<Vec3d> normalBuffer(numVertices(),Vec3d(0,0,0));
  vector<unsigned int> normalCount(numVertices(),0);

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      // the three vertices
      unsigned int index0 = face.vertex(0).position();
      unsigned int index1 = face.vertex(1).position();
      unsigned int index2 = face.vertex(2).position();

      Vec3d pos0 = vertexPositions_[index0];
      Vec3d pos1 = vertexPositions_[index1];
      Vec3d pos2 = vertexPositions_[index2];

      Vec3d normal = norm(cross(pos1-pos0,pos2-pos0));
      // this works even for non-triangle meshes

      normalBuffer[index0] += normal;
      normalBuffer[index1] += normal;
      normalBuffer[index2] += normal;

      normalCount[index0]++;
      normalCount[index1]++;
      normalCount[index2]++;

    }
  }

  bool errorMessageSeen=false;
  // normalize the normals
  for (unsigned int i=0; i < numVertices(); i++)
  {
    if (normalCount[i] == 0)
    {
      if (!errorMessageSeen)
        cout << "Warning: encountered a vertex not belonging to any triangle (suppressing further warnings)" << endl;
      errorMessageSeen = true;
      normalBuffer[i] = Vec3d(1,0,0); // assign some bogus normal
    }
    else
      normalBuffer[i] = norm(normalBuffer[i]);
  }

  // register new normals with the objfile data structure
  normals_.clear();
  for (unsigned int i=0; i < numVertices(); i++)
  {
    addVertexNormal(normalBuffer[i]);     
  }  

  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    Group * group = &(groups_[i]);
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      Face * face = group->faceHandle(iFace);

      if (face->vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for (unsigned k=0; k<face->vertexCount(); k++)
      {
        Vertex * vertex = face->vertexHandle(k);
        vertex->setNormal(vertex->position());
      }
    }
  }

}


unsigned int ObjFile::nearestVertex(Vec3d queryPos, double & dist)
{
  double closestDist2 = DBL_MAX;
  double candidateDist2;
  unsigned int indexClosest = 0;
  for(unsigned int i=0; i< numVertices(); i++)
  {
    Vec3d relPos = vertexPosition(i)-queryPos;
    if ((candidateDist2 = dot(relPos,relPos)) < closestDist2)
    {
      closestDist2 = candidateDist2;
      indexClosest = i;
    } 
  }

  dist = sqrt(closestDist2);
  return indexClosest;
}

double ObjFile::minEdgeLength()
{
  double minLength = -1;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for (unsigned k=0; k<face.vertexCount(); k++)
      {
        Vec3d pos0 = vertexPositions_[face.vertex(k).position()];
        Vec3d pos1 = vertexPositions_[face.vertex((k+1) % face.vertexCount()).position()];
        double length = len(pos1-pos0);

        if (minLength < 0) // only the first time
          minLength = length;
        else if (length < minLength)
          minLength = length;


      }

    }
  }

  return minLength;

}

double ObjFile::medianEdgeLength()
{
  vector<double> lengths;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for (unsigned k=0; k<face.vertexCount(); k++)
      {
        Vec3d pos0 = vertexPositions_[face.vertex(k).position()];
        Vec3d pos1 = vertexPositions_[face.vertex((k+1) % face.vertexCount()).position()];
        double length = len(pos1-pos0);

        lengths.push_back(length);
      }

    }
  }

  sort(lengths.begin(),lengths.end());

  return lengths[lengths.size() / 2];

}

double ObjFile::maxEdgeLength()
{
  double maxLength = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for (unsigned k=0; k<face.vertexCount(); k++)
      {
        Vec3d pos0 = vertexPositions_[face.vertex(k).position()];
        Vec3d pos1 = vertexPositions_[face.vertex((k+1) % face.vertexCount()).position()];
        double length = len(pos1-pos0);

        if (length > maxLength)
          maxLength = length;
      }
    }
  }

  return maxLength;
}

double ObjFile::minEdgeLength(int * vtxa, int * vtxb)
{
  *vtxa = *vtxb = -1;

  double minLength = -1;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for (unsigned k=0; k<face.vertexCount(); k++)
      {
        Vec3d pos0 = vertexPositions_[face.vertex(k).position()];
        Vec3d pos1 = vertexPositions_[face.vertex((k+1) % face.vertexCount()).position()];
        double length = len(pos1-pos0);

        if (minLength < 0) // only the first time
          minLength = length;
        else if (length < minLength)
          minLength = length;

        *vtxa = face.vertex(k).position();
        *vtxb = face.vertex((k+1) % face.vertexCount()).position();
      }
    }
  }

  return minLength;
}

double ObjFile::maxEdgeLength(int * vtxa, int * vtxb)
{
  *vtxa = *vtxb = -1;

  double maxLength = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      if (face.vertexCount() < 3)
        cout << "Warning: encountered a face (group=" << i << ",face=" << iFace << ") with fewer than 3 vertices." << endl;

      for (unsigned k=0; k<face.vertexCount(); k++)
      {
        Vec3d pos0 = vertexPositions_[face.vertex(k).position()];
        Vec3d pos1 = vertexPositions_[face.vertex((k+1) % face.vertexCount()).position()];
        double length = len(pos1-pos0);

        if (length > maxLength)
          maxLength = length;

        *vtxa = face.vertex(k).position();
        *vtxb = face.vertex((k+1) % face.vertexCount()).position();
      }
    }
  }

  return maxLength;
}

unsigned int ObjFile::numFaces()
{
  unsigned int counter = 0;
  for (unsigned int i=0; i < groups_.size(); i++)
    counter += groups_[i].faceCount();    

  return counter;
}  

void ObjFile::setNormalsToPseudoNormals()
{
  // nuke any previous normals
  normals_.clear();

  // registers pseudonormals as the new normals
  for(unsigned int i=0; i < numVertices(); i++)
  {
    normals_.push_back(pseudoNormals_[i]);
  }

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face * face = groups_[i].faceHandle(iFace); // get face whose number is iFace

      // over all vertices of the face
      for (unsigned k=0; k<face->vertexCount(); k++)
      {
        Vertex * vertex = face->vertexHandle(k);
        vertex->setNormal(vertex->position());
      }
    }
  }

}

void ObjFile::computePseudoNormals()
{
  vector<int> vertexDegree(numVertices());
  for(unsigned int i=0; i<numVertices(); i++)
    pseudoNormals_.push_back(vl_zero);

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace

      // over all vertices
      for (unsigned k=0; k<face.vertexCount(); k++)
      {
        // compute angle at that vertex in radians
        Vec3d pos = vertexPositions_[face.vertex(k).position()];
        Vec3d posNext, posPrev;
        if (k != face.vertexCount() - 1)
          posNext = vertexPositions_[face.vertex(k+1).position()];
        else
          posNext = vertexPositions_[face.vertex(0).position()];

        if (k != 0)
          posPrev = vertexPositions_[face.vertex(k-1).position()];
        else
          posPrev = vertexPositions_[face.vertex(face.vertexCount() - 1).position()];

        double lenNext = len(posNext-pos);
        double lenPrev = len(posPrev-pos);

        double angle = acos(dot(posNext-pos,posPrev-pos)/lenNext/lenPrev);
        Vec3d normal = norm(cross(posNext-pos, posPrev-pos));

        if (my_isnan(normal[0]) || my_isnan(normal[1]) || my_isnan(normal[2]))
        {
          cout << "Error (when computing vertex pseudonormals): nan encountered (face with zero surface area)." << endl;
          cout << "Group: " << i << " Face: " << iFace << " " << endl;
          exit(1);
          //cout << "  vtx0: " << index0 << " vtx1: " << index1 << " vtx2: " << index2 << endl;
          //cout << "  "  << p0 << endl;
          //cout << "  "  << p1 << endl;
          //cout << "  "  << p2 << endl;
          //cout << "Feature: " << normali << endl;
        }
        else
        {
          if ((lenNext == 0) || (lenPrev == 0) || my_isnan(angle))
          {
            cout << "Warning (when computing vertex pseudonormals): encountered zero-length edge" << endl;
            cout << "  lenNext: " << lenNext << " lenPrev: " << lenPrev << " angle: " << angle << endl;
          }
          else  
          {
            //if(face.vertex(k).position() == 159553)
            //{
            //printf("angle: %G normal: %G %G %G\n", angle, normal[0], normal[1], normal[2]);
            //}
            pseudoNormals_[face.vertex(k).position()] += angle * normal;
            vertexDegree[face.vertex(k).position()]++;
          }
        }
      }
    }
  }

  for(unsigned int i=0; i<numVertices(); i++)
  {
    if (vertexDegree[i] != 0)
    {
      Vec3d pseudoNormalRaw = pseudoNormals_[i];
      pseudoNormals_[i] = norm(pseudoNormalRaw);
      if (my_isnan(pseudoNormals_[i][0]) || my_isnan(pseudoNormals_[i][1]) || my_isnan(pseudoNormals_[i][2]))
      {
        cout << "Error (when computing vertex pseudonormals): nan encountered." << endl;
        cout << "Vertex: " << i << " pseudoNormal=" << pseudoNormals_[i][0] << " " << pseudoNormals_[i][1] << " " << pseudoNormals_[i][2] << endl;
        cout << "  Pseudonormal before normalization=";
        printf("%G %G %G\n", pseudoNormalRaw[0], pseudoNormalRaw[1], pseudoNormalRaw[2]); 
        cout << "  Vertex degree=" << vertexDegree[i] << endl;
        exit(1);
      }
    }
    else
      pseudoNormals_[i] = vl_zero;
  }
}

void ObjFile::computeEdgePseudoNormals()
{
  edgePseudoNormals_.clear();

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace
      Vec3d pos0 = vertexPositions_[face.vertex(0).position()];
      Vec3d pos1 = vertexPositions_[face.vertex(1).position()];
      Vec3d pos2 = vertexPositions_[face.vertex(2).position()];
      Vec3d normal = norm(cross(pos1-pos0, pos2-pos0));

      if (my_isnan(normal[0]) || my_isnan(normal[1]) || my_isnan(normal[2]))
      {
        cout << "Error: nan encountered (face with zero surface area)." << endl;
        cout << "Group: " << i << " Face: " << iFace << " " << endl;
        exit(1);
      }

      // over all edges at the face
      for (unsigned k=0; k<face.vertexCount(); k++)
      {

        unsigned int startVertex = face.vertex(k).position();
        unsigned int endVertex = face.vertex( (k+1) % face.vertexCount()).position();

        pair<unsigned int, unsigned int> edge;
        if (startVertex < endVertex)
        {
          edge.first = startVertex;
          edge.second = endVertex;
        }
        else
        {
          edge.first = endVertex;
          edge.second = startVertex;
        }

        map< pair<unsigned int, unsigned int>, Vec3d > :: iterator iter = edgePseudoNormals_.find(edge);

        if (iter == edgePseudoNormals_.end())
        {
          edgePseudoNormals_.insert(make_pair(edge,normal));
        }
        else
        {
          iter->second += normal;
        }
      }
    }
  }

  // normalize normals
  map< pair<unsigned int, unsigned int>, Vec3d > :: iterator iter;
  for(iter = edgePseudoNormals_.begin(); iter != edgePseudoNormals_.end(); ++iter)
  {
    Vec3d normal = norm(iter->second);
    if (my_isnan(normal[0]) || my_isnan(normal[1]) || my_isnan(normal[2]))
    {
      cout << "Warning (while computing edge pseudonormals): nan encountered (face with zero surface area)." << endl;
      normal[0] = 1; normal[1] = 0; normal[2] = 0;
    }
    iter->second = normal;
  }
}

int ObjFile::removeZeroAreaFaces()
{
  int numZeroAreaFaces = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      ObjFile::Face face = groups_[i].face(iFace); // get face whose number is iFace
      Vec3d pos0 = vertexPositions_[face.vertex(0).position()];
      Vec3d pos1 = vertexPositions_[face.vertex(1).position()];
      Vec3d pos2 = vertexPositions_[face.vertex(2).position()];
      Vec3d normal = norm(cross(pos1-pos0, pos2-pos0));

      bool identicalVertex = false;
      for(unsigned int jj=0; jj< face.vertexCount(); jj++)
        for(unsigned int kk=jj+1; kk< face.vertexCount(); kk++)
        {
          if (face.vertex(jj).position() == face.vertex(kk).position())
            identicalVertex = true;
        }

      if (my_isnan(normal[0]) || my_isnan(normal[1]) || my_isnan(normal[2]) || identicalVertex)
      {
        groups_[i].removeFace(iFace);
        iFace--;
        numZeroAreaFaces++;
      }
    }
  }

  return numZeroAreaFaces;
}

// returns 1 on success, 0 otherwise
int ObjFile::edgePseudoNormal(unsigned int i, unsigned int j, Vec3d * pseudoNormal) 
{ 
  pair<unsigned int,unsigned int> edge;

  if (i < j)
  {
    edge.first = i;
    edge.second = j;
  }
  else
  {
    edge.first = j;
    edge.second = i;
  }

  map< pair<unsigned int, unsigned int> , Vec3d > :: iterator iter = edgePseudoNormals_.find(edge);

  if (iter != edgePseudoNormals_.end())
  {
    *pseudoNormal = iter->second;
    return 0;
  }

  return 1;
}

double ObjFile::faceSurfaceArea(Face & face)
{ 
  double faceSurfaceArea_ = 0;

  // base vertex
  unsigned int indexBase = face.vertex(0).position();
  Vec3d basePos = vertexPositions_[indexBase];

  for ( unsigned int iVertex = 1; iVertex < face.vertexCount()-1; ++iVertex )
  {
    unsigned int index1 = face.vertex(iVertex).position();
    unsigned int index2 = face.vertex(iVertex+1).position();

    Vec3d pos1 = vertexPositions_[index1];
    Vec3d pos2 = vertexPositions_[index2];

    faceSurfaceArea_ += fabs(triangleSurfaceArea(basePos,pos1,pos2));
  }

  return faceSurfaceArea_;
}

Vec3d ObjFile::faceCentroid(Face & face)
{ 
  Vec3d centroid = vl_zero;

  for ( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
    centroid += vertexPositions_[face.vertex(iVertex).position()];

  centroid *= 1.0 / face.vertexCount();

  return centroid;
}

void ObjFile::computeSpecificInertiaTensor(Vec3d & v0, Vec3d & v1, Vec3d & v2, double t[6])
{

  t[0] = (v0[1]*v0[1] + v0[2]*v0[2] + v1[1]*v1[1] + v1[2]*v1[2] + 
      v1[1]*v2[1] + v2[1]*v2[1] + v0[1]*(v1[1] + v2[1]) + 
      v1[2]*v2[2] + v2[2]*v2[2] + v0[2]*(v1[2] + v2[2]))/6; 

  t[1] = (-2*v1[0]*v1[1] - v1[1]*v2[0] - v0[1]*(v1[0] + v2[0])
      - v1[0]*v2[1] - 2*v2[0]*v2[1] - 
      v0[0]*(2*v0[1] + v1[1] + v2[1]))/12; 

  t[2] = (-2*v1[0]*v1[2] - v1[2]*v2[0] - v0[2]*(v1[0] + v2[0]) - 
      v1[0]*v2[2] - 2*v2[0]*v2[2] - 
      v0[0]*(2*v0[2] + v1[2] + v2[2]))/12; 

  t[3] =  (v0[0]*v0[0] + v0[2]*v0[2] + v1[0]*v1[0] + v1[2]*v1[2] + 
      v1[0]*v2[0] + v2[0]*v2[0] + v0[0]*(v1[0] + v2[0]) + 
      v1[2]*v2[2] + v2[2]*v2[2] + v0[2]*(v1[2] + v2[2]))/6; 

  t[4] = (-2*v1[1]*v1[2] - v1[2]*v2[1] - 
      v0[2]*(v1[1] + v2[1]) - v1[1]*v2[2] - 2*v2[1]*v2[2] - 
      v0[1]*(2*v0[2] + v1[2] + v2[2]))/12; 

  t[5] = (v0[0]*v0[0] + v0[1]*v0[1] + v1[0]*v1[0] + v1[1]*v1[1] + 
      v1[0]*v2[0] + v2[0]*v2[0] + v0[0]*(v1[0] + v2[0]) + 
      v1[1]*v2[1] + v2[1]*v2[1] + v0[1]*(v1[1] + v2[1]))/6;
}

void ObjFile::dirname(char * path, char * result)
{
  // seek for last '/' or '\'
  char * ch = path;
  int lastPos = -1;
  int pos=0;

  while (*ch != 0)
  {
    if (*ch == '\\')
      lastPos = pos;

    if (*ch == '/')
      lastPos = pos;

    ch++;
    pos++;
  }

  if (lastPos != -1)
  {
    memcpy(result,path,sizeof(char)*lastPos);
    result[lastPos] = 0;
  }
  else
  {
    result[0] = '.';
    result[1] = 0;
  }

}

void ObjFile::parseMaterials(char * objFilename, char * materialFilename)
{
  FILE* file;
  char  buf[128];
  unsigned int nummaterials;

  char objFilenameCopy[4096];
  strcpy(objFilenameCopy,objFilename);

  char dir[4096];
  dirname(objFilenameCopy,dir);
  char filename[4096];
  strcpy(filename, dir);
  strcat(filename, "/");
  strcat(filename, materialFilename);

  file = fopen(filename, "r");
  if (!file) 
  {
    fprintf(stderr, "glmReadMTL() failed: can't open material file %s.\n", filename);
    std::string message = "Failed to open material file '";
    message.append( filename );
    message.append( "'" );
    throw ObjFileException( message );
  }

  /* count the number of materials in the file */
  nummaterials = 1;
  while(fscanf(file, "%s", buf) != EOF) {
    switch(buf[0]) {
      case '#':               /* comment */
        /* eat up rest of line */
        fgets(buf, sizeof(buf), file);
        break;
      case 'n':               /* newmtl */
        fgets(buf, sizeof(buf), file);
        nummaterials++;
        sscanf(buf, "%s %s", buf, buf);
        break;
      default:
        /* eat up rest of line */
        fgets(buf, sizeof(buf), file);
        break;
    }
  }

  rewind(file);

  double Ka[3];
  double Kd[3];
  double Ks[3];
  double shininess;
  string matName;

  /* now, read in the data */
  nummaterials = 0;
  while(fscanf(file, "%s", buf) != EOF) 
  {
    switch(buf[0]) {
      case '#':               /* comment */
        /* eat up rest of line */
        fgets(buf, sizeof(buf), file);
        break;
      case 'n':               /* newmtl */

        if (nummaterials >= 1) // flush previous material
          addMaterial(matName,Vec3d(Ka[0],Ka[1],Ka[2]),Vec3d(Kd[0],Kd[1],Kd[2]),Vec3d(Ks[0],Ks[1],Ks[2]),shininess);

        // reset to default
        Ka[0] = 0.1; Ka[1] = 0.1; Ka[2] = 0.1;
        Kd[0] = 0.5; Kd[1] = 0.5; Kd[2] = 0.5;
        Ks[0] = 0.0; Ks[1] = 0.0; Ks[2] = 0.0;
        shininess = 65;

        fgets(buf, sizeof(buf), file);
        sscanf(buf, "%s %s", buf, buf);
        nummaterials++;
        matName = string(buf);
        break;

      case 'N':
        if (buf[1] == 's')
        {
          fscanf(file, "%lf", &shininess);
          // wavefront shininess is from [0, 1000], so scale for OpenGL 
          shininess *= 128.0 / 1000.0;
        }
        else
          fgets(buf, sizeof(buf), file); //eat rest of line

        break;

      case 'K':
        switch(buf[1]) {
          case 'd':
            fscanf(file, "%lf %lf %lf", &Kd[0], &Kd[1], &Kd[2]);
            break;
          case 's':
            fscanf(file, "%lf %lf %lf",  &Ks[0], &Ks[1], &Ks[2]);
            break;
          case 'a':
            fscanf(file, "%lf %lf %lf",  &Ka[0], &Ka[1], &Ka[2]);
            break;
          default:
            /* eat up rest of line */
            fgets(buf, sizeof(buf), file);
            break;
        }
        break;

        /*
           case 'm': // texture file, allow only one texture file, has to be ppm format

        // different materials can have different textures
        // however, each single material can have only one texture

        // if texture file not already assigned it, assign it, 
        // and then load the texture file into main memory and assign texture name to it

        fgets(buf, sizeof(buf), file);
        sscanf(buf, "%s %s", buf, buf);
        model->materials[nummaterials].textureFile = strdup(buf);
         */


        break;

      default:
        /* eat up rest of line */
        fgets(buf, sizeof(buf), file);
        break;
    }
  }

  if (nummaterials >= 1) // flush last material
    addMaterial(matName,Vec3d(Ka[0],Ka[1],Ka[2]),Vec3d(Kd[0],Kd[1],Kd[2]),Vec3d(Ks[0],Ks[1],Ks[2]),shininess);

}

void ObjFile::initSurfaceSampling()
{
  if(!isTriangularMesh())
  {
    printf("Error in init surface sampling: surface not triangular.\n");
    throw ObjFileException("Error in init surface sampling: surface not triangular.");
  }

  double totalSurfaceArea = computeSurfaceArea();
  double area = 0;

  // over all faces
  for(unsigned int i = 0; i < groups_.size(); i++ )
  {
    for(unsigned int iFace = 0; iFace < groups_[i].faceCount(); ++iFace )
    {
      surfaceSamplingAreas.push_back(
          make_pair(area,groups_[i].faceHandle(iFace)));
      ObjFile::Face face = groups_[i].face(iFace); 
      area += faceSurfaceArea(face) / totalSurfaceArea;
    }
  }

}

Vec3d ObjFile::sampleSurfacePosition(double sample)
{
  unsigned int facePosition;
  for(facePosition=0; facePosition< surfaceSamplingAreas.size()-1; facePosition++)
  {
    if ((surfaceSamplingAreas[facePosition].first <= sample) && 
        (surfaceSamplingAreas[facePosition+1].first > sample))
      break;
  }

  // facePosition now contains the index of the face to sample from
  Face * face = surfaceSamplingAreas[facePosition].second;

  // sample at random on the face
  double alpha, beta;
  do
  {
    alpha = 1.0 * rand() / RAND_MAX;
    beta = 1.0 * rand() / RAND_MAX;
  }  
  while (alpha + beta > 1);

  double gamma = 1 - alpha - beta;

  Vec3d v0 = vertexPositions_[face->vertex(0).position()];
  Vec3d v1 = vertexPositions_[face->vertex(1).position()];
  Vec3d v2 = vertexPositions_[face->vertex(2).position()];

  Vec3d sampledPos = alpha * v0 + beta * v1 + gamma * v2;
  return sampledPos;
}


void ObjFile::Group::removeFace(unsigned int i)
{
  faces_.erase(faces_.begin() + i);
}

ObjFileException::ObjFileException( const std::string& reason )
{
  std::ostringstream oss;
  oss << "Error:  " << reason;
  std::cout << std::endl << oss.str() << std::endl;
  reason_ = oss.str();
}

ObjFileException::ObjFileException( const std::string& reason, const std::string& filename,
    unsigned int line)
{
  std::ostringstream oss;
  oss << "Error in file '" << filename
    << "', line " << line << ": "
    << reason;
  std::cout << std::endl << oss.str() << std::endl;
  reason_ = oss.str();
}

bool ObjFile::Material::operator==(Material & mat2)
{
  if (len(Ka_ - mat2.Ka_) > 1e-7)
    return false;

  if (len(Kd_ - mat2.Kd_) > 1e-7)
    return false;

  if (len(Ks_ - mat2.Ks_) > 1e-7)
    return false;

  if (fabs(shininess_ - mat2.shininess_) > 1e-7)
    return false;

  return true;
}

int ObjFile::removeDuplicatedMaterials()
{
  unsigned int numMaterials_ = numMaterials();
  vector<int> reNumberVector(numMaterials_);  

  vector<Material> newMaterials;

  // detected duplicated materials
  for(unsigned int i=0; i<numMaterials_; i++)
  {
    bool newMaterial = true; 
    for(unsigned int j=0; j<newMaterials.size(); j++)
    {
      if (newMaterials[j] == materials_[i])
      {
        newMaterial = false;
        reNumberVector[i] = j;
        break;
      }
    }

    if (newMaterial)
    {
      newMaterials.push_back(materials_[i]);
      reNumberVector[i] = newMaterials.size() - 1;
    }
  } 

  materials_ = newMaterials;

  // correct the groups
  for(unsigned int i=0; i<numGroups(); i++)
    groups_[i].setMaterialIndex(reNumberVector[groups_[i].materialIndex()]);

  return materials_.size();
}

void ObjFile::dumpGeometry(int * numVertices, double ** vertices, int * numTriangles
    , int ** triangles)
{
  // set vertices
  *numVertices = vertexPositions_.size();
  *vertices = (double*) malloc (sizeof(double) * 3 * *numVertices);
  for(int i=0; i< *numVertices; i++)
  {
    Vec3d vtx = vertexPosition(i);
    (*vertices)[3*i+0] = vtx[0];
    (*vertices)[3*i+1] = vtx[1];
    (*vertices)[3*i+2] = vtx[2];
  }

  if (numTriangles == NULL)
  {
    printf("Dumped %d vertices.\n", *numVertices);
    return;
  }

  // set triangles
  *numTriangles = 0;
  for(unsigned int i=0; i < groups_.size(); i++) // over all groups
    for (unsigned int j=0; j < groups_[i].faceCount(); j++) // over all faces
    {
      Face * face = groups_[i].faceHandle(j);
      if (face->vertexCount() < 3)
        continue;
      *numTriangles += face->vertexCount() - 2;
    }

  *triangles = (int*) malloc (sizeof(int) * 3 * *numTriangles);

  int tri = 0;
  for(unsigned int i=0; i < groups_.size(); i++) // over all groups
    for (unsigned int j=0; j < groups_[i].faceCount(); j++) // over all faces
    {
      Face * face = groups_[i].faceHandle(j);
      if (face->vertexCount() < 3)
      {
        printf("Warning: encountered a face with fewer than 3 vertices.\n");
        continue;
      }

      unsigned int faceDegree = face->vertexCount();

      // triangulate the face

      // get the vertices:
      vector<Vertex> vertices;
      for(unsigned int k=0; k<face->vertexCount(); k++)
        vertices.push_back(face->vertex(k));

      (*triangles)[3*tri+0] = vertices[0].position();
      (*triangles)[3*tri+1] = vertices[1].position();
      (*triangles)[3*tri+2] = vertices[2].position();
      tri++;

      for(unsigned int k=2; k<faceDegree-1; k++)
      {
        (*triangles)[3*tri+0] = vertices[0].position();
        (*triangles)[3*tri+1] = vertices[k].position();
        (*triangles)[3*tri+2] = vertices[k+1].position();
        tri++;
      }
    }

  printf("Dumped %d vertices and %d triangles.\n", *numVertices, *numTriangles);
}
