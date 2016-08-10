#include <fstream> 
#include <io/FBemReader.h> 

//##############################################################################
// This function is directly copied from
//  modec/src/multipole/fast_fit_multipole.cpp:fbemInToGeometry
//##############################################################################
bool FBemReader::
_ReadFBemInputToGeometry(const char *fBemInputFile, std::vector<Point3d> &verts, std::vector<Tuple3ui> &tris)
{
    // open file
    std::ifstream fileObj;
    fileObj.open(fBemInputFile, std::ios_base::in);
    if(!fileObj.good())
    {
        std::cerr << "FAIL: opening FBEM in file" << std::endl;
        fileObj.close();
        return false;
    }

    std::string lineString;
    int numElements;
    int numVerts;

    // read about num elements and verts
    for(int i = 0; i < 4; i++)
    {
        std::getline(fileObj, lineString);
    }

    std::stringstream sizesStream(lineString);

    sizesStream >> numElements >> numVerts;
    std::cout << "num elements: " << numElements << std::endl;
    std::cout << "num verts: " << numVerts << std::endl;

    verts.resize(numVerts);
    tris.resize(numElements);

    // strip material before "$ Nodes:"
    do
    {
        std::getline(fileObj, lineString);
    }
    while(lineString.size() < 1 || lineString[0] != '$');

    // read verts
    int vid;
    Point3d vert;

    for(int v = 0; v < numVerts; v++)
    {
        std::getline(fileObj, lineString);
        std::stringstream lineStream(lineString);

        lineStream >> vid >> vert.x >> vert.y >> vert.z;
        if(vid - 1 >= static_cast<int>(verts.size()))
        {
            std::cerr << "FAIL: vert " << vid - 1 << " >= "
                      << verts.size() << std::endl;
            fileObj.close();
            return false;
        }
        verts[vid - 1] = vert;
    }

    // strip start of elements section
    std::getline(fileObj, lineString);
    if(lineString.size() < 1 || lineString[0] != '$')
    {
        std::cerr << "FAIL: expected start of elements; got \""
                  << lineString << "\"" << std::endl;
        fileObj.close();
        return false;
    }

    // read elements
    int eid;
    Tuple3ui tri;
    for(int e = 0; e < numElements; e++)
    {
        std::getline(fileObj, lineString);
        std::stringstream lineStream(lineString);

        lineStream >> eid >> tri.x >> tri.y >> tri.z;
        if(eid - 1 >= static_cast<int>(tris.size()))
        {
            std::cerr << "FAIL: tri " << eid - 1 << " >= "
                      << tris.size() << std::endl;
            fileObj.close();
            return false;
        }

        tri.x--; tri.y--; tri.z--;
        tris[eid - 1] = tri;
    }

    // close and return
    if(!fileObj.good())
    {
        std::cerr << "FAIL: error at some point" << std::endl;
        fileObj.close();
        return false;
    }

    fileObj.close();

    std::cout << "FBEM input loaded successfully." << std::endl;
    return true;
}

//##############################################################################
//##############################################################################
bool FBemReader::
CheckFBemInputAgainstMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh, const std::string &fBemInputFile)
{
    std::cout << "check input\n";
    std::vector<Point3d> vertices; 
    std::vector<Tuple3ui> triangles; 
    const bool readSuccess = _ReadFBemInputToGeometry(fBemInputFile.c_str(), vertices, triangles); 
    std::cout << "Check if FBEM input has the same mesh as the current one using tolerance: " << _meshCheckTolerance << std::endl;
}

//##############################################################################
//##############################################################################
bool FBemReader::
ReadFBemOutputToInfo(std::shared_ptr<BEMSolutionMode> &solution, const std::string &fBemOutputFile)
{
    std::cout << "check output\n";
}
