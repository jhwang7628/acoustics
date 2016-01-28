#include "objfileOctree.h"
#include "triple.h"

#if 0
#include <mesh/Triangle.h>
#include <mesh/TriMesh.h>
#endif

// builds an octree on top of the obj file

template<class TriangleClass>
const double ObjFileOctree<TriangleClass>::bboxExpansionRatio = 1.05;

    template<class TriangleClass>
ObjFileOctree<TriangleClass>::ObjFileOctree( const std::string& filename, int maxTriCount_g, int maxDepth_g, int printInfo )
{
    ObjFile * objFile = new ObjFile(filename);
    ObjFileOctree(objFile, maxTriCount_g, maxDepth_g, printInfo);
    delete(objFile);
}

template<class TriangleClass>
ObjFileOctree<TriangleClass>::ObjFileOctree( ObjFile * objFileIn, int maxTriCount_g, int maxDepth_g, int printInfo ) : 
    maxTriCount(maxTriCount_g), maxDepth(maxDepth_g)
{
    // create triangles structure
    //ObjFile * objFile = new ObjFile(filename);
    ObjFile * objFile = new ObjFile(*objFileIn);

    cout << "Checking if mesh is triangular... ";
    if (!objFile->isTriangularMesh())
    {
        cout << "mesh was not triangular: triangulating... ";
        objFile->triangulate();
        cout << "done" << endl;
    }
    else
        cout << "yes" << endl;

    int triangleIndex = 0;
    triangles.clear();
    for(unsigned int i=0; i < objFile->numGroups(); i++) // over all groups
        for (unsigned int j=0; j < (objFile->groupHandle(i))->faceCount(); j++) // over all faces
        {
            Vec3d p0 = objFile->vertexPosition(objFile->vertexIndex(i,j,0));
            Vec3d p1 = objFile->vertexPosition(objFile->vertexIndex(i,j,1));
            Vec3d p2 = objFile->vertexPosition(objFile->vertexIndex(i,j,2));
            TriangleClass triangle(p0,p1,p2);
            triangle.setIndex(triangleIndex); // 0-indexed
            triangle.setVertexIndex(objFile->vertexIndex(i,j,0), objFile->vertexIndex(i,j,1), objFile->vertexIndex(i,j,2));
            triangleIndex++;
            triangles.push_back(triangle);
        }

    cout << "Total number of triangles is: " << triangles.size() << endl;

    // build the octree
    cout << "Building the octree data structure... " << endl;

    root = new Octree<TriangleClass>(maxDepth); 

    root->setBuildPrintInfo(printInfo);

    //root->build(triangles,bboxOctree,maxTriCount);
    root->build(triangles,maxTriCount);

    if(printInfo == 1)
    {
        int numMaxDepthExceededCases;
        int numMaxTriInDepthExceededCases;
        root->getBuildInfo(&numMaxDepthExceededCases, &numMaxTriInDepthExceededCases);
        printf("Total number of cells with more than %d triangles: %d. Max triangles in such cells: %d.\n", maxTriCount, numMaxDepthExceededCases,
                numMaxTriInDepthExceededCases);
    }

    triangles.clear(); // release memory
    delete(objFile);

    //cout << "done" << endl;
}

template<class TriangleClass>
ObjFileOctree<TriangleClass>::ObjFileOctree( TriangleMesh<REAL> *mesh, int maxTriCount_g,
                                             int maxDepth_g, int printInfo )
    : maxTriCount(maxTriCount_g), maxDepth(maxDepth_g)
{
    int triangleIndex = 0;
    triangles.clear();
#if 0
    for(unsigned int i=0; i < objFile->numGroups(); i++) // over all groups
        for (unsigned int j=0; j < (objFile->groupHandle(i))->faceCount(); j++) // over all faces
        {
            Vec3d p0 = objFile->vertexPosition(objFile->vertexIndex(i,j,0));
            Vec3d p1 = objFile->vertexPosition(objFile->vertexIndex(i,j,1));
            Vec3d p2 = objFile->vertexPosition(objFile->vertexIndex(i,j,2));
            TriangleClass triangle(p0,p1,p2);
            triangle.setIndex(triangleIndex); // 0-indexed
            triangleIndex++;
            triangles.push_back(triangle);
        }
#endif
    for (unsigned int i = 0; i < mesh->triangles().size(); i++)
    {
        const Tuple3ui &t = mesh->triangle_ids( i );
#if 0
        Triangle *t = meshTriangles.at(i);
        const VEC3F &x0 = t->getX(0);
        const VEC3F &x1 = t->getX(1);
        const VEC3F &x2 = t->getX(2);
#endif
        const Point3d &x0 = mesh->vertex( t[ 0 ] );
        const Point3d &x1 = mesh->vertex( t[ 1 ] );
        const Point3d &x2 = mesh->vertex( t[ 2 ] );

        Vec3d p0( x0[0], x0[1], x0[2] );
        Vec3d p1( x1[0], x1[1], x1[2] );
        Vec3d p2( x2[0], x2[1], x2[2] );

        TriangleClass triangle(p0,p1,p2);
        triangle.setIndex(triangleIndex); // 0-indexed

        // Set the vertex indices for future reference.  These
        // refer to indices in the TriMesh object
#if 0
        triangle.setVertexIndex( t->getIndex(0), t->getIndex(1), t->getIndex(2) );
#endif
        triangle.setVertexIndex( t[ 0 ], t[ 1 ], t[ 2 ] );

        triangleIndex++;
        triangles.push_back(triangle);
    }

    cout << "Total number of triangles is: " << triangles.size() << endl;

    // build the octree
    cout << "Building the octree data structure... " << endl;

    root = new Octree<TriangleClass>(maxDepth); 

    root->setBuildPrintInfo(printInfo);

    //root->build(triangles,bboxOctree,maxTriCount);
    root->build(triangles,maxTriCount);

    if(printInfo == 1)
    {
        int numMaxDepthExceededCases;
        int numMaxTriInDepthExceededCases;
        root->getBuildInfo(&numMaxDepthExceededCases, &numMaxTriInDepthExceededCases);
        printf("Total number of cells with more than %d triangles: %d. Max triangles in such cells: %d.\n", maxTriCount, numMaxDepthExceededCases,
                numMaxTriInDepthExceededCases);
    }

    triangles.clear(); // release memory
    //cout << "done" << endl;
}

template<>
ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals>::ObjFileOctree( ObjFile * objFileIn, int maxTriCount_g, int maxDepth_g, int printInfo ) : 
    maxTriCount(maxTriCount_g), maxDepth(maxDepth_g)
{
    // create triangles structure
    ObjFile * objFile = new ObjFile(*objFileIn);
    //ObjFile * objFile = new ObjFile(filename,0);

    cout << "Checking if mesh is triangular... ";
    if (!objFile->isTriangularMesh())
    {
        cout << "mesh was not triangular: triangulating... ";
        objFile->triangulate();
        cout << "done" << endl;
    }
    else
        cout << "yes" << endl;

    ObjFile * pseudoNormalObjFile = objFile;
    //let's not hack
    //ObjFile * pseudoNormalObjFile = new ObjFile("/nfs/hn18/sim/barbic/vps/twoHoses/deformed.mc.obj");

    int numZeroAreaFaces = pseudoNormalObjFile->removeZeroAreaFaces();
    cout << "Encountered and removed " << numZeroAreaFaces << " faces with zero surface area." << endl;

    triangles.clear();


    //unsigned int n = pseudoNormalObjFile->numVertices();
    pseudoNormalObjFile->computePseudoNormals();
    pseudoNormalObjFile->computeEdgePseudoNormals();

    // compute face pseudonormals
    Vec3d pseudoNormals[7];
    for(unsigned int i=0; i < pseudoNormalObjFile->numGroups(); i++) // over all groups
    {
        for (unsigned int j=0; j < (pseudoNormalObjFile->groupHandle(i))->faceCount(); j++) // over all faces
        {
            // vertices

            unsigned int index0 = pseudoNormalObjFile->vertexIndex(i,j,0);
            unsigned int index1 = pseudoNormalObjFile->vertexIndex(i,j,1);
            unsigned int index2 = pseudoNormalObjFile->vertexIndex(i,j,2);

            pseudoNormals[0] = pseudoNormalObjFile->pseudoNormal(index0);
            pseudoNormals[1] = pseudoNormalObjFile->pseudoNormal(index1);
            pseudoNormals[2] = pseudoNormalObjFile->pseudoNormal(index2);

            // edge pseudo normal
            if (pseudoNormalObjFile->edgePseudoNormal(index0,index1,&pseudoNormals[3]) != 0)
            {
                cout << "Error: encountered an edge without a pseudonormal. Degenerate face? Vertices: " << index0 << " " << index1 << endl;
                exit(1);
            }
            if (pseudoNormalObjFile->edgePseudoNormal(index1,index2,&pseudoNormals[4]) != 0)
            {
                cout << "Error: encountered an edge without a pseudonormal. Degenerate face? Vertices: " << index1 << " " << index2 << endl;
                exit(1);
            }
            if (pseudoNormalObjFile->edgePseudoNormal(index2,index0,&pseudoNormals[5]) != 0)
            {
                cout << "Error: encountered an edge without a pseudonormal. Degenerate face? Vertices: " << index2 << " " << index0 << endl;
                exit(1);
            }

            // face pseudo normal
            Vec3d p0 = pseudoNormalObjFile->vertexPosition(index0);
            Vec3d p1 = pseudoNormalObjFile->vertexPosition(index1);
            Vec3d p2 = pseudoNormalObjFile->vertexPosition(index2);

            pseudoNormals[6] = norm(cross(p1-p0,p2-p0)); 

            for(int normali=0; normali < 7; normali++)
            {
                if (my_isnan(pseudoNormals[normali][0]) || 
                        my_isnan(pseudoNormals[normali][1]) || 
                        my_isnan(pseudoNormals[normali][2]))
                {
                    cout << "Error: nan encountered: " << pseudoNormals[normali][0] << " " << pseudoNormals[normali][1] << " " << pseudoNormals[normali][2] << endl;
                    cout << "Group: " << i << " Triangle: " << j << " " << endl;
                    cout << "  vtx0: " << index0 << " vtx1: " << index1 << " vtx2: " << index2 << endl;
                    cout << "  "  << p0 << endl;
                    cout << "  "  << p1 << endl;
                    cout << "  "  << p2 << endl;
                    cout << "Feature: " << normali << endl;
                    exit(1);
                }
            }

            // hack
            //p0 = objFile->vertexPosition(index0);
            //p1 = objFile->vertexPosition(index1);
            //p2 = objFile->vertexPosition(index2);
            TriangleWithCollisionInfoAndPseudoNormals triangle(p0,p1,p2,pseudoNormals);
            triangle.setVertexIndex(objFile->vertexIndex(i,j,0), objFile->vertexIndex(i,j,1), objFile->vertexIndex(i,j,2));
            triangles.push_back(triangle);
        }
    }

    cout << "Total number of triangles is: " << triangles.size() << endl;

    // build the octree

    Vec3d bmin, bmax;

    objFile->boundingBoxCube(1.0,bmin,bmax);

    cout << "The model tight-fitting cube-shaped bounding box is: " << endl;
    cout << "xmin: " << bmin[0] << " xmax: " << bmax[0] << endl;
    cout << "ymin: " << bmin[1] << " ymax: " << bmax[1] << endl;
    cout << "zmin: " << bmin[2] << " zmax: " << bmax[2] << endl;

    FieldBoundingBox bboxOctree(bmin,bmax);
    bboxOctree.expand(bboxExpansionRatio);

    cout << "Starting the octree creation algorithm..." << endl;

    root = new Octree<TriangleWithCollisionInfoAndPseudoNormals>(maxDepth); 

    root->setBuildPrintInfo(printInfo);
    root->build(triangles,bboxOctree,maxTriCount);

    if(printInfo == 1)
    {
        int numMaxDepthExceededCases;
        int numMaxTriInDepthExceededCases;
        root->getBuildInfo(&numMaxDepthExceededCases, &numMaxTriInDepthExceededCases);
        printf("Total number of cells with more than %d triangles: %d. Max triangles in such cells: %d.\n", maxTriCount, numMaxDepthExceededCases,
                numMaxTriInDepthExceededCases);
    }

    triangles.clear(); // release memory
    delete(objFile);

    cout << "Octree creation completed successfully." << endl;

}

template ObjFileOctree<TriangleBasic>::ObjFileOctree( ObjFile * objFile, int maxTriCount_g, int maxDepth_g, int printInfo );  
template ObjFileOctree<TriangleWithCollisionInfo>::ObjFileOctree( ObjFile * objFile, int maxTriCount_g, int maxDepth_g, int printInfo );  
template ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals>::ObjFileOctree( ObjFile * objFile, int maxTriCount_g, int maxDepth_g, int printInfo );  

template ObjFileOctree<TriangleBasic>::ObjFileOctree( const std::string& filename, int maxTriCount_g, int maxDepth_g, int printInfo );  
template ObjFileOctree<TriangleWithCollisionInfo>::ObjFileOctree( const std::string& filename, int maxTriCount_g, int maxDepth_g, int printInfo );  
template ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals>::ObjFileOctree( const std::string& filename, int maxTriCount_g, int maxDepth_g, int printInfo );  

template ObjFileOctree<TriangleBasic>::ObjFileOctree( TriangleMesh<REAL> *mesh,
        int maxTriCount_g, int maxDepth_g, int printInfo );  
template ObjFileOctree<TriangleWithCollisionInfo>::ObjFileOctree(
        TriangleMesh<REAL> *mesh, int maxTriCount_g, int maxDepth_g,
        int printInfo );  
