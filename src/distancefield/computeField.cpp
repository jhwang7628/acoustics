#include <float.h>

#include <geometry/TriangleMesh.hpp>

#include "distanceField.h"

#ifdef COMPUTE_SIGNED_FIELD
#ifdef COMPUTE_CLOSEST_POINT
int ClosestPointField::computeSignedField(const std::string& filename, int resolution_g, int maxTriCount_g, int maxDepth_g)
#else
int DistanceField::computeSignedField(const std::string& filename, int resolution_g, int maxTriCount_g, int maxDepth_g)
#endif
#else
#ifdef COMPUTE_CLOSEST_POINT
int ClosestPointField::computeUnsignedField(const std::string& filename, int resolution_g, int maxTriCount_g, int maxDepth_g)
#else
int DistanceField::computeUnsignedField(const std::string& filename, int resolution_g, int maxTriCount_g, int maxDepth_g)
#endif
#endif
{
    ObjFile * objFile = new ObjFile(filename);

#ifdef COMPUTE_SIGNED_FIELD
#ifdef COMPUTE_CLOSEST_POINT
    int code =  computeSignedField(objFile, resolution_g, maxTriCount_g, maxDepth_g);
#else
    int code = computeSignedField(objFile, resolution_g, maxTriCount_g, maxDepth_g);
#endif
#else
#ifdef COMPUTE_CLOSEST_POINT
    int code = computeUnsignedField(objFile, resolution_g, maxTriCount_g, maxDepth_g);
#else
    int code = computeUnsignedField(objFile, resolution_g, maxTriCount_g, maxDepth_g);
#endif
#endif

    delete(objFile);
    return code;
}

#ifdef COMPUTE_SIGNED_FIELD
#ifdef COMPUTE_CLOSEST_POINT
int ClosestPointField::computeSignedField(ObjFile * objFileIn, int resolution_g, int maxTriCount_g, int maxDepth_g)
#else
int DistanceField::computeSignedField(ObjFile * objFileIn, int resolution_g, int maxTriCount_g, int maxDepth_g)
#endif
#else
#ifdef COMPUTE_CLOSEST_POINT
int ClosestPointField::computeUnsignedField(ObjFile * objFileIn, int resolution_g, int maxTriCount_g, int maxDepth_g)
#else
int DistanceField::computeUnsignedField(ObjFile * objFileIn, int resolution_g, int maxTriCount_g, int maxDepth_g)
#endif
#endif
{
    ObjFile objFile(*objFileIn);

    maxTriCount = maxTriCount_g; 
    maxDepth = maxDepth_g;
    resolution = resolution_g;

#ifdef COMPUTE_SIGNED_FIELD
    // === check if closed mesh
    if (!(objFile.isTriangularMesh()))
    {
        printf("Mesh was not triangular. Triangulating..."); fflush(NULL);
        objFile.triangulate();
        printf(" done\n");
    }
    ObjFileOrientable * objFileOrientable = new ObjFileOrientable(&objFile, 0);

    int numOrientationFlips = objFileOrientable->GenerateHalfEdgeDataStructure();

    cout << "Number of distinct connected components: " << objFileOrientable->nConnectedComponents() << endl;

    cout << "Checking if mesh has no boundary..." << endl;
    // check if mesh has no boundary
    if (objFileOrientable->hasBoundary())
    {
        cout << "Error: mesh has boundary. Signed distance field is ill-defined." << endl;
        cout << "  Num boundary edges: " << objFileOrientable->numBoundaryEdges() << endl;
        int edge = objFileOrientable->boundaryEdge(0);
        cout << "  A boundary edge: " << objFile.vertexPosition(objFileOrientable->halfEdge(edge).startVertex()) << " " <<
            objFile.vertexPosition(objFileOrientable->halfEdge(edge).endVertex()) << endl;

        return 1;
    }

    cout << "Mesh has no boundary (i.e. is closed surface)." << endl;

    cout << "Checking if input mesh is oriented consistently..."; 
    if (numOrientationFlips != 0)
    {
        cout << " no." << endl;
        cout << "Error: triangles in the input mesh are not oriented consistently." << endl;
        return 2;
    }
    else
        cout << " yes." << endl;

    delete(objFileOrientable);
#endif

    // === build octree

    printf("Preparing to build the octree. Max triangle count per cell: %d . Max depth: %d .\n", maxTriCount, maxDepth);fflush(NULL);

#ifdef COMPUTE_SIGNED_FIELD
    ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals> * objFileOctree = 
        new ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals>( &objFile, maxTriCount, maxDepth );
#else
    ObjFileOctree<TriangleWithCollisionInfo> * objFileOctree =
        new ObjFileOctree<TriangleWithCollisionInfo>( &objFile, maxTriCount, maxDepth );
#endif

    if (useAutomaticBox)
    {
        FieldBoundingBox bboxTight = objFileOctree->boundingBox();
        if (allBoxSidesEqual)
            bboxTight.regularize();
        bboxTight.expand(expansionRatio);
        bmax_ = bboxTight.bmax();
        bmin_ = bboxTight.bmin();
    }

    side = bmax_ - bmin_;
    gridX = side[0] / resolution;
    gridY = side[1] / resolution;
    gridZ = side[2] / resolution;
    setInvGridXYZ();

    double radius = len(side);

    cout << "Computing the distance field..." << endl;
    cout << "Corners of the bounding box are:" << endl;
    cout << "  " << bmin_ << endl;
    cout << "  " << bmax_ << endl;

    // do the zig-zag
    int i=0, j=0, k=0;

    int diri=1;
    int dirj=1;

    distanceData = (float*) malloc (sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));

    // triangleList will hold the triangles retrieved by range query
#ifdef COMPUTE_SIGNED_FIELD
    vector< TriangleWithCollisionInfoAndPseudoNormals* > triangleList;
#else
    vector< TriangleWithCollisionInfo* > triangleList;
#endif

    // debug
#ifdef GENERATE_DEBUG_DATA
    pseudoData = (float*) malloc (6*sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
#endif

#ifdef COMPUTE_CLOSEST_POINT

    closestPointData = (float*) malloc (3*sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
    this->closestFeatureData = new Feature[(resolution+1)*(resolution+1)*(resolution+1)];
#endif

    do
    {
        Vec3d currentPosition
            (bmin_[0] + 1.0 * i * gridX,
             bmin_[1] + 1.0 * j * gridY, 
             bmin_[2] + 1.0 * k * gridZ);

        //printf("%d %d %d\n",i,j,k);
        // process grid point i,j,k 
        triangleList.clear();

        while (triangleList.size() <= 0) // we should only go through this loop exactly once
        {
            //printf("radius: %f pos: %f %f %f\n", radius, currentPosition[0], currentPosition[1], currentPosition[2]);
            objFileOctree->rangeQuery(triangleList, Sphere(currentPosition, radius));

            if (triangleList.size() <= 0) 
            { // should never happen... but to stay robust
                cout << "Warning: range query didn't find any triangles. Incresing radius by a factor of 2 and re-trying." << endl;
                radius *= 2;
            }

            //if (triangleList.size() > 5)
            //makeUniqueList(triangleList);
        }    

        // find closest triangle among the retrieved ones
        // initialization:
        double closestDistance2 = 2.0 * radius * radius; // there will be somebody within that radius (actually even without the factor "2")
        // (factor "2" added to account for numerical round-off)

        int closestFeature = -1;
        int closestTriangle = -1;
        //int indexClosestTriangle = -1;

#ifdef COMPUTE_CLOSEST_POINT
        Vec3d closestPosition(DBL_MAX, DBL_MAX, DBL_MAX);
        Feature da_feature;
#endif

        for (unsigned int l=0; l<triangleList.size(); l++)
        {
            int closestLocalFeature = -1;
#ifdef COMPUTE_CLOSEST_POINT
            double alpha, beta, gamma;
            double d2 = triangleList[l]->distanceToPoint2(currentPosition, &closestLocalFeature, &alpha, &beta, &gamma);        
            if (d2 < closestDistance2)
            {
                // printf("%d -> (%d %d %d)\n", l,
                //                         triangleList.at(l)->firstIndex(),
                //                         triangleList.at(l)->secondIndex(),
                //                         triangleList.at(l)->thirdIndex());
                closestDistance2 = d2;
                closestPosition = triangleList[l]->getBarycentricLocation(alpha, beta, gamma);          
                closestFeature = closestLocalFeature; 
                closestTriangle = l;
                //indexClosestTriangle =  triangleList[l]->index();

                //Compute the feature
                da_feature.index1 = triangleList.at(l)->firstIndex();
                da_feature.index2 = triangleList.at(l)->secondIndex();
                da_feature.index3 = triangleList.at(l)->thirdIndex();
                da_feature.alpha = alpha;
                da_feature.beta = beta;
                da_feature.gamma = gamma;
            }
#else
            double d2 = triangleList[l]->distanceToPoint2(currentPosition,&closestLocalFeature);
            if (d2 < closestDistance2)
            {
                closestDistance2 = d2;
                closestFeature = closestLocalFeature;
                closestTriangle = l;
                //indexClosestTriangle =  triangleList[l]->index();          
            }
#endif
        }

        if (closestFeature < 0) // should never happen
        {
            cout << "Internal error: did not find any triangle within the guaranteed radius." << endl;
            return 3;
        }

//#ifdef GENERATE_DEBUG_DATA      
        //if ((i==34) && (j==17) && ((k==44) || (k==45)))      
        //{
        //    printf("Closest index: %d Feature: %d Min dist2: %.15f\n\n", indexClosestTriangle, closestFeature, closestDistance2);
        //}      
        //if ((i==34) && (j==17) && ((k==44) || (k==45)))        
        //    printf("Grid location: %.15f %.15f %.15f\n", currentPosition[0], currentPosition[1], currentPosition[2]);    
//#endif

        // square root...
        float closestDistance = sqrt(closestDistance2);

#ifdef COMPUTE_SIGNED_FIELD
        // determine sign, as in [Baerentzen 2002]
        Vec3d pseudoNormal = triangleList[closestTriangle]->pseudoNormal(closestFeature);
        Vec3d pseudoClosestPosition = triangleList[closestTriangle]->pseudoClosestPosition(closestFeature);

        if (dot(pseudoNormal,currentPosition-pseudoClosestPosition) < 0) // interior, must flip sign
            closestDistance *= -1;

        //printf("%G ", closestDistance);
#else 
        closestTriangle ++; // just to get rid of the set but unused warning..
        closestTriangle --;
#endif

        // register result
        distanceData[k * (resolution+1) * (resolution+1) + j * (resolution+1) + i] = closestDistance;


#ifdef COMPUTE_CLOSEST_POINT
        int dataIndex = 3 * (k * (resolution+1) * (resolution+1) + j * (resolution+1) + i);
        closestPointData[dataIndex+0] = closestPosition[0];
        closestPointData[dataIndex+1] = closestPosition[1];
        closestPointData[dataIndex+2] = closestPosition[2];
        this->closestFeatureData[(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)] = da_feature;
#endif

        // store debug info to disk file
#ifdef GENERATE_DEBUG_DATA
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+0] = pseudoNormal[0];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+1] = pseudoNormal[1];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+2] = pseudoNormal[2];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+3] = pseudoClosestPosition[0];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+4] = pseudoClosestPosition[1];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+5] = pseudoClosestPosition[2];
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+0] = 0;
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+1] = 0;
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+2] = 0;
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+3] = closestPosition[0];
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+4] = closestPosition[1];
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+5] = closestPosition[2];
#endif

        int oldi = i;
        int oldj = j;
        int oldk = k;

        // increment i,j,k
        i += diri; 

        if (i>resolution)
        {
            diri = -1; // start going to the left
            i = resolution;
            j += dirj;
        }

        if (i<0)
        {
            diri = +1; // start going to the right
            i = 0;
            j += dirj;
        }

        if (j>resolution)
        {
            dirj = -1; // start going down
            j = resolution;
            cout << k << " " << flush;
            k += 1;
        }

        if (j<0)
        {
            dirj = +1; // start going  up
            j = 0;
            cout << k << " " << flush;
            k += 1;
        }

        // temporal coherence: update radius for next iteration
        double delta=0.0; // this will be modified below
        if (i != oldi)
            delta = 1.01 * gridX;

        if (j != oldj)
            delta = 1.01 * gridY;

        if (k != oldk)
            delta = 1.01 * gridZ;

        radius = fabs(closestDistance) + delta;

    }
    while (k <= resolution);

#ifdef COMPUTE_SIGNED_FIELD
    cout << endl << "Signed distance field successfully computed..." << endl;
#else
    cout << endl << "Unsigned distance field successfully computed..." << endl;
#endif

    return 0;
}


#ifndef COMPUTE_SIGNED_FIELD
#ifdef COMPUTE_CLOSEST_POINT
int ClosestPointField::computeUnsignedField(TriangleMesh<REAL> *mesh, int resolution_g,
                                            int maxTriCount_g, int maxDepth_g)
#else
int DistanceField::computeUnsignedField(TriangleMesh<REAL> *mesh, int resolution_g,
                                        int maxTriCount_g, int maxDepth_g)
#endif
{
    //ObjFile objFile(*objFileIn);

    maxTriCount = maxTriCount_g; 
    maxDepth = maxDepth_g;
    resolution = resolution_g;

    // === build octree

    printf("Preparing to build the octree. Max triangle count per cell: %d . Max depth: %d .\n", maxTriCount, maxDepth);fflush(NULL);

    ObjFileOctree<TriangleWithCollisionInfo> * objFileOctree =
        new ObjFileOctree<TriangleWithCollisionInfo>( mesh, maxTriCount, maxDepth );

    if (useAutomaticBox)
    {
        FieldBoundingBox bboxTight = objFileOctree->boundingBox();
        if (allBoxSidesEqual)
            bboxTight.regularize();
        bboxTight.expand(expansionRatio);
        bmax_ = bboxTight.bmax();
        bmin_ = bboxTight.bmin();
    }

    side = bmax_ - bmin_;
    gridX = side[0] / resolution;
    gridY = side[1] / resolution;
    gridZ = side[2] / resolution;
    setInvGridXYZ();

    double radius = len(side);

    cout << "Computing the distance field..." << endl;
    cout << "Corners of the bounding box are:" << endl;
    cout << "  " << bmin_ << endl;
    cout << "  	" << bmax_ << endl;

    // do the zig-zag
    int i=0, j=0, k=0;

    int diri=1;
    int dirj=1;

    distanceData = (float*) malloc (sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));

    // triangleList will hold the triangles retrieved by range query
    vector< TriangleWithCollisionInfo* > triangleList;

    // debug
#ifdef GENERATE_DEBUG_DATA
    pseudoData = (float*) malloc (6*sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
#endif

#ifdef COMPUTE_CLOSEST_POINT
    closestPointData = (float*) malloc
        (3*sizeof(float)*(resolution+1)*(resolution+1)*(resolution+1));
    closestFeatureData
        = new Feature[(resolution+1)*(resolution+1)*(resolution+1)];
#endif

    do
    {
        Vec3d currentPosition
            (bmin_[0] + 1.0 * i * gridX,
             bmin_[1] + 1.0 * j * gridY, 
             bmin_[2] + 1.0 * k * gridZ);

        //printf("%d %d %d\n",i,j,k);
        // process grid point i,j,k 
        triangleList.clear();

        while (triangleList.size() <= 0) // we should only go through this loop exactly once
        {
            //printf("radius: %f pos: %f %f %f\n", radius, currentPosition[0], currentPosition[1], currentPosition[2]);
            objFileOctree->rangeQuery(triangleList, Sphere(currentPosition, radius));

            if (triangleList.size() <= 0) 
            { // should never happen... but to stay robust
                cout << "Warning: range query didn't find any triangles. Incresing radius by a factor of 2 and re-trying." << endl;
                radius *= 2;
            }

            //if (triangleList.size() > 5)
            //makeUniqueList(triangleList);
        }		

        // find closest triangle among the retrieved ones
        // initialization:
        double closestDistance2 = 2.0 * radius * radius; // there will be somebody within that radius (actually even without the factor "2")
        // (factor "2" added to account for numerical round-off)

        int closestFeature = -1;
        //int closestTriangle = -1;
        //int indexClosestTriangle = -1;

#ifdef COMPUTE_CLOSEST_POINT
        Vec3d closestPosition(DBL_MAX, DBL_MAX, DBL_MAX);
        ClosestPointField::Feature closestBaryFeature;
#endif

        for (unsigned int l=0; l<triangleList.size(); l++)
        {
            int closestLocalFeature = -1;				

#ifdef COMPUTE_CLOSEST_POINT
            double alpha, beta, gamma;
            double d2 = triangleList[l]->distanceToPoint2(currentPosition, &closestLocalFeature, &alpha, &beta, &gamma);				
            if (d2 < closestDistance2)
            {					
                closestDistance2 = d2;
                closestPosition = triangleList[l]->getBarycentricLocation(alpha, beta, gamma);					
                closestFeature = closestLocalFeature; 
                //closestTriangle = l;
                //indexClosestTriangle =	triangleList[l]->index();					

                closestBaryFeature.index1 = triangleList.at(l)->firstIndex();
                closestBaryFeature.index2 = triangleList.at(l)->secondIndex();
                closestBaryFeature.index3 = triangleList.at(l)->thirdIndex();
                closestBaryFeature.alpha = alpha;
                closestBaryFeature.beta = beta;
                closestBaryFeature.gamma = gamma;
            }
#else
            double d2 = triangleList[l]->distanceToPoint2(currentPosition,&closestLocalFeature);
            if (d2 < closestDistance2)
            {
                closestDistance2 = d2;
                closestFeature = closestLocalFeature;
                //closestTriangle = l;
                //indexClosestTriangle =	triangleList[l]->index();					
            }
#endif
        }

        if (closestFeature < 0) // should never happen
        {
            cout << "Internal error: did not find any triangle within the guaranteed radius." << endl;
            return 3;
        }

//#ifdef GENERATE_DEBUG_DATA			
//        if ((i==34) && (j==17) && ((k==44) || (k==45)))			
//        {
//            printf("Closest index: %d Feature: %d Min dist2: %.15f\n\n", indexClosestTriangle, closestFeature, closestDistance2);
//        }			
//        if ((i==34) && (j==17) && ((k==44) || (k==45)))				
//            printf("Grid location: %.15f %.15f %.15f\n", currentPosition[0], currentPosition[1], currentPosition[2]);		
//#endif

        // square root...
        float closestDistance = sqrt(closestDistance2);

        // register result
        distanceData[k * (resolution+1) * (resolution+1) + j * (resolution+1) + i] = closestDistance;


#ifdef COMPUTE_CLOSEST_POINT
        int dataIndex = 3 * (k * (resolution+1) * (resolution+1) + j * (resolution+1) + i);
        closestPointData[dataIndex+0] = closestPosition[0];
        closestPointData[dataIndex+1] = closestPosition[1];
        closestPointData[dataIndex+2] = closestPosition[2];

        int featureIndex = (k * (resolution+1) * (resolution+1) + j * (resolution+1) + i);
        closestFeatureData[featureIndex] = closestBaryFeature;
#endif

        // store debug info to disk file
#ifdef GENERATE_DEBUG_DATA
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+0] = pseudoNormal[0];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+1] = pseudoNormal[1];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+2] = pseudoNormal[2];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+3] = pseudoClosestPosition[0];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+4] = pseudoClosestPosition[1];
        //pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+5] = pseudoClosestPosition[2];
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+0] = 0;
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+1] = 0;
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+2] = 0;
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+3] = closestPosition[0];
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+4] = closestPosition[1];
        pseudoData[6*(k * (resolution+1) * (resolution+1) + j * (resolution+1) + i)+5] = closestPosition[2];
#endif

        int oldi = i;
        int oldj = j;
        int oldk = k;

        // increment i,j,k
        i += diri; 

        if (i>resolution)
        {
            diri = -1; // start going to the left
            i = resolution;
            j += dirj;
        }

        if (i<0)
        {
            diri = +1; // start going to the right
            i = 0;
            j += dirj;
        }

        if (j>resolution)
        {
            dirj = -1; // start going down
            j = resolution;
            cout << k << " " << flush;
            k += 1;
        }

        if (j<0)
        {
            dirj = +1; // start going	up
            j = 0;
            cout << k << " " << flush;
            k += 1;
        }

        // temporal coherence: update radius for next iteration
        double delta=0.0; // this will be modified below
        if (i != oldi)
            delta = 1.01 * gridX;

        if (j != oldj)
            delta = 1.01 * gridY;

        if (k != oldk)
            delta = 1.01 * gridZ;

        radius = fabs(closestDistance) + delta;

    }
    while (k <= resolution);

    cout << endl << "Unsigned distance field successfully computed..." << endl;

    return 0;
}
#endif
