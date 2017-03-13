#include <cassert>
#include <cfloat>
#include <iostream>
#include <fstream>

#include "trilinearInterpolation.h"
#include "objfile.h"
#include "objfileOctree.h"
#include "AdaptiveDistanceField.h"

// for testing
#include "utils/SimpleTimer.h"
#if PNG
#include "lodepng.h"
#endif


double AdaptiveDistanceField::BuildNode::distance(const Vector3d& pt, const Vector3d& cellBmin, const Vector3d& cellSize) const
{
    double wx = (pt.x - cellBmin.x) / cellSize.x;
    double wy = (pt.y - cellBmin.y) / cellSize.y;
    double wz = (pt.z - cellBmin.z) / cellSize.z;
    return TRILINEAR_INTERPOLATION(wx, wy, wz,
                                   _closestFeatures[0]._distance, _closestFeatures[1]._distance,
                                   _closestFeatures[3]._distance, _closestFeatures[2]._distance,
                                   _closestFeatures[4]._distance, _closestFeatures[5]._distance,
                                   _closestFeatures[7]._distance, _closestFeatures[6]._distance);
}

/*Vector3d AdaptiveDistanceField::BuildNode::gradient(const Vector3d& pt, const Vector3d& cellBmin, const Vector3d& cellSize) const
{
    double wx = (pt.x - cellBmin.x) / cellSize.x;
    double wy = (pt.y - cellBmin.y) / cellSize.y;
    double wz = (pt.z - cellBmin.z) / cellSize.z;
    return Vector3d(
                    GRADIENT_COMPONENT_X(wx, wy, wz,
                                         _closestFeatures[0]._distance, _closestFeatures[1]._distance,
                                         _closestFeatures[3]._distance, _closestFeatures[2]._distance,
                                         _closestFeatures[4]._distance, _closestFeatures[5]._distance,
                                         _closestFeatures[7]._distance, _closestFeatures[6]._distance,
                                         cellSize.x),
                    GRADIENT_COMPONENT_Y(wx, wy, wz,
                                         _closestFeatures[0]._distance, _closestFeatures[1]._distance,
                                         _closestFeatures[3]._distance, _closestFeatures[2]._distance,
                                         _closestFeatures[4]._distance, _closestFeatures[5]._distance,
                                         _closestFeatures[7]._distance, _closestFeatures[6]._distance,
                                         cellSize.y),
                    GRADIENT_COMPONENT_Z(wx, wy, wz,
                                         _closestFeatures[0]._distance, _closestFeatures[1]._distance,
                                         _closestFeatures[3]._distance, _closestFeatures[2]._distance,
                                         _closestFeatures[4]._distance, _closestFeatures[5]._distance,
                                         _closestFeatures[7]._distance, _closestFeatures[6]._distance,
                                         cellSize.z)
                    );
}*/


double AdaptiveDistanceField::QueueNode::distance(const Vector3d& pt) const
{
    return _node->distance(pt, _cellBmin, _cellSize);
}


AdaptiveDistanceField::AdaptiveDistanceField()
 : 
     _numNodesInProgress(0), _numNodes(0)
 { }
 
 
void AdaptiveDistanceField::computeSignedField(const std::string& filename, double subdivideRadius,
                                                int maxAdfLevel, double maxError, int tileLevels,
                                                int maxTriCount, const int maxDepth, bool allBoxSidesEqual, double expansionRatio,
                                                int numWorkerThreads)
{
    ObjFile* of = new ObjFile(filename);
    ObjFile& objFile = *of;
    
    // === check if closed mesh
    if (!(objFile.isTriangularMesh()))
    {
        printf("Mesh was not triangular. Triangulating..."); fflush(NULL);
        objFile.triangulate();
        printf(" done\n");
    }
    ObjFileOrientable * objFileOrientable = new ObjFileOrientable(&objFile, 0);

    int numOrientationFlips = objFileOrientable->GenerateHalfEdgeDataStructure();

    std::cout << "Number of distinct connected components: " << objFileOrientable->nConnectedComponents() << std::endl;

    std::cout << "Checking if mesh has no boundary..." << std::endl;
    // check if mesh has no boundary
    if (objFileOrientable->hasBoundary())
    {
        std::cout << "Error: mesh has boundary. Signed distance field is ill-defined." << std::endl;
        std::cout << "  Num boundary edges: " << objFileOrientable->numBoundaryEdges() << std::endl;
        int edge = objFileOrientable->boundaryEdge(0);
        std::cout << "  A boundary edge: " << objFile.vertexPosition(objFileOrientable->halfEdge(edge).startVertex()) << " " <<
            objFile.vertexPosition(objFileOrientable->halfEdge(edge).endVertex()) << std::endl;

        //return 1;
    }

    std::cout << "Mesh has no boundary (i.e. is closed surface)." << std::endl;

    std::cout << "Checking if input mesh is oriented consistently..."; 
    if (numOrientationFlips != 0)
    {
        std::cout << " no." << std::endl;
        std::cout << "Error: triangles in the input mesh are not oriented consistently." << std::endl;
        //return 2;
    }
    else
        std::cout << " yes." << std::endl;

    delete(objFileOrientable);
    
    
    // octree
       
    printf("Preparing to build the octree. Max triangle count per cell: %d . Max depth: %d .\n", maxTriCount, maxDepth);fflush(NULL);

    ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals> * objFileOctree = 
        new ObjFileOctree<TriangleWithCollisionInfoAndPseudoNormals>( &objFile, maxTriCount, maxDepth );

    Vec3d bmin_, bmax_;
    //if (useAutomaticBox)
    //{
    FieldBoundingBox bboxTight = objFileOctree->boundingBox();
    if (allBoxSidesEqual)
        bboxTight.regularize();
    bboxTight.expand(expansionRatio);
    bmax_ = bboxTight.bmax();
    bmin_ = bboxTight.bmin();
    //}
        
    Vector3d bmin(bmin_[0], bmin_[1], bmin_[2]);
    Vector3d bmax(bmax_[0], bmax_[1], bmax_[2]);
    double bboxDiagonal = (bmax - bmin).length();
    
    ClosestFeatureFunc closestFunc = [objFileOctree, bboxDiagonal]
                                     (const Vector3d& point, Feature* closest, bool hintValid, const Vector3d& hintPoint, const Feature& hintClosest)
    {
        Vec3d currentPosition(point.x, point.y, point.z);
        
        std::vector< TriangleWithCollisionInfoAndPseudoNormals*> triangleList;
        
        double radius = hintValid ? 
                        (point - hintPoint).length() + std::abs(hintClosest._distance) :
                        bboxDiagonal;
        
        while (triangleList.size() <= 0) // we should only go through this loop exactly once
        {
            //printf("radius: %f pos: %f %f %f\n", radius, currentPosition[0], currentPosition[1], currentPosition[2]);
            objFileOctree->rangeQuery(triangleList, Sphere(currentPosition, radius));

            if (triangleList.size() <= 0) 
            { // should never happen... but to stay robust
                std::cout << "Warning: range query didn't find any triangles. Increasing radius by a factor of 2 and re-trying." << std::endl;
                radius *= 2;
            }
        }
                
        double closestDistance2 = 2.0 * radius * radius;
        int closestFeature = -1;
        int closestTriangle = -1;
        for (size_t l=0; l<triangleList.size(); l++)
        {
            int closestLocalFeature = -1;
            double d2 = triangleList[l]->distanceToPoint2(currentPosition,&closestLocalFeature);
            if (d2 < closestDistance2)
            {
                closestDistance2 = d2;
                closestFeature = closestLocalFeature;
                closestTriangle = l;
                //indexClosestTriangle =  triangleList[l]->index();          
            }
        }
        assert(closestFeature >= 0);    // at least one triangle should be found
        
        // square root...
        double closestDistance = sqrt(closestDistance2);
        
        // determine sign, as in [Baerentzen 2002]
        Vec3d pseudoNormal = triangleList[closestTriangle]->pseudoNormal(closestFeature);
        Vec3d pseudoClosestPosition = triangleList[closestTriangle]->pseudoClosestPosition(closestFeature);
        if (dot(pseudoNormal,currentPosition-pseudoClosestPosition) < 0) // interior, must flip sign
            closestDistance *= -1.0;
        //printf("%G ", closestDistance);

        // register result
        closest->_distance = closestDistance;
    };

    computeSignedField(closestFunc, bmin, bmax, subdivideRadius, maxAdfLevel, maxError, tileLevels, numWorkerThreads);
}
 
 
 void AdaptiveDistanceField::computeSignedField(ClosestFeatureFunc closestFunc,
                                               const Vector3d& bmin, const Vector3d& bmax, double subdivideRadius,
                                               int maxAdfLevel, double maxError, int tileLevels, int numWorkerThreads)
{
    assert (2 <= tileLevels && tileLevels <= 10 && 1 <= maxAdfLevel && maxAdfLevel <= 31);
    assert(bmax.x >= bmin.x && bmax.y >= bmin.y && bmax.z >= bmin.z);
    
    _adfNodes.clear();
    _nodeQueue = std::queue<QueueNode>();
    _numNodesInProgress = 0;
    
    _subdivideRadius = subdivideRadius;
    _tileLevels = tileLevels;
    _maxAdfLevel = maxAdfLevel;
    _maxError = maxError;
    _closestFunc = closestFunc;
    _rootBmin = bmin;
    _rootBmax = bmax;
 
    // tile cache
    _tileDim = (1<<_tileLevels) + 1;
        
 SimpleTimer timer;
 timer.Start();
    
    // create root node; push into node queue
    BuildNode* root = new AdaptiveDistanceField::BuildNode;
    Vector3d _rootCellSize = _rootBmax - _rootBmin;
    computeNodeClosestFeatures(root, _rootBmin, _rootCellSize);
    _nodeQueue.push(QueueNode(root, 0, _rootBmin, _rootCellSize));
    _numNodes = 1;
    
    std::vector<std::thread> nodeQueueThreads(numWorkerThreads);
    for (int tt = 0; tt < numWorkerThreads; tt++)
    {
        nodeQueueThreads[tt] = std::thread(&AdaptiveDistanceField::nodeQueueThreadFunc, this);
    }
    for (int tt = 0; tt < numWorkerThreads; tt++)
    {
        nodeQueueThreads[tt].join();
    }
    
    timer.Pause();
    printf("\ntotal nodes: %d\n", (int)_numNodes);
    printf("\ntreebuild: %lf\n", timer.Duration());
    
    timer = SimpleTimer();
    timer.Start();
    
    coalesce(root);
    delete root;
    root = nullptr;
    
    timer.Pause();
    printf("coalesce: %lf\n", timer.Duration());
    
#if PNG
    std::string xyz[3] = {"x", "y", "z"};
    for (int n =0; n < 3; n++)
        savePng(n, "tree_dist_"+xyz[n]+".png");
#endif
}

// thread-pool worker function
void AdaptiveDistanceField::nodeQueueThreadFunc()
{
    uint32 tileFlatLength = _tileDim * _tileDim * _tileDim;
    Feature* featuresTile = new Feature[tileFlatLength];
    bool* validFlagsTile = new bool[tileFlatLength];
    
    while (true)
    {
        // grab a node from the queue
        QueueNode queueNode;
        {
            std::unique_lock<std::mutex> ul(_nodeQueueLock);
            while (_nodeQueue.empty())
            {
                if (_numNodesInProgress == 0)
                {
                    delete[] featuresTile;
                    delete[] validFlagsTile;
                    return;
                }
                _nodeQueueCv.wait(ul);
            }
            queueNode = _nodeQueue.front();
            _nodeQueue.pop();
            ++_numNodesInProgress;
        }
        
        // process the node
        std::memset(featuresTile, 0, tileFlatLength*sizeof(Feature));
        std::memset(validFlagsTile, 0, tileFlatLength*sizeof(bool));
        
        setTileCornersFromNodeClosestFeatures(featuresTile, validFlagsTile, queueNode._node);
        
        // just use the (0,0,0) corner of this node cell as the hint
        QueryHint queryHint(queueNode._cellBmin, featuresTile[0]);
        recursiveSubdivToMaxLevel(queueNode, Vector3u(0,0,0), _tileDim-1, featuresTile, validFlagsTile, queryHint);
        
        _nodeQueueLock.lock();
        --_numNodesInProgress;
        _nodeQueueCv.notify_all();
        _nodeQueueLock.unlock();
    }
}

void AdaptiveDistanceField::computeNodeClosestFeatures(BuildNode* node, const Vector3d& cellBmin, const Vector3d& cellSize) const
{
    //bool queryHintValid = false;
    //QueryHint queryHint;
    #pragma omp parallel for schedule(static) collapse(3)
    for (int k = 0; k < 2; k++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < 2; i++)
            {
                Vector3d point(cellBmin.x + i * cellSize.x,
                               cellBmin.y + j * cellSize.y,
                               cellBmin.z + k * cellSize.z);
                Feature* closestFeature = &(node->_closestFeatures[(k*2 + j)*2 + i]);
                _closestFunc(point, closestFeature, false, Vector3d(), Feature());//queryHintValid, queryHint._point, queryHint._closest);
                               
                //queryHintValid = true;
                //queryHint._point = point;
                //queryHint._closest = *closestFeature;
            }
        }
    }
}

void AdaptiveDistanceField::setNodeClosestFeaturesFromTile(BuildNode* node, Feature* featuresTile, bool* validFlagsTile,
                                                           const Vector3u& cellTileBmin, uint32 cellTileWidth) const
{
    for (uint32 k = 0; k < 2; k++)
    {
        for (uint32 j = 0; j < 2; j++)
        {
            for (uint32 i = 0; i < 2; i++)
            {
                Vector3u pointTileIndices(cellTileBmin.x + i * cellTileWidth,
                                          cellTileBmin.y + j * cellTileWidth,
                                          cellTileBmin.z + k * cellTileWidth);
                uint32 pointTileFlatIndex = (pointTileIndices[2]*_tileDim + pointTileIndices[1])*_tileDim + pointTileIndices[0];
                assert(validFlagsTile[pointTileFlatIndex]);
                node->_closestFeatures[(k*2 + j)*2 + i] = featuresTile[pointTileFlatIndex];
            }
        }
    }
}

void AdaptiveDistanceField::setTileCornersFromNodeClosestFeatures(Feature* featuresTile, bool* validFlagsTile, BuildNode* node) const
{
    uint32 tileWidth = _tileDim - 1;
    for (uint32 k = 0; k < 2; k++)
    {
        for (uint32 j = 0; j < 2; j++)
        {
            for (uint32 i = 0; i < 2; i++)
            {
                Vector3u pointTileIndices(i * tileWidth,
                                          j * tileWidth,
                                          k * tileWidth);
                uint32 pointTileFlatIndex = (pointTileIndices[2]*_tileDim + pointTileIndices[1])*_tileDim + pointTileIndices[0];
                featuresTile[pointTileFlatIndex] = node->_closestFeatures[(k*2 + j)*2 + i];
                validFlagsTile[pointTileFlatIndex] = true;
            }
        }
    }
}
                   
void AdaptiveDistanceField::recursiveSubdivToMaxLevel(const QueueNode& queueNode, const Vector3u& cellTileBmin, uint32 cellTileWidth,
                                                      Feature* featuresTile, bool* validFlagsTile, const QueryHint& queryHint)
{
    // make sure cellTileWidth is power of 2 and at least 2
    // make sure cellTileBmin and bbox are reasonable
    assert(cellTileWidth >= 2 && (cellTileWidth&(cellTileWidth-1)) == 0);
    for (int i = 0; i < 3; i++)
        assert(cellTileBmin[i] % cellTileWidth == 0);
    
    if (queueNode._level >= _maxAdfLevel)
        return;
    
    
    Vector3d halfCellSize = 0.5 * queueNode._cellSize;
    uint32 halfCellTileWidth = cellTileWidth / 2;
    
    // compute and cache closest features at cell, face, and edge centers and compute max error
    static const Vector3u offsets[19] = { Vector3u(0,0,1),Vector3u(0,1,0),Vector3u(0,1,1),Vector3u(0,1,2),Vector3u(0,2,1),
                                          Vector3u(1,0,0),Vector3u(1,0,1),Vector3u(1,0,2),Vector3u(1,1,0),Vector3u(1,1,1),
                                          Vector3u(1,1,2),Vector3u(1,2,0),Vector3u(1,2,1),Vector3u(1,2,2),
                                          Vector3u(2,0,1),Vector3u(2,1,0),Vector3u(2,1,1),Vector3u(2,1,2),Vector3u(2,2,1) };
    double nodeError = 0.0;
    for (size_t i = 0; i < 19; i++)
    {
        Vector3d point(queueNode._cellBmin.x + offsets[i].x * halfCellSize.x,
                       queueNode._cellBmin.y + offsets[i].y * halfCellSize.y,
                       queueNode._cellBmin.z + offsets[i].z * halfCellSize.z);
        Vector3u pointTileIndices = cellTileBmin + halfCellTileWidth * offsets[i];
        uint32 pointTileFlatIndex = (pointTileIndices[2]*_tileDim + pointTileIndices[1])*_tileDim + pointTileIndices[0];
        
        Feature* pointClosestFeature = &featuresTile[pointTileFlatIndex];
        if (!validFlagsTile[pointTileFlatIndex])
        {
            _closestFunc(point, pointClosestFeature, true, queryHint._point, queryHint._closest);
            validFlagsTile[pointTileFlatIndex] = true;
        }
        
        double error = std::abs(pointClosestFeature->_distance - queueNode.distance(point));
        if (error > nodeError)
            nodeError = error;
    }

    if (nodeError < _maxError && queueNode._level > 0)  // tree root should always subdivide regardless of error
        return;
    
    // if this cell doesn't come within _subdivideRadius of the surface, do not subdivide.
    Vector3u cellCenterTileIndices = cellTileBmin + cellTileWidth / 2;
    uint32 cellCenterTileFlatIndex = (cellCenterTileIndices[2]*_tileDim + cellCenterTileIndices[1])*_tileDim + cellCenterTileIndices[0];
    assert(validFlagsTile[cellCenterTileFlatIndex]);
    Feature* cellCenterClosestFeature = &featuresTile[cellCenterTileFlatIndex];
    double cellCenterDistance = std::abs(cellCenterClosestFeature->_distance);
    if (cellCenterDistance > _subdivideRadius + 0.5*(queueNode._cellSize.length()))
        return;
    
    // use this cell's center as the hint for all children cell queries
    QueryHint childQueryHint(queueNode._cellBmin + halfCellSize, *cellCenterClosestFeature);
           
    // subdivide cell
    BuildNode* node = queueNode._node;
    node->_children = new BuildNode[8];
    
    //bool tileMaxLevelReached = ((queueNode._level + 1) % _tileLevels == 0);
    bool tileMaxLevelReached = ((_maxAdfLevel - queueNode._level) % _tileLevels == 0);
    
    for (uint32 k = 0; k < 2; k++)
    {
        for (uint32 j = 0; j < 2; j++)
        {
            for (uint32 i = 0; i < 2; i++)
            {
                BuildNode* childNode = &(node->_children[(k*2 + j)*2 + i]);
                
                Vector3u childCellTileBmin(cellTileBmin.x + i * halfCellTileWidth,
                                           cellTileBmin.y + j * halfCellTileWidth,
                                           cellTileBmin.z + k * halfCellTileWidth);
                Vector3d childCellBmin(queueNode._cellBmin.x + i * halfCellSize.x,
                                       queueNode._cellBmin.y + j * halfCellSize.y,
                                       queueNode._cellBmin.z + k * halfCellSize.z);
                                       
                setNodeClosestFeaturesFromTile(childNode, featuresTile, validFlagsTile, childCellTileBmin, halfCellTileWidth);
                                       
                QueueNode childQueueNode(childNode, queueNode._level+1, childCellBmin, halfCellSize);
                if (tileMaxLevelReached)
                {
                    _nodeQueueLock.lock();
                    _nodeQueue.push(childQueueNode);
                    _nodeQueueCv.notify_all();
                    _nodeQueueLock.unlock();
                }
                else
                {
                    recursiveSubdivToMaxLevel(childQueueNode, childCellTileBmin, halfCellTileWidth, featuresTile, validFlagsTile, childQueryHint);   
                }
            }
        }
    }
    _numNodes += 8;
    if (_numNodes / 10000 != (_numNodes-8) / 10000) {
        std::cout << '\r' << _numNodes << std::flush;
    }
}

AdaptiveDistanceField::FeatureBlock::FeatureBlock(const BuildNode* children)
{
    assert(children != nullptr);
    static const int childAndFeatureIndices[27][2] = {{0,0}, {0,1}, {1,1},
                                                      {0,2}, {0,3}, {1,3},
                                                      {2,2}, {2,3}, {3,3},
                                                      
                                                      {0,4}, {0,5}, {1,5},
                                                      {0,6}, {0,7}, {1,7},
                                                      {2,6}, {2,7}, {3,7},
                                                      
                                                      {4,4}, {4,5}, {5,5},
                                                      {4,6}, {4,7}, {5,7},
                                                      {6,6}, {6,7}, {7,7} };
     for (int i = 0; i < 27; i++)
     {
         int childId = childAndFeatureIndices[i][0];
         int featureIndex = childAndFeatureIndices[i][1];
         _features[i] = children[childId]._closestFeatures[featureIndex];
     }
}

double AdaptiveDistanceField::FeatureBlock::distance(int childId, const Vector3d& lerpWeights) const
{
    static const int magicNumbers[8] = {0,1,3,4,9,10,12,13};
    
    int base = magicNumbers[childId];
    return TRILINEAR_INTERPOLATION(lerpWeights.x, lerpWeights.y, lerpWeights.z,
                                   _features[base+magicNumbers[0]]._distance,
                                   _features[base+magicNumbers[1]]._distance,
                                   _features[base+magicNumbers[3]]._distance,
                                   _features[base+magicNumbers[2]]._distance,
                                   _features[base+magicNumbers[4]]._distance,
                                   _features[base+magicNumbers[5]]._distance,
                                   _features[base+magicNumbers[7]]._distance,
                                   _features[base+magicNumbers[6]]._distance);
}

Vector3d AdaptiveDistanceField::FeatureBlock::gradient(int childId, const Vector3d& lerpWeights, const Vector3d& cellSize) const
{
    static const int magicNumbers[8] = {0,1,3,4,9,10,12,13};
    
    int base = magicNumbers[childId];
    return Vector3d(
                GRADIENT_COMPONENT_X(lerpWeights.x, lerpWeights.y, lerpWeights.z,
                                     _features[base+magicNumbers[0]]._distance,
                                     _features[base+magicNumbers[1]]._distance,
                                     _features[base+magicNumbers[3]]._distance,
                                     _features[base+magicNumbers[2]]._distance,
                                     _features[base+magicNumbers[4]]._distance,
                                     _features[base+magicNumbers[5]]._distance,
                                     _features[base+magicNumbers[7]]._distance,
                                     _features[base+magicNumbers[6]]._distance,
                                     cellSize.x),
                GRADIENT_COMPONENT_Y(lerpWeights.x, lerpWeights.y, lerpWeights.z,
                                     _features[base+magicNumbers[0]]._distance,
                                     _features[base+magicNumbers[1]]._distance,
                                     _features[base+magicNumbers[3]]._distance,
                                     _features[base+magicNumbers[2]]._distance,
                                     _features[base+magicNumbers[4]]._distance,
                                     _features[base+magicNumbers[5]]._distance,
                                     _features[base+magicNumbers[7]]._distance,
                                     _features[base+magicNumbers[6]]._distance,
                                     cellSize.y),
                GRADIENT_COMPONENT_Z(lerpWeights.x, lerpWeights.y, lerpWeights.z,
                                     _features[base+magicNumbers[0]]._distance,
                                     _features[base+magicNumbers[1]]._distance,
                                     _features[base+magicNumbers[3]]._distance,
                                     _features[base+magicNumbers[2]]._distance,
                                     _features[base+magicNumbers[4]]._distance,
                                     _features[base+magicNumbers[5]]._distance,
                                     _features[base+magicNumbers[7]]._distance,
                                     _features[base+magicNumbers[6]]._distance,
                                     cellSize.z)
                );
}



void AdaptiveDistanceField::AdfNode::setFeatureIndex(int featureBlockIndex, int childId)
{
    assert(0 <= childId && childId < 8);
    _index = -((featureBlockIndex << 3) | (childId));
}

void AdaptiveDistanceField::AdfNode::setChildrenIndex(int childrenIndex)
{
    _index = childrenIndex;
}

bool AdaptiveDistanceField::AdfNode::isLeaf() const
{
    return (_index <= 0);
}

int AdaptiveDistanceField::AdfNode::getChildIndex(int childId) const
{
    assert(_index > 0);
    return _index + childId;
}

int AdaptiveDistanceField::AdfNode::getFeatureBlockIndex(int* childId) const
{
    assert(_index <= 0);
    int bitfield = -_index;
    *childId = bitfield & 7;        // bits [0,2]
    return bitfield >> 3;           // bits 3 and above
}


void AdaptiveDistanceField::coalesce(const BuildNode* root)
{        
    _adfNodes.reserve(_numNodes);
    
    // coalesce root node
    assert(!root->isLeaf());
    _adfNodes.emplace_back();
    
    // recursively coalesce tree 
    recursiveCoalesceToMaxLevel(&_adfNodes[0], root->_children);
    
    assert((int)_adfNodes.size() == _numNodes);
    assert(recursiveNodeCount(&_adfNodes[0]) == _numNodes);
}

void AdaptiveDistanceField::recursiveCoalesceToMaxLevel(AdfNode* coalescedNode, BuildNode* children)
{
    assert(children != nullptr);
    
    coalescedNode->setChildrenIndex(_adfNodes.size());
    for (int childId = 0; childId < 8; childId++)
        _adfNodes.emplace_back();
    
    bool blockAdded = false;
    //const FeatureBlock* block = nullptr;
    int blockIndex = -1;
    
    for (int childId = 0; childId < 8; childId++)
    {
        BuildNode* child = &children[childId];
        AdfNode* coalescedChildNode = &_adfNodes[coalescedNode->getChildIndex(childId)];
        
        if (child->isLeaf())
        {
            if (!blockAdded)
            {
                blockIndex = _featureBlocks.size();
                _featureBlocks.emplace_back(children);
                //block = &_featureBlocks.back();
                blockAdded = true;
            }
            coalescedChildNode->setFeatureIndex(blockIndex, childId);
        }
        else
        {
            recursiveCoalesceToMaxLevel(coalescedChildNode, child->_children);
        }
    }
    
    delete[] children;
    children = nullptr;
}

int AdaptiveDistanceField::recursiveNodeCount(const AdfNode* node)
{
    //if (!node->_children)
    if (node->isLeaf())
        return 1;
    int sum = 1;
    for (int childId = 0; childId < 8; childId++)
    {
        sum += recursiveNodeCount(&_adfNodes[node->getChildIndex(childId)]);//&(node->_children[i]));
    }
    return sum;
}


const AdaptiveDistanceField::AdfNode* AdaptiveDistanceField::containingNode(const Vector3d& point, Vector3d* lerpWeights, Vector3d* cellSize) const
{
    // find the leaf node containing this point
    const AdfNode* node = &_adfNodes[0];
    Vector3d nodeCellBmin = _rootBmin;
    Vector3d nodeCellSize = _rootBmax - _rootBmin;
    while (!node->isLeaf())
    {
        // make sure point falls within the root bbox
        //for (int i = 0; i < 3; i++)
        //    assert(nodeCellBmin[i] <= point[i] && point[i] <= nodeCellBmin[i]+nodeCellSize[i]);
        
        Vector3d cellCenter = nodeCellBmin + 0.5 * nodeCellSize;
        int childId = 0;
        if (point.x > cellCenter.x)
        {
            nodeCellBmin.x = cellCenter.x;
            childId += 1;
        }
        if (point.y > cellCenter.y)
        {
            nodeCellBmin.y = cellCenter.y;
            childId += 2;
        }
        if (point.z > cellCenter.z)
        {
            nodeCellBmin.z = cellCenter.z;
            childId += 4;
        }
        nodeCellSize = 0.5 * nodeCellSize;
        node = &_adfNodes[node->getChildIndex(childId)];  //&(node->_children[childId]);
    }
    *cellSize = nodeCellSize;
    lerpWeights->x = (point.x - nodeCellBmin.x) / nodeCellSize.x;
    lerpWeights->y = (point.y - nodeCellBmin.y) / nodeCellSize.y;
    lerpWeights->z = (point.z - nodeCellBmin.z) / nodeCellSize.z;
    return node;
}

double AdaptiveDistanceField::distance(const Vector3d& point) const
{
    if (!(_rootBmin.x<=point.x && point.x<=_rootBmax.x && _rootBmin.y<=point.y && point.y<=_rootBmax.y && _rootBmin.z<=point.z && point.z<=_rootBmax.z))
        return DBL_MAX;
        
    Vector3d lerpWeights, cellSize;
    const AdfNode* node = containingNode(point, &lerpWeights, &cellSize);
    int childId;
    int featureBlockIndex = node->getFeatureBlockIndex(&childId);
    return _featureBlocks[featureBlockIndex].distance(childId, lerpWeights);
}

Vector3d AdaptiveDistanceField::gradient(const Vector3d& point) const
{
    if (!(_rootBmin.x<=point.x && point.x<=_rootBmax.x && _rootBmin.y<=point.y && point.y<=_rootBmax.y && _rootBmin.z<=point.z && point.z<=_rootBmax.z))
        return Vector3d(0.0, 0.0, 0.0);
        
    Vector3d lerpWeights, cellSize;
    const AdfNode* node = containingNode(point, &lerpWeights, &cellSize);
    int childId;
    int featureBlockIndex = node->getFeatureBlockIndex(&childId);
    return _featureBlocks[featureBlockIndex].gradient(childId, lerpWeights, cellSize);
}


void AdaptiveDistanceField::FeatureBlock::save(std::ostream& os) const
{
    os.write(reinterpret_cast<const char*>(_features), 27*sizeof(Feature));
}

void AdaptiveDistanceField::FeatureBlock::load(std::istream& is)
{
    is.read(reinterpret_cast<char*>(_features), 27*sizeof(Feature));
}


void AdaptiveDistanceField::AdfNode::save(std::ostream& os) const
{
    os.write(reinterpret_cast<const char*>(&_index), sizeof(int));
}

void AdaptiveDistanceField::AdfNode::load(std::istream& is)
{
    is.read(reinterpret_cast<char*>(&_index), sizeof(int));
}



int AdaptiveDistanceField::save(const std::string& filename) const
{
    std::ofstream fout(filename.c_str(), ios::out | ios::trunc | ios::binary);
    if (!fout)
        return 1;
    
    fout.write(reinterpret_cast<const char*>(&_subdivideRadius), sizeof(double));
    fout.write(reinterpret_cast<const char*>(&_maxAdfLevel), sizeof(int));
    fout.write(reinterpret_cast<const char*>(&_maxError), sizeof(double));
    
    int numNodes = _adfNodes.size();
    fout.write(reinterpret_cast<const char*>(&numNodes), sizeof(int));
    fout.write(reinterpret_cast<const char*>(&_rootBmin), sizeof(Vector3d));
    fout.write(reinterpret_cast<const char*>(&_rootBmax), sizeof(Vector3d));
    
    // write feature blocks
    int numFeatureBlocks = _featureBlocks.size();
    fout.write(reinterpret_cast<const char*>(&numFeatureBlocks), sizeof(int));
    for (const FeatureBlock& block : _featureBlocks)
        block.save(fout);
    
    // write adf nodes
    for (const AdfNode& node : _adfNodes)
        node.save(fout);
    
    return 0;
}


int AdaptiveDistanceField::load(const std::string& filename)
{
    std::ifstream fin(filename.c_str(), ios::in | ios::binary);
    if (!fin)
        return 1;
 
    _tileLevels = -1;
    _tileDim = -1;
    _closestFunc = ClosestFeatureFunc();
    
    fin.read(reinterpret_cast<char*>(&_subdivideRadius), sizeof(double));
    fin.read(reinterpret_cast<char*>(&_maxAdfLevel), sizeof(int));
    fin.read(reinterpret_cast<char*>(&_maxError), sizeof(double));
        
    int numNodes;
    fin.read(reinterpret_cast<char*>(&numNodes), sizeof(int));
    fin.read(reinterpret_cast<char*>(&_rootBmin), sizeof(Vector3d));
    fin.read(reinterpret_cast<char*>(&_rootBmax), sizeof(Vector3d));
    
    _numNodes = numNodes;
    
    // read feature blocks
    int numFeatureBlocks;
    fin.read(reinterpret_cast<char*>(&numFeatureBlocks), sizeof(int));
    _featureBlocks.resize(numFeatureBlocks);
    for (int i = 0; i < numFeatureBlocks; i++)
        _featureBlocks[i].load(fin);

    // read adf nodes and verify validity of the pointer stored in each node
    _adfNodes.resize(numNodes);
    for (int i = 0; i < numNodes; i++)
    {
        AdfNode& node = _adfNodes[i];
        node.load(fin);
        if (node.isLeaf())
        {
            int childId;
            int featuresBlockIndex = node.getFeatureBlockIndex(&childId);
            if (featuresBlockIndex >= numFeatureBlocks)
            {
                printf("node %d has an invalid feature-block index of %d\n", i, featuresBlockIndex);
                return 3;
            }
        }
        else
        {
            int lastChildIndex = node.getChildIndex(7);
            if (lastChildIndex >= numNodes)
            {
                printf("node %d has an invalid first-child index of %d\n", i, node.getChildIndex(0));
                return 3;
            }
        }
    }
    
    if (!(recursiveNodeCount(&_adfNodes[0]) == numNodes))
    {
        printf("whole-tree traversal counted different number of nodes than numNodes!\n");
        return 2;
    }
    
#if PNG
    std::string xyz[3] = {"x", "y", "z"};
    for (int n =0; n < 3; n++)
        savePng(n, filename+"_"+xyz[n]+".png");
#endif

    return 0;
}

#if PNG
static const int IM_SIZE = 1024;
    
void AdaptiveDistanceField::recursiveColorDots(int dim, double sliceCenter, double sliceWidth,
                                               const AdfNode* node, const Vector3d& cellBmin, const Vector3d& cellSize,
                                               std::vector<unsigned char>& image) const
{
    static const Vector3u offsets[8] = {Vector3u(0,0,0),Vector3u(1,0,0),Vector3u(0,1,0),Vector3u(1,1,0),
                                        Vector3u(0,0,1),Vector3u(1,0,1),Vector3u(0,1,1),Vector3u(1,1,1)};
                                        
    for (int j = 0; j < 8; j++) {
        Vector3d point(cellBmin.x + offsets[j].x * cellSize.x,
                       cellBmin.y + offsets[j].y * cellSize.y,
                       cellBmin.z + offsets[j].z * cellSize.z);
        if (std::abs(point[dim] - sliceCenter) > 0.5*sliceWidth)
            continue;
            
        int indices[3];
        for (int n = 0; n < 3; n++) {
            indices[n] = std::min(IM_SIZE, (int)std::round((point[n]-_rootBmin[n])/(_rootBmax[n]-_rootBmin[n])*IM_SIZE));
        }
        int flatIndex = dim == 1 ?
                        (indices[2]*(IM_SIZE +1) + indices[0]) :                // xz
                        (indices[(dim+2)%3]*(IM_SIZE+1) + indices[(dim+1)%3]);  // xy or yz
        image[4*flatIndex] = 255;
        image[4*flatIndex+1] = 0;
        image[4*flatIndex+2] = 0;
        image[4*flatIndex+3] = 255;
    }
    
    if (!node->isLeaf())
    {
        Vector3d halfCellSize = 0.5*cellSize;
        for (int j = 0; j < 8; j++) {
            Vector3d childCellBmin(cellBmin.x + offsets[j].x * halfCellSize.x,
                                   cellBmin.y + offsets[j].y * halfCellSize.y,
                                   cellBmin.z + offsets[j].z * halfCellSize.z);
            recursiveColorDots(dim, sliceCenter, sliceWidth, &_adfNodes[node->getChildIndex(j)], childCellBmin, halfCellSize, image);
        }
    }
}

void AdaptiveDistanceField::savePng(int dim, const std::string& filename) const
{
    {
        double sliceCenter = 0.5*(_rootBmax[dim] + _rootBmin[dim]);
        double sliceWidth = (_rootBmax[dim] - _rootBmin[dim]) / (double)(1<<_maxAdfLevel) * 0.5;
        
        std::vector<unsigned char> image((IM_SIZE+1)*(IM_SIZE+1) * 4, 0);
        recursiveColorDots(dim, sliceCenter, sliceWidth,
                           &_adfNodes[0], _rootBmin, _rootBmax - _rootBmin, image);
       
        Vector3d bbox = _rootBmax - _rootBmin;                    
        for (int j = 0; j < IM_SIZE+1; j++) {
            for (int i = 0; i < IM_SIZE+1; i++) {
                
                int flatIndex = j*(IM_SIZE+1) + i;
                if (image[4*flatIndex+3] != 0)
                    continue;
                
                Vector3d point;

                if (dim == 1)
                {
                    point[0] = (double)i/(IM_SIZE+1)*bbox[0] + _rootBmin[0];
                    point[2] = (double)j/(IM_SIZE+1)*bbox[2] + _rootBmin[2];
                }
                else
                {
                    point[(dim+1)%3] = (double)i/(IM_SIZE+1)*bbox[(dim+1)%3] + _rootBmin[(dim+1)%3];
                    point[(dim+2)%3] = (double)j/(IM_SIZE+1)*bbox[(dim+2)%3] + _rootBmin[(dim+2)%3];
                }
                point[dim] = sliceCenter;

                double dist = std::abs(distance(point));
                double mag = 1.0 - std::min(dist/_subdivideRadius, 1.0);
                Vector3d color (mag, mag, mag);
                
                /*Vector3d grad = gradient(point);
                grad.normalize();
                Vector3d color(std::abs(grad.x),
                               std::abs(grad.y),
                               std::abs(grad.z));
                */
                image[4*flatIndex] = 255 * color.x;
                image[4*flatIndex+1] = 255 * color.y;
                image[4*flatIndex+2] = 255 * color.z;
                image[4*flatIndex+3] = 255;
            }
        }
        
        lodepng::encode(filename, image, IM_SIZE+1, IM_SIZE+1);
    }
}
#endif
