#ifndef _ADAPTIVEDISTANCEFIELD_H_
#define _ADAPTIVEDISTANCEFIELD_H_

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>
#include <functional>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include "linearalgebra/Vector3.hpp"


#define PNG 1   // for testing

struct Feature
{
    double _distance;
    //Vector3d _position;
    //int _indices[3];                // indices of triangle vertices
    //double _alpha, _beta, _gamma;   // barycentric coordinates
};

class AdaptiveDistanceField
{
    typedef uint_fast32_t uint32;
    typedef Vector3<uint32> Vector3u;
    typedef std::function<void(const Vector3d& point, Feature* closest,
                               bool hintValid, const Vector3d& hintPoint, const Feature& hintClosest)>
            ClosestFeatureFunc;

public:

    
private:
    struct BuildNode
    {
        BuildNode() : _children(nullptr) {}
        
        double distance(const Vector3d& pt, const Vector3d& cellBmin, const Vector3d& cellSize) const;
        //Vector3d gradient(const Vector3d& pt, const Vector3d& cellBmin, const Vector3d& cellSize) const;
        
        bool isLeaf() const { return (_children == nullptr); }
        
        Feature _closestFeatures[8];    // zyx: 000,001,...111
        BuildNode* _children;           // zyx: 000,001,...111
    };
    
    struct QueueNode
    {
        QueueNode() : _node(nullptr), _level(-1){}
        QueueNode(BuildNode* node, int level, const Vector3d& cellBmin, const Vector3d& cellSize) :
            _node(node), _level(level), _cellBmin(cellBmin), _cellSize(cellSize) {}
            
        double distance(const Vector3d& pt) const;
            
        BuildNode* _node;
        int _level;
        Vector3d _cellBmin;
        Vector3d _cellSize;
    };
    
    struct QueryHint
    {
        QueryHint() {}
        QueryHint(const Vector3d& point, const Feature& closest) : _point(point), _closest(closest) {}
        
        Vector3d _point;
        Feature _closest;
    };
    
    struct FeatureBlock
    {
        FeatureBlock() {}
        FeatureBlock(const BuildNode* children);
        //void fillWithChildren(const BuildNode* children);
        
        double distance(int childId, const Vector3d& lerpWeights) const;
        Vector3d gradient(int childId, const Vector3d& lerpWeights, const Vector3d& cellSize) const;
        
        void save(std::ostream& os) const;
        void load(std::istream& is);
        
    private:
        Feature _features[27];    // zyx: 000,001,002,...222
    };
    
    struct AdfNode
    {
        AdfNode() : _index(0) {}
        void setFeatureIndex(int featureBlockIndex, int childId);
        void setChildrenIndex(int childrenIndex);
        
        bool isLeaf() const;
        int getChildIndex(int childId) const;
        int getFeatureBlockIndex(int* childId) const;
        
        void save(std::ostream& os) const;
        void load(std::istream& is);
    
    private:
        int _index;  // if > 0: _adfNodes[index] is the first child.  if <= 0: abs(_index) is bitfield of featureBlockIndex and childId
    };
    
public:
    AdaptiveDistanceField();
    
    void computeSignedField(ClosestFeatureFunc closestFunc,
                            const Vector3d& bmin, const Vector3d& bmax, double subdivideRadius, 
                            int maxAdfLevel, double maxError, int tileLevels=4,
                            int numWorkerThreads=std::thread::hardware_concurrency());
                           
    void computeSignedField(const std::string& filename, double subdivideRadius,
                            int maxAdfLevel, double maxError, int tileLevels=4,
                            // objfileOctree params
                            int maxTriCount=15, int maxDepth=10, bool allBoxSidesEqual=true, double expansionRatio=1.5,
                            int numWorkerThreads=std::thread::hardware_concurrency()
                            );
                           
    double distance(const Vector3d& point) const;
    Vector3d gradient(const Vector3d& point) const;
    
    int save(const std::string& filename) const;
    int load(const std::string& filename);
    
private:    
    // prevent copying
    AdaptiveDistanceField(const AdaptiveDistanceField&);
    AdaptiveDistanceField& operator=(const AdaptiveDistanceField&);
    
    // thread-pool worker function
    void nodeQueueThreadFunc();
    
    void computeNodeClosestFeatures(BuildNode* node, const Vector3d& cellBmin, const Vector3d& cellSize) const;
    void setNodeClosestFeaturesFromTile(BuildNode* node, Feature* featuresTile, bool* validFlagsTile,
                                        const Vector3u& cellTileBmin, uint32 cellTileWidth) const;    
    void setTileCornersFromNodeClosestFeatures(Feature* featuresTile, bool* validFlagsTile, BuildNode* node) const;
    
    // helpers, as in (Perry & Frisken, 2001)
    void recursiveSubdivToMaxLevel(const QueueNode& queueNode, const Vector3u& cellTileBmin, uint32 cellTileWidth,
                                   Feature* featuresTile, bool* validFlagsTile, const QueryHint& queryHint);
                                   
    void coalesce(const BuildNode* root);                               
    void recursiveCoalesceToMaxLevel(AdfNode* coalescedNode, BuildNode* children);
    
    int recursiveNodeCount(const AdfNode* node);
    
    const AdfNode* containingNode(const Vector3d& point, Vector3d* lerpWeights, Vector3d* cellSize) const;
    
        
    // for visualizing the adf
#if PNG
    void savePng(int dim, const std::string& filename) const;
    void recursiveColorDots(int dim, double sliceCenter, double sliceWidth,
                            const AdfNode* node, const Vector3d& cellBmin, const Vector3d& cellSize,
                            std::vector<unsigned char>& image) const;
#endif

private:    
    double _subdivideRadius;
    int _tileLevels;      // L
    int _maxAdfLevel;
    double _maxError;
    ClosestFeatureFunc _closestFunc;
    
    uint32 _tileDim;
    Vector3d _rootBmin, _rootBmax;
    
    std::queue<QueueNode> _nodeQueue;
    int _numNodesInProgress;
    std::mutex _nodeQueueLock;
    std::condition_variable _nodeQueueCv;
    
    std::atomic<int> _numNodes;
    
    std::vector<AdfNode> _adfNodes;
    std::vector<FeatureBlock> _featureBlocks;
};

#endif
