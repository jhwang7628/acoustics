#ifndef RIGID_POINT_CLOUD_CONSTRAINT_HPP
#   define RIGID_POINT_CLOUD_CONSTRAINT_HPP

#include <vector>
#include <fstream>
#include "utils/print_msg.h"
#include "generic/null_type.hpp"
#include "geometry/PCOBBTree.hpp"
#include "LSCollisionDetect.hpp"
#include "CollisionConstraint.hpp"

/*
 * The point cloud is fixed in the space
 */
template <typename T, class TCollProc>
class PointCloudConstraint : public CollisionConstraint<T, TCollProc>
{
    public:
        typedef typename TCollProc::TTreeNode       TTreeNode;
        typedef PCOBBTree<T, carbine::NullType>     TPCTree;
        typedef PCOBBTreeNode<T, TPCTree>           TPCTreeNode;

        PointCloudConstraint():mu_(0), eps_(0), pcTree_(NULL) { }
        PointCloudConstraint(T mu, T eps):mu_(mu), eps_(eps), pcTree_(NULL) { }

        ~PointCloudConstraint()
        {   delete pcTree_; }

        void load_points_from_file(const char* file);

        bool is_fixed() const { return true; }
        const std::vector< Point3<T> >& points() const
        {   return pts_; }

        bool deepest_collision(TCollProc* cp, TTreeNode* node, CollisionRec<T>& outRec)
        {   return deepest_collision(cp, node, pcTree_->root(), outRec); }
        
        bool is_colliding(TCollProc* cp, TTreeNode* node)
        {   return is_colliding(cp, node, pcTree_->root()); }

    private:
        bool is_colliding(TCollProc* cp, TTreeNode* node, const TPCTreeNode* pcnode);
        bool deepest_collision(TCollProc* cp, TTreeNode* node, 
                const TPCTreeNode* pcnode, CollisionRec<T>& outRec);
        bool is_disjoint(TCollProc* cp, TTreeNode* node, const TPCTreeNode* pcnode);

    private:
        T           mu_;       // friction coefficient
        T           eps_;      // restitution coefficient
        std::vector< Point3<T> >    pts_;   // point cloud
        TPCTree*    pcTree_;
};

// --------------------------------------------------------------------------------------
template <typename T, class TCollProc>
void PointCloudConstraint<T, TCollProc>::load_points_from_file(const char* file)
{
    using namespace std;

    ifstream fin(file);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot open file: %s to read\n", file);
        exit(1);
    }

    int N;
    fin >> N;
    pts_.resize(N);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file %s correctly\n", file);
        exit(1);
    }
    for(int i = 0;i < N;++ i)
        fin >> pts_[i].x >> pts_[i].y >> pts_[i].z;
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file %s correctly\n", file);
        exit(1);
    }
    fin.close();

    printf("INFO: %d points loaded\n", N);
    pcTree_  = new TPCTree(&pts_[0], N);
}

template <typename T, class TCollProc>
bool PointCloudConstraint<T, TCollProc>::is_disjoint(TCollProc* cp, 
        TTreeNode* node, const TPCTreeNode* pcnode)
{
    TreeNodeData<T>* ndData = node->data();

    // make sure predc0, predR is updated
    if ( ndData->ts < cp->pred_timestamp() ) cp->update_tree_node_state(ndData, node);

    const Vector3<T> dirab = ndData->predc0 - pcnode->c();
    const Tuple3<T>& ra = node->r();
    const Tuple3<T>& rb = pcnode->r();
    const Matrix3<T>& axisB = pcnode->R();

    // the three principle dir of A (OBB tree node) as the projection dir
    for(int i = 0;i < 3;++ i)
    {
        T tl = fabs(dirab.dotProduct(ndData->predR.cols[i]));
        T tr = fabs(axisB.cols[0].dotProduct(ndData->predR.cols[i])*rb[0]) +
               fabs(axisB.cols[1].dotProduct(ndData->predR.cols[i])*rb[1]) + 
               fabs(axisB.cols[2].dotProduct(ndData->predR.cols[i])*rb[2]) +
               ra[i];
        if ( tl > tr ) return true;
    }
    
    // the three principle dir of B (PCOBB tree node) as the projection dir
    for(int i = 0;i < 3;++ i)
    {
        T tl = fabs(dirab.dotProduct(axisB.cols[i]));
        T tr = fabs(ndData->predR.cols[0].dotProduct(axisB.cols[i])*ra[0]) +
               fabs(ndData->predR.cols[1].dotProduct(axisB.cols[i])*ra[1]) +
               fabs(ndData->predR.cols[2].dotProduct(axisB.cols[i])*ra[2]) +
               rb[i];
        if ( tl > tr ) return true;
    }

    double R[3][3];
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
        R[i][j] = fabs(ndData->predR.cols[i].dotProduct(axisB.cols[j]));

    // cross product of the principle dir A and B as the projection dir
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
    {
        // pd = A_i x B_j
        Vector3<T> pd = ndData->predR.cols[i].crossProduct(axisB.cols[j]);
        if ( pd.lengthSqr() < 1E-18 ) continue;

        T tl = fabs(dirab.dotProduct(pd));
        T tr = 0;
        for(int k = 0;k < 3;++ k)
        {
            if ( k != i ) tr += ra[k]*R[3-i-k][j];
            if ( k != j ) tr += rb[k]*R[i][3-j-k];
        }
        if ( tl > tr ) return true;
    }
    return false;
}

template <typename T, class TCollProc>
bool PointCloudConstraint<T, TCollProc>::deepest_collision(
        TCollProc* cp, TTreeNode* node, const TPCTreeNode* pcnode, 
        CollisionRec<T>& outRec)
{
    //// check if the bounding box in prediction state is overlap with each other
    if ( is_disjoint(cp, node, pcnode) ) return false;

    //// if both bounding tree nodes are leaf nodes,
    //   do the level-set based collision detection
    if ( node->is_leaf() && pcnode->is_leaf() )
    {
        bool ret = false;
        Vector3<T> nml;
        const std::vector<int>& ptIds = pcTree_->point_indices();
        for(int i = pcnode->start_id();i < pcnode->end_id();++ i)
        {
            // transform the collision point to the object's initial state
            const Point3<T> pt = cp->rigid_body()->initial_predicted_position(pts_[ptIds[i]]);
            const T isoval = cp->levelset()->negative_dist_with_normal(pt, nml);    // nml is already normalized
            if ( isoval < 0. && isoval < outRec.depth )
            {
                // transform the normal into body's current configuration
                const Vector3<T> preNml = cp->rigid_body()->predicted_normal(nml);
                // check if the vertex and the object have nonseparating relative velocity
                // now the normal is point outward to the body
                const Vector3<T> vab = cp->rigid_body()->predicted_velocity(pts_[ptIds[i]]); // given input point should be in current/predicted configuration
                const T vnrel = vab.dotProduct(preNml);

                if ( vnrel > 0. )
                {
                    outRec.depth        = isoval;
                    outRec.pt           = pts_[ptIds[i]];
                    outRec.impulseDir   = -preNml;
                    outRec.vnrel        = -vnrel;
                    outRec.vrel         = vab;
                    outRec.eps          = fmin(eps_, cp->rigid_body()->rest_coeff());
                    outRec.mu           = fmax(mu_,  cp->rigid_body()->friction_coeff());
                    ret = true;
                }
            }
        }
        return ret;
    }

    bool b1, b2;
    if ( !node->is_leaf() )
    {
        // both A and B are not leaf
        if ( !pcnode->is_leaf() && node->size() < pcnode->size() )
        {
            b1 = deepest_collision(cp, node, pcnode->left_child(), outRec);
            b2 = deepest_collision(cp, node, pcnode->right_child(), outRec);
        } 
        else
        {
            b1 = deepest_collision(cp, node->left_child(), pcnode, outRec);
            b2 = deepest_collision(cp, node->right_child(), pcnode, outRec);
        }
    }
    else    // node is a leaf OBB tree node, now pcnode[Id] should not be a leaf
    {
        b1 = deepest_collision(cp, node, pcnode->left_child(), outRec);
        b2 = deepest_collision(cp, node, pcnode->right_child(), outRec);
    }
    return b1 || b2;
}

template <typename T, class TCollProc>
bool PointCloudConstraint<T, TCollProc>::is_colliding(
        TCollProc* cp, TTreeNode* node, const TPCTreeNode* pcnode)
{
    //// check if the bounding box in prediction state is overlap with each other
    if ( is_disjoint(cp, node, pcnode) ) return false;

    //// if both bounding tree nodes are leaf nodes,
    //   do the level-set based collision detection
    if ( node->is_leaf() && pcnode->is_leaf() )
    {
        const std::vector<int>& ptIds = pcTree_->point_indices();
        for(int i = pcnode->start_id();i < pcnode->end_id();++ i)
        {
            // transform the collision point to the object's initial state
            const Point3<T> pt = cp->rigid_body()->initial_predicted_position(pts_[ptIds[i]]);
            const T isoval = cp->levelset()->distance(pt);
            if ( isoval < 0. ) return true;
        }
        return false;
    }

    if ( !node->is_leaf() )
    {
        // both A and B are not leaf
        if ( !pcnode->is_leaf() && node->size() < pcnode->size() )
            return is_colliding(cp, node, pcnode->left_child()) ||
                   is_colliding(cp, node, pcnode->right_child());
        else
            return is_colliding(cp, node->left_child(), pcnode) ||
                   is_colliding(cp, node->right_child(), pcnode);
    }
    else    // node is a leaf OBB tree node, now pcnode[Id] should not be a leaf
        return is_colliding(cp, node, pcnode->left_child()) ||
               is_colliding(cp, node, pcnode->right_child());
}

#endif

