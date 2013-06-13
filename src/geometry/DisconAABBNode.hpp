/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  DisconAABBNode.hpp
 *
 *    Description:  The AABB Tree node with an indication about whether 
 *                  the triangles in this node are connected or not
 *
 *        Version:  1.0
 *        Created:  03/24/11 15:56:00
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef DISCON_AABB_NODE_INC
#   define DISCON_AABB_NODE_INC

#include <stack>
#include <limits>

#include "linearalgebra/Matrix3.hpp"
#include "geometry/TriangleMesh.hpp"
#include "linearalgebra/eig3.h"
#include "utils/tuple.hpp"
#include "utils/macros.h"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T, class TTree>
class DisconAABBNode
{
    public:
        typedef typename TTree::TMesh       TMesh;
        typedef typename TTree::TData       TData;
        typedef DisconAABBNode<T, TTree>    TSelf;

    public:
        DisconAABBNode(TTree* ptree, const TMesh* pmesh, TSelf* parent,
                       const std::set<int>& ts, int level, 
                       const std::vector< std::set<int> >& vtxTgls, 
                       bool check_conn);

        ~DisconAABBNode()
        {
            delete child_[0];
            delete child_[1];
        }

        bool is_connected_node() const
        {   return conn_; }
        bool is_leaf() const
        {  return isLeaf_; }
        int level() const
        {   return level_; }
        int num_triangles() const
        {   return tgl_[0].size() + tgl_[1].size(); }

        const TSelf* parent() const
        {   return parent_; }
        TSelf* parent()
        {   return parent_; }
        const TSelf* left_child() const
        {  return child_[0]; }
        TSelf* left_child() 
        {  return child_[0]; }
        const TSelf* right_child() const
        {  return child_[1]; }
        TSelf* right_child() 
        {  return child_[1]; }

        const std::vector<int>& triangles(int c) const
        {
            assert(c==0 || c==1);
            return tgl_[c];
        }

        TData* data() {  return &data_; }
        const TData* data() const {  return &data_; }

        const Point3<T>& max_bbox_pt() const
        {   return maxPt_; }

        const Point3<T>& min_bbox_pt() const
        {   return minPt_; }

        size_t size() const
        {   return tgl_[0].size() + tgl_[1].size(); }

#ifdef BVH_PRECOMP_MAT
        const std::vector<int>& vertices() const
        {   return vtx_; }
#endif

    private:
        void eigen_decomp(const Matrix3<T>& CM);

        bool sep_disconn_node(TTree* ptree, const TMesh* pmesh, 
                              const std::set<int>& ts, 
                              const std::vector< std::set<int> >& vtxTgls);

        void sep_children(TTree* ptree, const TMesh* pmesh,
                          const Point3<T>& meanPt, const std::set<int>& ts,
                          const std::vector< std::set<int> >& vtxTgls);

        void update_bounding_box(const TMesh* pmesh, const std::set<int>& ts);

    private:
        bool    conn_;      // whether or not all the triangles are connected
        int     level_;     // node level on the BVH

        /* Bounding box of this node */
        Point3<T>   maxPt_;
        Point3<T>   minPt_;

        bool                isLeaf_;
        TSelf*              parent_;
        TSelf*              child_[2];
        std::vector<int>    tgl_[2];    // the contained triangles, given by indices
#ifdef BVH_PRECOMP_MAT
        std::vector<int>    vtx_;       // the contained vertices, given by indices
#endif
        Vector3<T>          projVec_;
        TData               data_;
};

// ----------------------------------------------------------------------------

template <typename T, class TTree>
DisconAABBNode<T, TTree>::DisconAABBNode(TTree* ptree, const TMesh* pmesh, 
        TSelf* parent, const std::set<int>& ts, int level,
        const std::vector< std::set<int> >& vtxTgls, 
        bool check_conn):level_(level), parent_(parent)
{
#ifdef BVH_PRECOMP_MAT
    std::set<int>  vset;
    const std::vector<Tuple3ui>& indices = pmesh->surface_indices();
    const std::set<int>::const_iterator end = ts.end();
    for(std::set<int>::const_iterator it = ts.begin();it != end;++ it)
    {
        tgl_[0].push_back(*it);

        for(int j = 0;j < 3;++ j)
        {
            int cvid = indices[*it][j];
            if ( !vset.count(cvid) )
            {
                vset.insert(cvid);
                vtx_.push_back(cvid);
            }
        }
    }
#endif

    /* Compute bounding box */
    update_bounding_box(pmesh, ts);

    /* 
     * NOTE: the leaf node may not be connected. But we are not applying
     *       any certificate computation on it anyways
     */
    if ( (int)ts.size() <= TTree::LEAF_SIZE )   // Create leaf node, no children
    {
        conn_     = false;
        child_[0] = child_[1] = NULL;
        isLeaf_   = true;

        return;
    }

    /* ------ Create children ------ */
    if ( check_conn && sep_disconn_node(ptree, pmesh, ts, vtxTgls) ) return;
    
    conn_   = true;
    // for connected node, separate its triangles
    const std::vector<typename TTree::moment>& moments = ptree->triangle_moments();

    T           atot = 0;
    Point3<T>   wPtSum;
    Matrix3<T>  CM;

    /* eigen decomposition */
    const std::set<int>::const_iterator tsend = ts.end();
    for(std::set<int>::const_iterator it = ts.begin();it != tsend;++ it)
    {
        atot += moments[*it].area;
        wPtSum.scaleAdd(moments[*it].area, moments[*it].centroid);
        CM   += moments[*it].C;
    }

    const T invA = (T)1 / atot;
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
        CM.cols[i][j] -= wPtSum[j]*wPtSum[i]*invA;

    eigen_decomp(CM);   // find the projection direction
    wPtSum *= invA;     // the mean point

    /* create children node */
    sep_children(ptree, pmesh, wPtSum, ts, vtxTgls);
}

/*
 * - eigen decomposition
 */
template <typename T, class TTree>
void DisconAABBNode<T, TTree>::eigen_decomp(const Matrix3<T>& CM)
{
    double V[3][3], d[3];
    eigen_decomposition((const double(*)[3])&CM, V, d);

    int maxid = 0;
    if ( fabs(d[1]) > fabs(d[maxid]) ) maxid = 1;
    if ( fabs(d[2]) > fabs(d[maxid]) ) maxid = 2;

    projVec_.set(V[0][maxid], V[1][maxid], V[2][maxid]);
}

template <typename T, class TTree>
bool DisconAABBNode<T, TTree>::sep_disconn_node(TTree* ptree, 
        const TMesh* pmesh, const std::set<int>& ts, 
        const std::vector< std::set<int> >& vtxTgls)
{
    const std::vector<Tuple3ui>& tvids = pmesh->triangles();

    // check if this node is connected
    std::set<int>       ts1;    // first child, connected patch
    std::stack<int>     tstk;
    tstk.push(*(ts.begin()));   // BFS stack
    ts1.insert(*(ts.begin()));  // first child set

    while ( !tstk.empty() ) 
    {
        int tid = tstk.top();
        tstk.pop();

        for(int vid = 0;vid < 3;++ vid)
        {
            const int vtx_id = tvids[tid][vid];

            const std::set<int>::const_iterator end = vtxTgls[vtx_id].end();
            for(std::set<int>::const_iterator it = vtxTgls[vtx_id].begin();
                it != end;++ it)
                if ( !ts1.count(*it) && ts.count(*it) )
                {
                    ts1.insert(*it);
                    tstk.push(*it);
                }
        }
    }

    assert(ts1.size() <= ts.size());
    if ( ts1.size() < ts.size() ) 
    {
        std::set<int>   ts2;
        const std::set<int>::const_iterator end = ts.end();
        for(std::set<int>::const_iterator it = ts.begin();it != end;++ it)
            if ( !ts1.count(*it) ) 
            {
                ts2.insert(*it);    // ts2 = triangles not in ts1
                tgl_[1].push_back(*it);
            }

        const std::set<int>::const_iterator e1 = ts1.end();
        for(std::set<int>::const_iterator it = ts1.begin();it != e1;++ it)
            tgl_[0].push_back(*it); // triangles in ts1

        conn_     = false;
        isLeaf_   = false;
        child_[0] = new TSelf(ptree, pmesh, this, ts1, level_+1, vtxTgls, false);
        child_[1] = new TSelf(ptree, pmesh, this, ts2, level_+1, vtxTgls, true);

        return true;
    }

    return false;   // connected node
}

template <typename T, class TTree>
void DisconAABBNode<T, TTree>::update_bounding_box(
        const TMesh* pmesh, const std::set<int>& ts)
{
    maxPt_.set(-std::numeric_limits<T>::infinity(),
               -std::numeric_limits<T>::infinity(),
               -std::numeric_limits<T>::infinity());
    minPt_.set( std::numeric_limits<T>::infinity(),
                std::numeric_limits<T>::infinity(),
                std::numeric_limits<T>::infinity());

    const std::vector<Tuple3ui>& tglInd = pmesh->surface_indices();
    const std::vector< Point3<T> >& vtx = pmesh->vertices();

    const std::set<int>::const_iterator end = ts.end();
    for(std::set<int>::const_iterator it = ts.begin();it != end;++ it)
    for(int i = 0;i < 3;++ i)
    {
        minPt_.x = fmin(vtx[tglInd[*it][i]].x, minPt_.x);
        minPt_.y = fmin(vtx[tglInd[*it][i]].y, minPt_.y);
        minPt_.z = fmin(vtx[tglInd[*it][i]].z, minPt_.z);

        maxPt_.x = fmax(vtx[tglInd[*it][i]].x, maxPt_.x);
        maxPt_.y = fmax(vtx[tglInd[*it][i]].y, maxPt_.y);
        maxPt_.z = fmax(vtx[tglInd[*it][i]].z, maxPt_.z);
    }
}

template <typename T, class TTree>
void DisconAABBNode<T, TTree>::sep_children(
        TTree* ptree, const TMesh* pmesh,
        const Point3<T>& meanPt, const std::set<int>& ts,
        const std::vector< std::set<int> >& vtxTgls)
{
    const std::vector<typename TTree::moment>& moments = ptree->triangle_moments();

    std::set<int> ts1, ts2;
    tgl_[0].clear();
    tgl_[1].clear();
    const T axdmp = projVec_.dotProduct(meanPt);

    const std::set<int>::const_iterator end = ts.end();
    for(std::set<int>::const_iterator it = ts.begin(); it != end;++ it)
    {
        T td = projVec_.dotProduct(moments[*it].centroid);
        if ( td <= axdmp )
        {
            tgl_[0].push_back(*it);
            ts1.insert(*it);
        }
        else
        {
            tgl_[1].push_back(*it);
            ts2.insert(*it);
        }
    }

    if ( ts1.empty() || ts2.empty() ) 
    {
        isLeaf_ = true;
        child_[0] = child_[1] = NULL;
    }
    else
    {
        isLeaf_ = false;
        child_[0] = new TSelf(ptree, pmesh, this, ts1, level_+1, vtxTgls, true);
        child_[1] = new TSelf(ptree, pmesh, this, ts2, level_+1, vtxTgls, true);
    }
}

#ifdef USE_NAMESPACE
}
#endif

#endif
