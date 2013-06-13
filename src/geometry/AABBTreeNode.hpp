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
 *       Filename:  AABBTreeNode.hpp
 *
 *    Description:  class of AABB tree node
 *
 *        Version:  1.0
 *        Created:  03/10/11 18:40:02
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef AABB_TREE_NODE_INC
#   define AABB_TREE_NODE_INC

#include <set>
#include <limits>
#include <vector>
#include "geometry/Point3.hpp"
#include "utils/tuple.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T, class TTree>
class AABBTreeNode
{
    public:
        typedef typename TTree::TMesh   TMesh;
        typedef typename TTree::TData   TData;
        typedef AABBTreeNode<T, TTree>  TSelf;

    public:
        AABBTreeNode(const TMesh* pmesh, const std::vector<int>& tgls, 
                     int level);

        ~AABBTreeNode()
        {   
            delete child_[0];
            delete child_[1];
        }

        bool is_leaf() const
        {  return isLeaf_; }

        int level() const
        {   return level_; }
        int num_triangles() const
        {   return tgl_[0].size() + tgl_[1].size(); }

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

        size_t size() const
        {   return tgl_[0].size() + tgl_[1].size(); }

    private:
        int                 level_;     // tree level
        /* Bounding box of this node */
        Point3<T>           maxPt_;
        Point3<T>           minPt_;

        bool                isLeaf_;
        TSelf*              child_[2];
        std::vector<int>    tgl_[2];    // the contained triangles, given by indices
        std::set<int>       vtx_;       // the contained vertices, given by indices
                                        // only non-empty when this node is a leaf
        TData               data_;
};

// ----------------------------------------------------------------------------

template <typename T, class TTree>
AABBTreeNode<T, TTree>::AABBTreeNode(const TMesh* pmesh, 
        const std::vector<int>& tgls, int level):
        level_(level), 
        minPt_(std::numeric_limits<T>::infinity(),
               std::numeric_limits<T>::infinity(),
               std::numeric_limits<T>::infinity()),
        maxPt_(-std::numeric_limits<T>::infinity(),
               -std::numeric_limits<T>::infinity(),
               -std::numeric_limits<T>::infinity()),
        isLeaf_(false)
{
    /* find bounding box boundary */
    const std::vector<Tuple3ui>& tglInd = pmesh->surface_indices();
    const std::vector< Point3<T> >& vtx = pmesh->vertices();
    for(size_t i = 0;i < tgls.size();++ i)
    for(int j = 0;j < 3;++ j)
    {
        minPt_.x = fmin(vtx[tglInd[tgls[i]][j]].x, minPt_.x);
        minPt_.y = fmin(vtx[tglInd[tgls[i]][j]].y, minPt_.y);
        minPt_.z = fmin(vtx[tglInd[tgls[i]][j]].z, minPt_.z);

        maxPt_.x = fmax(vtx[tglInd[tgls[i]][j]].x, maxPt_.x);
        maxPt_.y = fmax(vtx[tglInd[tgls[i]][j]].y, maxPt_.y);
        maxPt_.z = fmax(vtx[tglInd[tgls[i]][j]].z, maxPt_.z);
    }

    /* Create child nodes */
    Vector3<T>  sz = maxPt_ - minPt_;
    const int    d = max_in_triple((T*)&sz);
    const T    mid = (maxPt_[d] + minPt_[d])*(T)0.5;
    for(size_t i = 0;i < tgls.size();++ i)
    {
        if ( (vtx[tglInd[tgls[i]][0]][d]<=mid) &&
             (vtx[tglInd[tgls[i]][1]][d]<=mid) &&
             (vtx[tglInd[tgls[i]][2]][d]<=mid) )
            tgl_[0].push_back(tgls[i]);
        else
            tgl_[1].push_back(tgls[i]);
    }

    if ( tgls.size() < TTree::LEAF_SIZE ||
         tgl_[0].empty() || tgl_[1].empty() )
    {
        isLeaf_ = true;
        child_[0] = child_[1] = NULL;
        /* cache all vertices */
        for(size_t i = 0;i < tgls.size();++ i)
        {
            vtx_.insert(tglInd[tgls[i]].x);
            vtx_.insert(tglInd[tgls[i]].y);
            vtx_.insert(tglInd[tgls[i]].z);
        }
    }
    else
    {
        isLeaf_ = false;
        child_[0] = new TSelf(pmesh, tgl_[0], level+1);
        child_[1] = new TSelf(pmesh, tgl_[1], level+1);
    }
}

#ifdef USE_NAMESPACE
}
#endif
#endif
