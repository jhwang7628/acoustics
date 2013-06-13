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
 *       Filename:  DisconAABBTree.hpp
 *
 *    Description:  The AABB Tree with indications about whether 
 *                  the triangles in tree nodes are connected or not
 *
 *        Version:  1.0
 *        Created:  03/25/11 22:47:30
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef DISCON_AABB_TREE_INC
#   define DISCON_AABB_TREE_INC

#include <vector>
#include <set>
#include "DisconAABBNode.hpp"
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T, class TMesh, class TData>
class DisconAABBTree
{
    public:
        typedef TMesh                           TMesh;
        typedef TData                           TData;
        typedef DisconAABBTree<T, TMesh, TData> TSelf;
        typedef DisconAABBNode<T, TSelf>        TNode;

        static const int LEAF_SIZE;

        struct moment
        {
            T           area;
            Point3<T>   centroid;
            Matrix3<T>  C;
        };

    public:
        DisconAABBTree(const TMesh* pmesh)
        {   init(pmesh); }

        ~DisconAABBTree()
        {   delete root_; }

        const std::vector<moment>& triangle_moments() const
        {   return moments_; }

        void init(const TMesh* pmesh);

        TNode* root()
        {   return root_; }

    private:
        TNode*              root_;
        std::vector<moment> moments_;
};

// ----------------------------------------------------------------------------
template <typename T, class TMesh, class TData>
const int DisconAABBTree<T, TMesh, TData>::LEAF_SIZE = 16;

template <typename T, class TMesh, class TData>
void DisconAABBTree<T, TMesh, TData>::init(const TMesh* pmesh)
{
    if ( !pmesh )
    {
        root_ = NULL;
        return;
    }

    const std::vector< Tuple3ui >&  indices = pmesh->surface_indices();
    const std::vector< Point3<T> >& vtx     = pmesh->vertices();
    std::set<int> ts;
    moments_.resize(indices.size());

    /* compute moments */
    const T TS  = 1./3.;
    const T TS2 = 1./12.;
    for(size_t tid = 0;tid < indices.size();++ tid)
    {
        ts.insert((int)tid);

        const Point3<T>& v0 = vtx[indices[tid][0]]; 
        const Point3<T>& v1 = vtx[indices[tid][1]]; 
        const Point3<T>& v2 = vtx[indices[tid][2]];

        moments_[tid].area = Triangle<T>::area(v0, v1, v2);
        moments_[tid].centroid = (v0 + v1 + v2) * TS;

        if ( moments_[tid].area < 1E-18 ) 
        {
            for(int k = 0;k < 3;++ k)
            {
                moments_[tid].C.cols[k][k] =
                    M_SQR(v0[k]) + M_SQR(v1[k]) + M_SQR(v2[k]);
                for(int j = k+1;j < 3;++ j)
                    moments_[tid].C.cols[k][j] = moments_[tid].C.cols[j][k] =
                        v0[k]*v0[j] + v1[k]*v1[j] + v2[k]*v2[j];
            }
        }
        else
        {
            for(int k = 0;k < 3;++ k)
            {
                moments_[tid].C.cols[k][k] = moments_[tid].area * TS2 *
                    (9.*M_SQR(moments_[tid].centroid[k]) + 
                     M_SQR(v0[k]) + M_SQR(v1[k]) + M_SQR(v2[k]));
                for(int j = k+1;j < 3;++ j)
                {
                    moments_[tid].C.cols[k][j] = moments_[tid].C.cols[j][k] = 
                        moments_[tid].area * TS2 *
                        (9.*moments_[tid].centroid[k]*moments_[tid].centroid[j] +
                         v0[k]*v0[j] + v1[k]*v1[j] + v2[k]*v2[j]);
                }
            }
        }
    } //// end for

    std::vector< std::set<int> > vtxTgls;
    pmesh->get_vtx_tgls(vtxTgls);
    root_ = new DisconAABBNode<T, TSelf>(this, pmesh, NULL, ts, 0, vtxTgls, false);
}

#ifdef USE_NAMESPACE
}
#endif

#endif
