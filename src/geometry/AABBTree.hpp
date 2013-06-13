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
 *       Filename:  AABBTree.hpp
 *
 *    Description:  Implement the AABB Tree
 *
 *        Version:  1.0
 *        Created:  03/10/11 18:33:55
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef AABB_TREE_INC
#   define AABB_TREE_INC

#include "AABBTreeNode.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

/*!
 * The AABB Tree
 */
template <typename T, class TMesh, class TData>
class AABBTree
{
    public:
        typedef TMesh                       TMesh;
        typedef TData                       TData;
        typedef AABBTree<T, TMesh, TData>   TSelf;
        typedef AABBTreeNode<T, TSelf>      TNode;

        static const int LEAF_SIZE;

    public:
        AABBTree(const TMesh* pmesh)
        {   init(pmesh); }

        ~AABBTree()
        {   delete root_; }

        TNode* root()
        {   return root_; }

        void init(const TMesh* pmesh);

    private:
        TNode*          root_;
};

// ----------------------------------------------------------------------------
template <typename T, class TMesh, class TData>
const int AABBTree<T, TMesh, TData>::LEAF_SIZE = 16;

template <typename T, class TMesh, class TData>
void AABBTree<T, TMesh, TData>::init(const TMesh* pmesh)
{
    if ( !pmesh )
    {
        root_ = NULL;
        return;
    }

    std::vector<int> tgls(pmesh->num_triangles());
    for(int i = 0;i < tgls.size();++ i) 
        tgls[i] = i;
    root_ = new TNode(pmesh, tgls, 0);
}

#ifdef USE_NAMESPACE
}
#endif
#endif
