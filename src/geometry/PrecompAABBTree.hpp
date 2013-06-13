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
 *       Filename:  PrecompAABBTree.hpp
 *
 *    Description:  Precomputed AABB Tree. The BVH is precomputed and can be loaded
 *                  from files at runtime
 *
 *        Version:  1.0
 *        Created:  10/13/11 16:38:43
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef PRECOMP_AABB_TREE_INC
#   define PRECOMP_AABB_TREE_INC

#include <ifstream>
#include <cctype>
#include <cstring>
#include "utils/macros.h"
#include "utils/strings.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <class TNode, class TData>
class PrecompAABBTree
{
    public:
        typedef TData       TData;
        typedef TNode       TNode;

    public:
        /* Load the precomputed AABB tree from a file */
        PrecompAABBTree(const char* filename);

    private:
        TNode*                  root_;  // the root node
        std::vector<TNode*>     nodes_; // all the nodes on the tree
};

// ----------------------------------------------------------------------------

template <class TNode, class TData>
PrecompAABBTree::PrecompAABBTree(const char* filename)
{
    using namespace std;

    char text[1024];
    vector<int> intval(3);
    int zone = 0;
    ifstream fin(filename);
    if ( fin.fail() ) SHOULD_NEVER_HAPPEN(2);

    // --- process each line of the file ---
    fin.getline(text, 1024);
    while ( !fin.fail() ) 
    {
        // --- trim left ---
        int len = strlen(text);
        int tbegin = 0;
        for(;tbegin < len && isspace(text[tbegin]);++ tbegin);
        
        // --- empty line OR comments ---
        if ( tbegin >= len || text[tbegin] == '#' ) goto NEXT_LINE;

        if ( zone == 0 ) 
        {
        }

NEXT_LINE:
        fin.getline(text, 1024);
    }
    fin.close();
}

#ifdef USE_NAMESPACE
}
#endif
#endif
