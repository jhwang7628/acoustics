/*
 * =====================================================================================
 *
 *       Filename:  TrianglePatch.hpp
 *
 *    Description:  A subset of triangle mesh
 *
 *        Version:  1.0
 *        Created:  09/09/10 15:59:00
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  TRIANGLE_PATCH_HPP_INC
#define  TRIANGLE_PATCH_HPP_INC

#include <assert.h>
#include <set>
#include "geometry/TriangleMesh.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class TrianglePatch
{
    public:
        TrianglePatch(const TriangleMesh<T>* mesh):pmesh_(mesh)
        { }

        /* Add all triangles and vertices in the current mesh 
         * into the patch
         */
        void add_all() 
        {
            for(int i = 0;i < pmesh_->num_vertices();++ i) vs_.insert(i);
            for(int i = 0;i < pmesh_->num_triangles();++ i) ts_.insert(i);
        }

        void add_triangle(int tid) 
        {
            assert(tid >= 0 && tid < pmesh_->num_triangles());

            const Tuple3ui& ts = pmesh_->triangle_ids(tid);
            ts_.insert(tid);
            vs_.insert(ts.x);
            vs_.insert(ts.y);
            vs_.insert(ts.z);
        }

        /*
         * Generate an independent mesh based on this patch
         *   return the generated triangle mesh
         *   when it returns, vmap[i] is the i-th vertex ID in
         *   original mesh
         */
        TriangleMesh<T>* to_mesh(int* vmap) const;

        const TriangleMesh<T>* mesh() const
        {   return pmesh_; }

        int num_triangles() const
        {   return ts_.size(); }
        int num_verticies() const
        {   return vs_.size(); }
        const std::set<int>& vertices() const
        {   return vs_; }
        const std::set<int>& triangles() const
        {   return ts_; }

        /*!  
         * Check if the given triangle is in the mesh
         * \param trid  the id of the triangle on the original mesh
         */
        bool is_triangle_in(int trid) const
        {   return ts_.count(trid); }
        /*!  
         * Check if the given vertex is in the mesh
         * \param trid  the id of the vertex on the original mesh
         */
        bool is_vertex_in(int vid) const
        {   return vs_.count(vid); }

    private:
        const TriangleMesh<T>*  pmesh_;

        std::set<int>           vs_;    // vertex set
        std::set<int>           ts_;    // triangle set
};

// ----------------------------------------------------------------------------

template <typename T>
TriangleMesh<T>* TrianglePatch<T>::to_mesh(int* vmap) const
{
    using namespace std;
#ifdef USE_UNORDERED_MAP
    using namespace std::tr1;
    unordered_map< int, int > hash;
#else
    using namespace __gnu_cxx;
    hash_map< int, int > hash;
#endif

    TriangleMesh<T>* ret = new TriangleMesh<T>;

    /* Add vertices */
    const vector< Point3<T> >& vtx = pmesh_->vertices();
    set<int>::const_iterator end = vs_.end();
    for(set<int>::const_iterator it = vs_.begin();it != end;++ it)
    {
        int vid = ret->add_vertex(vtx[*it]);
        if ( vmap ) vmap[vid] = *it;
        hash[*it] = vid;
    }

    /* Add triangles */
    const vector< Tuple3ui >& tgl = pmesh_->triangles();
    end = ts_.end();
    for(set<int>::const_iterator it = ts_.begin();it != end;++ it)
    {
        assert(hash.count(tgl[*it].x) &&
               hash.count(tgl[*it].y) &&
               hash.count(tgl[*it].z)); 
        ret->add_triangle(hash[tgl[*it].x], hash[tgl[*it].y],
                          hash[tgl[*it].z]);
    }

    return ret;
}

#ifdef USE_NAMESPACE
}
#endif

#endif   /* ----- #ifndef TRIANGLE_PATCH_HPP_INC  ----- */

