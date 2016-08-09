#ifndef DIFF_DEFINE
/******************************************************************************
 *  File: TriangleMesh.hpp
 *
 *  This file is part of isostuffer
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
 ******************************************************************************/
#endif /* ! DIFF_DEFINE */
#ifndef GEOMETRY_TRIANGLE_MESH_HPP
#   define GEOMETRY_TRIANGLE_MESH_HPP

#include "config.h"
#include <assert.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <valarray>
#ifdef DIFF_DEFINE
#include <limits>
#include <set>
#endif /* DIFF_DEFINE */

#ifdef USE_HASH_MAP
#   include <ext/hash_map>
#else
#   include <unordered_map>
#endif

#include "Triangle.hpp"
#include "Point3.hpp"
#include "linearalgebra/Vector3.hpp"
#include <stdexcept>
#include <utils/STL_Wrapper.h>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class TriangleMesh
{
    public:
        struct NeighborRec
        {
            int num;    // how many neighbor faces
            int id[3];
        };

        TriangleMesh()
            : m_totArea(0)
        {
        }

        // needed since now it can be a base class to mesh class that can
        // search for nearest neighbours.
        virtual ~TriangleMesh(){}

        void clear()
        {
            m_vertices.clear();
            m_normals.clear();
            m_triangles.clear();
            m_totArea = 0;
        }

        //! Return the ID of the added vertices
        int add_vertex(const Point3<T>& vtx)
        {
            m_vertices.push_back(vtx);
            return m_vertices.size() - 1;
        }

        //! Return the ID of added vertices
        int add_vertex_normal(const Point3<T>& vtx, const Vector3<T>& nml)
        {
            m_vertices.push_back(vtx);
            m_normals.push_back(nml);
            return m_vertices.size() - 1;
        }

        //! Return the current number of triangles
        int add_triangle(unsigned int v0, unsigned int v1, unsigned int v2)
        {
            assert(v0 != v1 && v0 != v2 && v1 != v2);
            assert(v0 < m_vertices.size() &&
                   v1 < m_vertices.size() &&
                   v2 < m_vertices.size());
            m_triangles.push_back(Tuple3ui(v0, v1, v2));
            return m_triangles.size();
        }

        int num_vertices() const { return m_vertices.size(); }
        int num_triangles() const { return m_triangles.size(); }

#ifdef DIFF_DEFINE
        std::vector< Point3<T> >& vertices() 
        {  return m_vertices; }

        std::vector<Tuple3ui>& triangles()
        {  return m_triangles; }

#endif /* DIFF_DEFINE */
        const std::vector< Point3<T> >& vertices() const
        {  return m_vertices; }

        const std::vector<Tuple3ui>& surface_indices() const
        {  return m_triangles; }

        const std::vector<Tuple3ui>& triangles() const
        {  return m_triangles; }

        const std::vector< Vector3<T> >& normals() const
        {  return m_normals; }

        const std::valarray<T>& vertex_areas() const
        {  return m_vtxAreas; }

#ifdef DIFF_DEFINE
        const Tuple3ui& triangle_ids(int tid) const
        {
            assert(tid < m_triangles.size());
            return m_triangles[tid];
        }

#endif /* DIFF_DEFINE */
        const Point3<T>& vertex(size_t vid) const
        {  
            assert(vid < m_vertices.size());
            return m_vertices[vid];
        }

#ifdef DIFF_DEFINE
        const Vector3<T>& normal(size_t vid) const
        {
            assert(vid < m_vertices.size() && m_normals.size() == m_vertices.size());
            return m_normals[vid];
        }

        const Vector3<T>* normal_ptr(size_t vid) const
        {
            assert(vid < m_vertices.size() && m_normals.size() == m_vertices.size());
            return &(m_normals[vid]);
        }

#endif /* DIFF_DEFINE */
        bool has_normals() const 
        {
            return !m_normals.empty() && m_vertices.size() == m_normals.size();
        }

        bool empty() const 
        {  return m_vertices.empty(); }

        double total_area() const 
        {  return m_totArea; }

        // Find the radius of the bounding sphere enclosing this mesh and centered
        // at the given point
        T boundingSphereRadius( const Point3<T> &center ) const;

#ifdef DIFF_DEFINE
        void bounding_box(Point3<T>& low, Point3<T>& up) const;

#endif /* DIFF_DEFINE */
        void generate_normals();
#ifndef DIFF_DEFINE
        //! save the current mesh to a file
        int save_mesh(const char* file) const;
        int load_mesh(const char* file);
#else /* DIFF_DEFINE */
        void generate_pseudo_normals();
        //! Save the current mesh to a file
        int save_mesh_txt(const char* file) const;
        int load_mesh_txt(const char* file);
#endif /* DIFF_DEFINE */
        void update_vertex_areas();

        void translate_x(T dx);
        void translate_y(T dy);
        void translate_z(T dz);

#ifndef DIFF_DEFINE
        //! get the neighborship of faces
#else /* DIFF_DEFINE */
        //! Get the triangle face neighbors; for each triangle, return
        //  a list of its neighbor triangles
#endif /* DIFF_DEFINE */
        void get_face_neighborship(std::vector<NeighborRec>&) const;
#ifdef DIFF_DEFINE
        //! Get the table of neighbors of all vertices
        void get_vtx_neighborship(std::vector< std::set<int> >&) const;
        //! Get the adjacent triangles of each vertex
        void get_vtx_tgls(std::vector< std::set<int> >&) const;
#endif /* DIFF_DEFINE */

        // abstract nearest triangles lookup that should be specific to data
        // structures of choice. 
        virtual REAL FindKNearestTriangles(const int &k, const Vector3d &point, std::vector<int> &triangleIndices)
        { throw std::runtime_error("**ERROR** Nearest neighbour search for class TriangleMesh not implemented."); }
        virtual REAL FindNearestTriangle(const Vector3d &point, int &triangleIndex)
        { throw std::runtime_error("**ERROR** Nearest neighbour search for class TriangleMesh not implemented."); }

        // need FindKNearestTriangles()
        REAL ComputeClosestPointOnMesh(const Vector3d &queryPoint, Vector3d &closestPoint, int &closestTriangle, Vector3d &projectedPoint); 

    protected:
        double                      m_totArea;
        std::vector< Point3<T> >    m_vertices;
        std::vector< Vector3<T> >   m_normals;
        std::vector<Tuple3ui>       m_triangles;    // indices of triangle vertices
        std::valarray<T>            m_vtxAreas;     // area of each triangles
};

template <typename T>
void TriangleMesh<T>::translate_x(T dx)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(dx)
    for(int i = 0;i < (int)m_vertices.size();++ i)
#else
    for(size_t i = 0;i < m_vertices.size();++ i)
#endif
        m_vertices[i].x += dx;
}

template <typename T>
void TriangleMesh<T>::translate_y(T dy)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(dy)
    for(int i = 0;i < (int)m_vertices.size();++ i)
#else
    for(size_t i = 0;i < m_vertices.size();++ i)
#endif
        m_vertices[i].y += dy;
}

template <typename T>
void TriangleMesh<T>::translate_z(T dz)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(dz)
    for(int i = 0;i < (int)m_vertices.size();++ i)
#else
    for(size_t i = 0;i < m_vertices.size();++ i)
#endif
        m_vertices[i].z += dz;
}

template <typename T>
void TriangleMesh<T>::get_face_neighborship(std::vector<NeighborRec>& neighbors) const
{
#ifdef USE_UNORDERED_MAP
    //using namespace std::tr1;
    using namespace std;
    unordered_map< int, unordered_map<int, int> > hash;
#else
    using namespace __gnu_cxx;
    hash_map< int, hash_map<int, int> > hash;
#endif
    //using namespace std;
    //map< int, map<int, int> > hash;
    int v0, v1;

    neighbors.resize(m_triangles.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(neighbors)
    for(int i = 0;i < (int)m_triangles.size();++ i)
#else
    for(size_t i = 0;i < m_triangles.size();++ i)
#endif
        neighbors[i].num = 0;

    for(size_t i = 0;i < m_triangles.size();++ i)
    {
        // edge 0
        v0 = m_triangles[i][0];
        v1 = m_triangles[i][1];
        if ( v0 > v1 ) std::swap(v0, v1);   // v0 < v1
        if ( hash.count(v0) && hash[v0].count(v1) )
        {
            int fid = hash[v0][v1];
            neighbors[i].id[neighbors[i].num ++] = fid;
            neighbors[fid].id[neighbors[fid].num ++] = i;
        }
        else
            hash[v0][v1] = i;

        // edge 1
        v0 = m_triangles[i][1];
        v1 = m_triangles[i][2];
        if ( v0 > v1 ) std::swap(v0, v1);
        if ( hash.count(v0) && hash[v0].count(v1) )
        {
            int fid = hash[v0][v1];
            neighbors[i].id[neighbors[i].num ++] = fid;
            neighbors[fid].id[neighbors[fid].num ++] = i;
        }
        else
            hash[v0][v1] = i;

        // edge 2
        v0 = m_triangles[i][2];
        v1 = m_triangles[i][0];
        if ( v0 > v1 ) std::swap(v0, v1);
        if ( hash.count(v0) && hash[v0].count(v1) )
        {
            int fid = hash[v0][v1];
            neighbors[i].id[neighbors[i].num ++] = fid;
            neighbors[fid].id[neighbors[fid].num ++] = i;
        }
        else
            hash[v0][v1] = i;
    }
}

template <typename T>
void TriangleMesh<T>::update_vertex_areas()
{
    const T s = (T)1 / (T)3;
    m_totArea = 0;
    m_vtxAreas.resize(m_vertices.size(), (T)0);
    for(int i = 0;i < (int)m_triangles.size();++ i)
    {
        T area = Triangle<T>::area(
                m_vertices[m_triangles[i][0]],
                m_vertices[m_triangles[i][1]],
                m_vertices[m_triangles[i][2]]); 
        m_totArea += area;
        area *= s;

        for(int j = 0;j < 3;++ j)
            m_vtxAreas[m_triangles[i][j]] += area;
    }
}

template <typename T>
#ifdef DIFF_DEFINE
void TriangleMesh<T>::generate_pseudo_normals()
{
    std::vector< Vector3<T> > tglnmls(m_triangles.size());
    for(size_t i = 0;i < m_triangles.size();++ i)
    {
        tglnmls[i] = Triangle<T>::normal(
                m_vertices[m_triangles[i][0]],
                m_vertices[m_triangles[i][1]],
                m_vertices[m_triangles[i][2]]);
        if ( tglnmls[i].lengthSqr() < 1E-24 )
        {
            fprintf(stderr, "ERROR: triangle has zero area: %.30g\n",
                    tglnmls[i].lengthSqr());
            exit(1);
        }
        tglnmls[i].normalize();
    }

    m_normals.resize(m_vertices.size());
    memset(&m_normals[0], 0, sizeof(Vector3<T>)*m_vertices.size());
    for(size_t i = 0;i < m_triangles.size();++ i)
    {
        const Vector3<T>& nml = tglnmls[i];
        m_normals[m_triangles[i][0]] += nml * Triangle<T>::angle(
                m_vertices[m_triangles[i][2]],
                m_vertices[m_triangles[i][0]],
                m_vertices[m_triangles[i][1]]);
        m_normals[m_triangles[i][1]] += nml * Triangle<T>::angle(
                m_vertices[m_triangles[i][0]],
                m_vertices[m_triangles[i][1]],
                m_vertices[m_triangles[i][2]]);
        m_normals[m_triangles[i][2]] += nml * Triangle<T>::angle(
                m_vertices[m_triangles[i][1]],
                m_vertices[m_triangles[i][2]],
                m_vertices[m_triangles[i][0]]);
    }

    for(size_t i = 0;i < m_vertices.size();++ i)
        m_normals[i].normalize();
}

template <typename T>
#endif /* DIFF_DEFINE */
void TriangleMesh<T>::generate_normals()
{
    m_vtxAreas.resize(m_vertices.size(), (T)0);
    m_normals.resize(m_vertices.size());
    memset(&m_normals[0], 0, sizeof(Vector3<T>)*m_vertices.size());
    m_totArea = 0;

    for(int i = 0;i < (int)m_triangles.size();++ i)
    {
        Point3<T>& p0 = m_vertices[m_triangles[i][0]];
        Point3<T>& p1 = m_vertices[m_triangles[i][1]];
        Point3<T>& p2 = m_vertices[m_triangles[i][2]];
        Vector3<T> nml = (p1 - p0).crossProduct(p2 - p0) * 0.5;  // weight = area
        T area = nml.length();
        m_totArea += area;

        for(int j = 0;j < 3;++ j)
        {
            m_vtxAreas[m_triangles[i][j]] += area;
            m_normals[m_triangles[i][j]] += nml;
        }
    }

    const T alpha = (T)1 / (T)3;
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 10000)
    for(int i = 0;i < (int)m_vertices.size();++ i)
#else
    for(int i = 0;i < (int)m_vertices.size();++ i)
#endif
    {
#ifndef DIFF_DEFINE
        m_normals[i] /= m_vtxAreas[i];
#endif /* ! DIFF_DEFINE */
        m_normals[i].normalize();
        m_vtxAreas[i] *= alpha;
    }
}

/*!
 * number of vertices [N/F] number of triangles
 * v0.x v0.y v0.z
 * v1.x v1.y v1.z
 * ......
 * triangle0.v0 triangle0.v1 triangle0.v2
 * triangle1.v0 triangle1.v1 triangle1.v2
 * triangle2.v0 triangle2.v1 triangle2.v2
 * ......
 */
template <typename T>
#ifndef DIFF_DEFINE
int TriangleMesh<T>::save_mesh(const char* file) const
#else /* DIFF_DEFINE */
int TriangleMesh<T>::save_mesh_txt(const char* file) const
#endif /* DIFF_DEFINE */
{
    using namespace std;

    ofstream fout(file);
    if ( !fout.good() ) return -1;

    bool nml_out;
    fout << m_vertices.size();
    if ( !m_vertices.empty() && m_vertices.size() == m_normals.size() )
    {
        fout << " N ";
        nml_out = true;
    }
    else
    {
        fout << " F ";
        nml_out = false;
    }
    fout << m_triangles.size() << endl;

    fout << setprecision(16);
    for(size_t i = 0;i < m_vertices.size();++ i)
        fout << m_vertices[i].x 
             << ' ' << m_vertices[i].y 
             << ' ' << m_vertices[i].z << endl;
    if ( nml_out )
        for(size_t i = 0;i < m_normals.size();++ i)
            fout << m_normals[i].x << ' ' 
                 << m_normals[i].y << ' ' 
                 << m_normals[i].z << endl;
    for(size_t i = 0;i < m_triangles.size();++ i)
        fout << m_triangles[i][0] << ' ' 
             << m_triangles[i][1] << ' ' 
             << m_triangles[i][2] << endl;

    fout.flush();
    fout.close();
    return 0;
}

template <typename T>
#ifndef DIFF_DEFINE
int TriangleMesh<T>::load_mesh(const char* file)
#else /* DIFF_DEFINE */
int TriangleMesh<T>::load_mesh_txt(const char* file)
#endif /* DIFF_DEFINE */
{
    using namespace std;

    ifstream fin(file);
    if ( !fin.good() ) return -1;

    int nvtx, ntgl;
    char C;
    fin >> nvtx >> C >> ntgl;
    clear();
    T a, b, c;
    unsigned int v0, v1, v2;
    for(int i = 0;i < nvtx;++ i)
    {
        fin >> a >> b >> c;
        m_vertices.push_back(Point3<T>(a, b, c));
    }
    if ( fin.fail() ) return -1;
    if ( C == 'N' )
    {
        for(int i = 0;i < nvtx;++ i)
        {
            fin >> a >> b >> c;
            m_normals.push_back(Vector3<T>(a, b, c));
        }
        if ( fin.fail() ) return -1;
    }
    for(int i = 0;i < ntgl;++ i)
    {
        fin >> v0 >> v1 >> v2;
        m_triangles.push_back(Tuple3ui(v0, v1, v2));
    }
    if ( fin.fail() ) return -1;

    fin.close();
    return 0;
#ifdef DIFF_DEFINE
}

// Find the radius of the bounding sphere enclosing this mesh and centered
// at the given point
template <typename T>
T TriangleMesh<T>::boundingSphereRadius( const Point3<T> &center ) const
{
    Vector3<T>               diff;
    T                        maxRadius = 0.0;

    for ( int vertex_id = 0; vertex_id < m_vertices.size(); vertex_id++ ) {
        diff = m_vertices[ vertex_id ] - center;

        maxRadius = max( maxRadius, diff.length() );
    }

    return maxRadius;
}

template <typename T>
void TriangleMesh<T>::bounding_box(Point3<T>& low, Point3<T>& up) const
{
    //const T inf = std::numeric_limits<T>::infinity();
    const T maxValue = std::numeric_limits<T>::max(); 
    const T minValue = std::numeric_limits<T>::lowest(); 
    low = Point3<T>(maxValue, maxValue, maxValue);
    up  = Point3<T>(minValue, minValue, minValue);

    const typename std::vector< Point3<T> >::const_iterator end = m_vertices.end();
    for(typename std::vector< Point3<T> >::const_iterator it = m_vertices.begin();
            it != end;++ it)
    {
        low.x = min(low.x, it->x);
        low.y = min(low.y, it->y);
        low.z = min(low.z, it->z);

        up.x = max(up.x, it->x);
        up.y = max(up.y, it->y);
        up.z = max(up.z, it->z);
    }
}

template <typename T>
void TriangleMesh<T>::get_vtx_neighborship(std::vector< std::set<int> >& tbl) const
{
    tbl.resize(m_vertices.size());
    for(size_t i = 0;i < m_vertices.size();++ i) tbl[i].clear();

    const std::vector<Tuple3ui>::const_iterator end = m_triangles.end();
    for(std::vector<Tuple3ui>::const_iterator it = m_triangles.begin();
            it != end;++ it)
    {
        tbl[it->x].insert(it->y); tbl[it->x].insert(it->z);
        tbl[it->y].insert(it->x); tbl[it->y].insert(it->z);
        tbl[it->z].insert(it->x); tbl[it->z].insert(it->y);
    }
}

/*
 * Return a list of triangles adjacent to each vertex
 */
template <typename T>
void TriangleMesh<T>::get_vtx_tgls(std::vector< std::set<int> >& vtx_ts) const
{
    vtx_ts.resize(m_vertices.size());
    for(size_t i = 0;i < m_triangles.size();++ i)
    {
        vtx_ts[ m_triangles[i].x ].insert(i);
        vtx_ts[ m_triangles[i].y ].insert(i);
        vtx_ts[ m_triangles[i].z ].insert(i);
    }
#endif /* DIFF_DEFINE */
}

/* 
 * This function finds a point on the mesh that is closest to the queryPoint. 
 *
 * First a list of N_neighbours triangles that have distance closest to the queryPoint 
 * is retrived. How distance metric is defined depends on the KNN search, if 
 * TriangleMeshKDTree class is used then its the quadratic L2 distance between triangle 
 * centroid and queryPoint. Next, the point is projected onto the plane defined by the 
 * triangle vertices and its barycentric coordinates (u, v) were computed. We can then 
 * reason about the closest feature (vertex, edge, or interior) using the (u, v) 
 * coordinates. Notations are similar to the website:
 *
 * www.blackpawn.com/texts/pointinpoly
 * 
 */
template <typename T> 
REAL TriangleMesh<T>::
ComputeClosestPointOnMesh(const Vector3d &queryPoint, Vector3d &closestPoint, int &closestTriangle, Vector3d &projectedPoint)
{
    const int N_neighbours = 10; 
    std::vector<int> triangleIndices; 
    FindKNearestTriangles(N_neighbours, queryPoint, triangleIndices); 
    assert(triangleIndices.size() > 0);

    REAL minDistance = std::numeric_limits<REAL>::max(); 
    Vector3d closestPointBuffer, projectedPointBuffer;
    for (const int t_idx : triangleIndices)
    {
        // points CCW ordering
        const Tuple3ui &triangle = this->triangle_ids(t_idx); 
        const Point3<T> &p0 = this->vertex(triangle.x); 
        const Point3<T> &p1 = this->vertex(triangle.y); 
        const Point3<T> &p2 = this->vertex(triangle.z); 

        // edges w.r.t p0
        const Vector3<T> e0 = p2 - p0; 
        const Vector3<T> e1 = p1 - p0; 

        // dot products
        const REAL dot00 = e0.dotProduct(e0); 
        const REAL dot01 = e0.dotProduct(e1); 
        const REAL dot11 = e1.dotProduct(e1); 

        // normal n corresponds to CCW point ordering and plane spanned by three
        // points n_x x + n_y y + n_z z = d, for all (x, y, z) in this plane.
        const Vector3<T> triangleNormal = e1.crossProduct(e0); 
        const REAL d = triangleNormal.x * p0.x + triangleNormal.y * p0.y + triangleNormal.z * p0.z; 

        // to find projection of query point on this plane, first find the time
        // travelled t for the ray emitted from the query point to this plane.
        const REAL t = (d - triangleNormal.dotProduct(queryPoint)) / (triangleNormal.lengthSqr()); 
        projectedPointBuffer = queryPoint + triangleNormal * t;

        // compute barycentric coordinates of this projected point (u, v)
        // projectedPoint = p0 + u*(p2 - p0) + v*(p1 - p0)
        const Vector3<T> e2 = projectedPointBuffer - p0; 
        const REAL dot02 = e0.dotProduct(e2); 
        const REAL dot12 = e1.dotProduct(e2); 
        const REAL invDenom = 1. / (dot00 * dot11 - dot01 * dot01); 
        const REAL u = (dot11 * dot02 - dot01 * dot12) * invDenom; 
        const REAL v = (dot00 * dot12 - dot01 * dot02) * invDenom; 

        // distinguish cases and compute distance, closestPoint
        // There are six cases other than lie inside triangle
        // closest to points p0, p1, p2, 
        // or edges (p0, p1), (p1, p2), (p0, p2)
        if (u>=0 && v>=0 && (u+v)<1) // inside triangle
        {
            closestPointBuffer = Vector3<T>(projectedPointBuffer); 
        }
        else if (u<0 && v>=0 && (u+v)<1) // closest to edge e = e1(p0, p1)
        {
            const REAL l = std::max<REAL>(dot12, 0.0);
            const REAL e1Sqr = e1.lengthSqr(); 
            if (l > e1Sqr)
                closestPointBuffer = projectedPointBuffer + (e1 - e2);
            else
                closestPointBuffer = projectedPointBuffer + (e1 * l/e1Sqr - e2);
        }
        else if (u>=0 && v<0 && (u+v)<1) // closest to edge e = e0(p0, p2)
        {
            const REAL l = std::max<REAL>(dot02, 0.0);
            const REAL e0Sqr = e0.lengthSqr(); 
            if (l > e0Sqr)
                closestPointBuffer = projectedPointBuffer + (e0 - e2);
            else
                closestPointBuffer = projectedPointBuffer + (e0 * l/e0Sqr - e2);
        }
        else if (u>=0 && v>=0 && (u+v)>=1) // closest to edge e = e(p1, p2)
        {
            const Vector3<T> e12 = p2 - p1; 
            const Vector3<T> o12 = projectedPointBuffer - p1; 
            const REAL l = std::max<REAL>(o12.dotProduct(e12), 0.0);
            const REAL e12Sqr = e12.lengthSqr(); 
            if (l > e12Sqr)
                closestPointBuffer = projectedPointBuffer + (e12 - o12);
            else
                closestPointBuffer = projectedPointBuffer + (e12 * l/e12Sqr - o12);
        }
        else if (u<0 && v<0 && (u+v)<1) // closest to p0
        {
            closestPointBuffer = Vector3<T>(p0); 
        }
        else if (u<0 && v>=0 && (u+v)>=1) // closest to p1
        {
            closestPointBuffer = Vector3<T>(p1); 
        }
        else if (u>=0 && v<0 && (u+v)>=1) // closest to p2
        {
            closestPointBuffer = Vector3<T>(p2); 
        }
        else // remaining case is u<0, v<0, (u+v)>=1 but this case is impossible
        {
            throw std::runtime_error("**ERROR** Barycentric coordinates computation yields impossible case");
        }
        const REAL distance = (queryPoint - closestPointBuffer).length();

        // keep closest
        if (distance < minDistance) 
        {
            minDistance = distance; 
            closestTriangle = t_idx; 
            closestPoint = closestPointBuffer;
            projectedPoint = projectedPointBuffer; 
        }
    }
    return minDistance; 
}

#ifdef USE_NAMESPACE
}
#endif

#endif
