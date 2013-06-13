#ifndef DIFF_DEFINE
/******************************************************************************
 *  File: TglMeshWriter.hpp
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
#ifndef GEOMETRIC_MESH_WRITER
#   define GEOMETRIC_MESH_WRITER

#include <fstream>
#include "geometry/TriangleMesh.hpp"
#include "utils/macros.h"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

/*
 * Write into RTGI MDL format
 */
class MeshMdlWriter
{
    public:
        template <typename T>
        static int write(const TriangleMesh<T>& mesh, const char* file);
};

template <typename T>
int MeshMdlWriter::write(const TriangleMesh<T>& mesh, const char* file)
{
    using namespace std;

    ofstream fout(file);
    if ( !fout.good() ) return ERROR_RETURN;

    fout << "msh \"$NAME$\"" << endl;
    fout << "  $MTRL$" << endl;

    const vector< Point3<T> >& vtx = mesh.vertices();
    fout << "  vrtxPstn" << endl;
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "    " << vtx[i].x << ' ' << vtx[i].z << ' ' << -vtx[i].y << endl;
    fout << "  end" << endl;

    if ( mesh.has_normals() )
    {
        const vector< Vector3<T> >& nml = mesh.normals();
        fout << "  vrtxNrml" << endl;
        for(size_t i = 0;i < nml.size();++ i)
            fout << "    " << nml[i].x << ' ' << nml[i].z << ' ' << -nml[i].y << endl;
        fout << "  end" << endl;
    }
    
    const vector<Tuple3ui>& tgl = mesh.triangles();
    fout << "  trngl" << endl;
    for(size_t i = 0;i < tgl.size();++ i)
        fout << "    " << tgl[i].x << ' ' << tgl[i].y << ' ' << tgl[i].z << endl;
    fout << "  end" << endl;
    fout << "end" << endl;
    fout.close();
    return SUCC_RETURN;
}

/*
 * write a mesh as a wavefront obj file
 */
class MeshObjWriter
{
    public:
        template <typename T>
        static int write(const TriangleMesh<T>& mesh, const char* file);

        /*
         * write the area associated with each vertex also into the OBJ file.
         *
         * This is not a "standard" obj file format
         */
        template <typename T>
        static int write_with_vtx_weight(const TriangleMesh<T>& mesh, 
                                         const char* file);
};

template <typename T>
int MeshObjWriter::write_with_vtx_weight(const TriangleMesh<T>& mesh, const char* file)
{
    using namespace std;

    ofstream fout(file);
    fout << setprecision(12);
    fout << "# OBJ Mesh: Generated by MeshObjWriter (Changxi Zheng)" << endl;
    const vector< Point3<T> >& vtx = mesh.vertices();
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "v " << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z << endl;
    if ( mesh.has_normals() )
    {
        const vector< Vector3<T> >& nml = mesh.normals();
        for(size_t i = 0;i < nml.size();++ i)
            fout << "vn " << nml[i].x << ' ' << nml[i].y << ' ' << nml[i].z << endl;
        fout << "g all" << endl;
        fout << "s 0" << endl;
    }
    else
    {
        fout << "g all" << endl;
        fout << "s 1" << endl;
    }
    const vector<Tuple3ui>& tgl = mesh.triangles();
    for(size_t i = 0;i < tgl.size();++ i)
        fout << "f " << tgl[i].x+1 << ' ' << tgl[i].y+1 << ' ' << tgl[i].z+1 << endl;

    const valarray<T>& areas = mesh.vertex_areas();
    for(size_t i = 0;i < areas.size();++ i)
        fout << "w " << areas[i] << endl;

    fout.close();
    return SUCC_RETURN;
}

template <typename T>
int MeshObjWriter::write(const TriangleMesh<T>& mesh, const char* file)
{
    using namespace std;

    ofstream fout(file);
    if ( fout.fail() )
    {
        cerr << "ERROR: cannot open file: " << file << " for writing" << endl;
        return ERROR_RETURN;
    }

    fout << "# OBJ Mesh: Generated by MeshObjWriter (Changxi Zheng)" << endl;
    fout << setprecision(12);
    const vector< Point3<T> >& vtx = mesh.vertices();

    for(size_t i = 0;i < vtx.size();++ i)
        fout << "v " << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z << endl;

    if ( mesh.has_normals() )
    {
        const vector< Vector3<T> >& nml = mesh.normals();
        for(size_t i = 0;i < nml.size();++ i)
            fout << "vn " << nml[i].x << ' ' << nml[i].y << ' ' << nml[i].z << endl;
        fout << "g all" << endl;
        fout << "s 0" << endl;
    }
    else
    {
        fout << "g all" << endl;
        fout << "s 1" << endl;
    }
    const vector<Tuple3ui>& tgl = mesh.triangles();
    for(size_t i = 0;i < tgl.size();++ i)
        fout << "f " << tgl[i].x+1 << ' ' << tgl[i].y+1 << ' ' << tgl[i].z+1 << endl;

    fout.close();
    return SUCC_RETURN;
}

///////////////////////////////////////////////////////////////////////////////
/*
 * Write a mesh as the renderman RIB file
 */
class MeshRibWriter
{
    public:
        template <typename T>
        static int write(const TriangleMesh<T>& mesh, const char* file);

        template <typename T>
        static int write(const TriangleMesh<T>& mesh, std::ofstream& fout);

        /* 
         * write mesh with constant normal
         * the given normal would override the normals in the given mesh
         */
        template <typename T>
        static int write(const TriangleMesh<T>& mesh, const Vector3<T>& nml, const char* file);
};

template <typename T>
int MeshRibWriter::write(const TriangleMesh<T>& mesh, std::ofstream& fout)
{
    using namespace std;

    const vector< Point3<T> >& vtx = mesh.vertices();
    const vector<Tuple3ui>&    tgl = mesh.triangles();
    fout << "PointsPolygons [";
    for(size_t i = 0;i < tgl.size();++ i)
    {
        if ( i % 40 == 0 ) fout << endl << ' ';
        fout << " 3";
    }
    fout << ']' << endl;
    fout << '[' << endl;
    for(size_t i = 0;i < tgl.size();++ i)
        fout << "  " << tgl[i][0] << ' ' << tgl[i][1] << ' ' << tgl[i][2] << endl;
    fout << ']' << endl;

    //// vertices' position
    fout << "\"P\" [" << endl;
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "  " << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z << endl;
    fout << ']' << endl;

    //// output normals
    if ( mesh.has_normals() )
    {
        fout << "\"N\" [" << endl;
        const vector< Vector3<T> >& nml = mesh.normals();
        for(size_t i = 0;i < nml.size();++ i)
            fout << "  " << nml[i].x << ' ' << nml[i].y << ' ' << nml[i].z << endl;
        fout << ']' << endl;
    }

    return SUCC_RETURN;
}

template <typename T>
int MeshRibWriter::write(const TriangleMesh<T>& mesh, const char* file)
{
    using namespace std;

    ofstream fout(file);
    fout << setprecision(10);
    fout << "## Renderman RIB  Generated by MeshRibWriter (Changxi Zheng)" << endl;
    fout << "version 3.04" << endl;
    fout << "Declare \"N\" \"varying normal\"" << endl;

    const vector< Point3<T> >& vtx = mesh.vertices();
    const vector<Tuple3ui>&    tgl = mesh.triangles();
    fout << "PointsPolygons [";
    for(size_t i = 0;i < tgl.size();++ i)
    {
        if ( i % 40 == 0 ) fout << endl << ' ';
        fout << " 3";
    }
    fout << ']' << endl;
    fout << '[' << endl;
    for(size_t i = 0;i < tgl.size();++ i)
        fout << "  " << tgl[i][0] << ' ' << tgl[i][1] << ' ' << tgl[i][2] << endl;
    fout << ']' << endl;

    //// vertices' position
    fout << "\"P\" [" << endl;
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "  " << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z << endl;
    fout << ']' << endl;

    //// output normals
    if ( mesh.has_normals() )
    {
        fout << "\"N\" [" << endl;
        const vector< Vector3<T> >& nml = mesh.normals();
        for(size_t i = 0;i < nml.size();++ i)
            fout << "  " << nml[i].x << ' ' << nml[i].y << ' ' << nml[i].z << endl;
        fout << ']' << endl;
    }

    fout.close();
    return SUCC_RETURN;
}

template <typename T>
int MeshRibWriter::write(const TriangleMesh<T>& mesh, const Vector3<T>& nml, const char* file)
{
    using namespace std;

    ofstream fout(file);
    fout << setprecision(10);
    fout << "## Renderman RIB  Generated by MeshRibWriter (Changxi Zheng)" << endl;
    fout << "version 3.04" << endl;
    fout << "Declare \"N\" \"varying normal\"" << endl;

    const vector< Point3<T> >& vtx = mesh.vertices();
    const vector<Tuple3ui>&    tgl = mesh.triangles();
    fout << "PointsPolygons [";
    for(size_t i = 0;i < tgl.size();++ i)
    {
        if ( i % 40 == 0 ) fout << endl << ' ';
        fout << " 3";
    }
    fout << ']' << endl;
    fout << '[' << endl;
    for(size_t i = 0;i < tgl.size();++ i)
        fout << "  " << tgl[i][0] << ' ' << tgl[i][1] << ' ' << tgl[i][2] << endl;
    fout << ']' << endl;

    //// vertices' position
    fout << "\"P\" [" << endl;
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "  " << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z << endl;
    fout << ']' << endl;

    //// output normals
    fout << "\"N\" [" << endl;
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "  " << nml.x << ' ' << nml.y << ' ' << nml.z << endl;
    fout << ']' << endl;

    fout.close();
    return SUCC_RETURN;
}
#ifdef USE_NAMESPACE
}
#endif
#endif
