#ifndef RENDERER_TET_MESH_HPP
#   define RENDERER_TET_MESH_HPP

#include <GL/gl.h>
#include "geometry/TetMesh.hpp"

/*!
 * Display tetrahedron meshes using OpenGL
 */
class TetMeshRenderer
{
    public:
        TetMeshRenderer():m_selected(-1), m_color(0,1,0) { }

        template <typename T>
        void render(const TetMesh<T>& mesh) const;

        template <typename T>
        void render_with_names(const TetMesh<T>& mesh) const;

        void set_selected(int s) 
        {
            m_selected = s;
        }

        int selected() const
        { return m_selected; }

        inline void set_color(const Tuple3f& c)
        {
            m_color = c;
            m_color.clamp(Tuple3f(0, 0, 0), Tuple3f(1, 1, 1));
        }

        void set_color(float r, float g, float b)
        {
            m_color.set(r, g, b);
            m_color.clamp(Tuple3f(0, 0, 0), Tuple3f(1, 1, 1));
        }

    private:
        void load_vertices(const Tuple3d* vtxPtr, const Tuple3d* nmlPtr) const
        {
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);

            glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)vtxPtr);
            glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)nmlPtr);
        }
        void load_vertices(const Tuple3f* ptr) const;

        /*! selected triangles by mouse click */
        int         m_selected;
        Vector3f    m_color;
};

///////////////////////////////////////////////////////////////////////////////

/*
 * Render double precison tet mesh
 */
template <typename T>
void TetMeshRenderer::render(const TetMesh<T>& mesh) const
{
    using namespace std;

    const vector< Point3<T> >&  vtx = mesh.vertices();
    const vector<Tuple3ui>&     tgl = mesh.surface_indices();
    const vector< Vector3<T> >& nml = mesh.normals();
    load_vertices(&vtx[0], &nml[0]);

    glColor3fv(m_color);
    glEnable(GL_LIGHTING);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_TRIANGLES, 
            tgl.size()*3, 
            GL_UNSIGNED_INT,
            (const GLvoid*)&tgl[0]);

    if ( m_selected >= 0 && m_selected < (int)tgl.size() )
    {
        glColor3f(1, 0, 0);
        GLint polyMode;
        glGetIntegerv(GL_POLYGON_MODE, &polyMode);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT,
                (const GLvoid*)&tgl[m_selected]);
        glPolygonMode(GL_FRONT_AND_BACK, polyMode);
    }
}

template <typename T>
void TetMeshRenderer::render_with_names(const TetMesh<T>& mesh) const
{
    using namespace std;

    const vector< Point3<T> >&  vtx = mesh.vertices();
    const vector<Tuple3ui>&     tgl = mesh.surface_indices();
    const vector< Vector3<T> >& nml = mesh.normals();
    load_vertices(&vtx[0], &nml[0]);
    
    for(size_t i = 0;i < tgl.size();++ i)
    {
        glPushName(i);
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT,
                (const GLvoid*)&tgl[i]);
        glPopName();
    }
}

#endif

