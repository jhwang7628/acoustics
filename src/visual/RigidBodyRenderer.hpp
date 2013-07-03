#ifndef RIGID_BODY_RENDERER_HPP
#   define RIGID_BODY_RENDERER_HPP

#include <GL/gl.h>
#include "colormap/ColorMap.h"
#include "rigid/LSCollisionRigidBody.hpp"
#include "geometry/TetMesh.hpp"

class RigidBodyRenderer
{
    public:
        typedef JetColorMap                          TColorMap;

    public:
        RigidBodyRenderer():color_(0.2, 0.2, 0.2) 
        { }

        void set_color(float r, float g, float b)
        {   color_.set(r, g, b); }

        template <class TMesh>
        void render(const LSCollisionRigidBody<double, TMesh>* rbody) const;

        template <class TMesh>
        void render(const TetMesh<double>* mesh,
                    const LSCollisionRigidBody<double, TMesh>* rbody) const;
        
        template <class TMesh>
        void render(const LSCollisionRigidBody<double, TMesh>* rbody,
                    const Vector3<double>& x, const Quaternion<double>& rot) const;

    private:
        Vector3f    color_;
        TColorMap   cmap_;
};

template <class TMesh>
void RigidBodyRenderer::render(const TetMesh<double>* mesh,
        const LSCollisionRigidBody<double, TMesh>* rbody) const
{
    using namespace std;

    const vector< Point3d >&  vtx = mesh->vertices();
    const vector<Tuple3ui>&   tgl = mesh->surface_indices();
    const vector< Vector3d >& nml = mesh->normals();

    const Point3<double>&     x0 = rbody->initial_mass_center();
    const Point3<double>&     x  = rbody->mass_center();
    const Quaternion<double>& rq = rbody->rotation();
    Vector3<double> rvec;
    REAL rd = rq.toAxisRotD(rvec);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)&vtx[0]);
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)&nml[0]);

    glColor3fv(color_);
    glEnable(GL_LIGHTING);

    glPushMatrix();
    glTranslated(x.x, x.y, x.z);
    glRotated(rd, rvec.x, rvec.y, rvec.z);
    glTranslated(-x0.x, -x0.y, -x0.z);
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT,
                   (const GLvoid*)&tgl[0]);
    glPopMatrix();
}

template <class TMesh>
void RigidBodyRenderer::render(const LSCollisionRigidBody<double, TMesh>* rbody) const
{
    using namespace std;

    const TriangleMesh<double>*    pmesh = rbody->collision_processor()->surface_mesh();    // surface triangles
    const vector<Tuple3ui>&     tgl = pmesh->surface_indices();                             // surface triangles
    const vector< Point3<double> >&  vtx = pmesh->vertices();                               // surface triangles
    const vector< Vector3<double> >& nml = rbody->collision_processor()->vtx_pseudo_nml();  // surface triangles
    const Point3<double>&     x0 = rbody->initial_mass_center();
    const Point3<double>&     x  = rbody->mass_center();
    const Quaternion<double>& rq = rbody->rotation();
    Vector3<double> rvec;
    REAL rd = rq.toAxisRotD(rvec);
    glEnable(GL_LIGHTING);

    glColor3fv(color_);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)&vtx[0]);
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)&nml[0]);

    glPushMatrix();
    glTranslated(x.x, x.y, x.z);
    glRotated(rd, rvec.x, rvec.y, rvec.z);
    glTranslated(-x0.x, -x0.y, -x0.z);

    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
    glPopMatrix();
}

template <class TMesh>
void RigidBodyRenderer::render(const LSCollisionRigidBody<double, TMesh>* rbody,
        const Vector3<double>& x, const Quaternion<double>& rot) const
{
    using namespace std;

    const TriangleMesh<double>*    pmesh = rbody->collision_processor()->surface_mesh();    // surface triangles
    const vector<Tuple3ui>&     tgl = pmesh->surface_indices();                             // surface triangles
    const vector< Point3<double> >&  vtx = pmesh->vertices();                               // surface triangles
    const vector< Vector3<double> >& nml = rbody->collision_processor()->vtx_pseudo_nml();  // surface triangles
    const Point3<double>&     x0 = rbody->initial_mass_center();
    Vector3<double> rvec;
    REAL rd = rot.toAxisRotD(rvec);
    glEnable(GL_LIGHTING);

    glColor3fv(color_);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)&vtx[0]);
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)&nml[0]);

    glPushMatrix();
    glTranslated(x.x, x.y, x.z);
    glRotated(rd, rvec.x, rvec.y, rvec.z);
    glTranslated(-x0.x, -x0.y, -x0.z);

    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
    glPopMatrix();
}

#endif

