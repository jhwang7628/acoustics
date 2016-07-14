#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <geometry/Point3.hpp>
#include <ui/ModalViewer.h>
#include <utils/STL_Wrapper.h>

using namespace qglviewer;
using namespace std;

//##############################################################################
//##############################################################################
void ModalViewer::init()
{
    restoreStateFromFile();
    glDisable(GL_LIGHTING);
    glPointSize(3.0);
    setGridIsDrawn();
    help();
    startAnimation();
}

//##############################################################################
//##############################################################################
void ModalViewer::draw()
{
    // draw rigid sound object mesh
    std::shared_ptr<TriangleMesh<REAL> > meshPtr = _rigidSoundObject->GetMeshPtr();
    const std::vector<Point3<REAL> >  &vertices = meshPtr->vertices(); 
    const std::vector<Tuple3ui>       &triangles = meshPtr->triangles(); 
    const std::vector<Vector3<REAL> > &normals = meshPtr->normals();  // defined on vertices
    const int N_vertices = vertices.size(); 
    const int N_triangles = triangles.size(); 
    const int N_normals = normals.size(); 
    const REAL offsetEpsilon = 1E-5;

    // plot points
    glBegin(GL_POINTS);
    glColor3f(0.6f, 0.6f, 0.6f); 
    for (int v_idx=0; v_idx<N_vertices; ++v_idx)
    {
        // offset the points to make it more visible
        const Point3<REAL> &vertex = vertices.at(v_idx); 
        const Vector3<REAL> &normal = normals.at(v_idx); 
        Point3<REAL> offsetVertex = vertex + normal.normalized() * offsetEpsilon;
        glVertex3f(offsetVertex.x, offsetVertex.y, offsetVertex.z);
    }
    glEnd();

    // plot triangles
    glBegin(GL_TRIANGLES); 
    glColor3f(1, 1, 1); 
    for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
    {
        const Tuple3ui &triangle = triangles.at(t_idx); 
        const Point3<REAL> &x = vertices.at(triangle.x); 
        const Point3<REAL> &y = vertices.at(triangle.y); 
        const Point3<REAL> &z = vertices.at(triangle.z); 
        glVertex3f(x.x, x.y, x.z); 
        glVertex3f(y.x, y.y, y.z); 
        glVertex3f(z.x, z.y, z.z); 
    }
    glEnd(); 

    // plot edges of the triangles
    glBegin(GL_LINES); 
    glColor3f(0.6f, 0.6f, 0.6f); 
    for (int t_idx=0; t_idx<N_triangles; ++t_idx) 
    {
        const Tuple3ui &triangle = triangles.at(t_idx); 
        Point3<REAL> x = vertices.at(triangle.x); 
        Point3<REAL> y = vertices.at(triangle.y); 
        Point3<REAL> z = vertices.at(triangle.z); 
        const Vector3<REAL> &normal_x = normals.at(triangle.x); 
        const Vector3<REAL> &normal_y = normals.at(triangle.y); 
        const Vector3<REAL> &normal_z = normals.at(triangle.z); 
        x = x + normal_x.normalized() * offsetEpsilon; 
        y = y + normal_y.normalized() * offsetEpsilon; 
        z = z + normal_z.normalized() * offsetEpsilon; 
        glVertex3f(x.x, x.y, x.z); 
        glVertex3f(y.x, y.y, y.z); 

        glVertex3f(y.x, y.y, y.z); 
        glVertex3f(z.x, z.y, z.z); 

        glVertex3f(z.x, z.y, z.z); 
        glVertex3f(x.x, x.y, x.z); 
    }
    glEnd(); 


}

//##############################################################################
//##############################################################################
void ModalViewer::animate()
{
}

//##############################################################################
//##############################################################################
QString ModalViewer::helpString() const
{
    QString text("<h2>Modal ModalViewer</h2>");
    text += "Used for debugging FDTD_RigidSoundObject and associated classes";
    return text;
}

