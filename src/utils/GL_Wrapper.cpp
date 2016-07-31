#include <math.h>
#include <utils/GL_Wrapper.h>
#include <QGLViewer/qglviewer.h>

namespace GL_Wrapper 
{
//##############################################################################
//##############################################################################
void DrawSphere(double r, int lats, int longs) 
{
    int i, j;
    for(i = 0; i <= lats; i++) 
    {
        double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
        double z0  = sin(lat0);
        double zr0 =  cos(lat0);

        double lat1 = M_PI * (-0.5 + (double) i / lats);
        double z1 = sin(lat1);
        double zr1 = cos(lat1);

        glBegin(GL_QUAD_STRIP);
        for(j = 0; j <= longs; j++) {
            double lng = 2 * M_PI * (double) (j - 1) / longs;
            double x = cos(lng);
            double y = sin(lng);
            glNormal3f(r* x * zr0, r* y * zr0, r*z0);
            glVertex3f(r* x * zr0, r* y * zr0, r*z0);
            glNormal3f(r* x * zr1, r* y * zr1, r*z1);
            glVertex3f(r* x * zr1, r* y * zr1, r*z1);
        }
        glEnd();
    }
}

//##############################################################################
//##############################################################################
void DrawWireBox(const double *const minBound, const double *const maxBound) 
{
    glBegin(GL_LINE_LOOP); 
    glVertex3f(minBound[0], minBound[1], minBound[2]); 
    glVertex3f(minBound[0], maxBound[1], minBound[2]); 
    glVertex3f(maxBound[0], maxBound[1], minBound[2]); 
    glVertex3f(maxBound[0], minBound[1], minBound[2]); 
    glEnd();

    glBegin(GL_LINE_LOOP); 
    glVertex3f(minBound[0], minBound[1], maxBound[2]); 
    glVertex3f(minBound[0], maxBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], maxBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], minBound[1], maxBound[2]); 
    glEnd();

    glBegin(GL_LINE_LOOP); 
    glVertex3f(minBound[0], minBound[1], minBound[2]); 
    glVertex3f(minBound[0], minBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], minBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], minBound[1], minBound[2]); 
    glEnd();

    glBegin(GL_LINE_LOOP); 
    glVertex3f(minBound[0], maxBound[1], minBound[2]); 
    glVertex3f(minBound[0], maxBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], maxBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], maxBound[1], minBound[2]); 
    glEnd();

    glBegin(GL_LINE_LOOP); 
    glVertex3f(minBound[0], minBound[1], minBound[2]); 
    glVertex3f(minBound[0], maxBound[1], minBound[2]); 
    glVertex3f(minBound[0], maxBound[1], maxBound[2]); 
    glVertex3f(minBound[0], minBound[1], maxBound[2]); 
    glEnd();

    glBegin(GL_LINE_LOOP); 
    glVertex3f(maxBound[0], minBound[1], minBound[2]); 
    glVertex3f(maxBound[0], maxBound[1], minBound[2]); 
    glVertex3f(maxBound[0], maxBound[1], maxBound[2]); 
    glVertex3f(maxBound[0], minBound[1], maxBound[2]); 
    glEnd();
}
}; // namespace GL_Wrapper
