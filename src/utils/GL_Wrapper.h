#ifndef GL_WRAPPER_H 
#define GL_WRAPPER_H 

namespace GL_Wrapper
{
void DrawSphere(double r, int lats, int longs) {
    int i, j;
    for(i = 0; i <= lats; i++) {
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
};  // namespace GLWrapper

#endif 
