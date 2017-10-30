#ifndef GL_WRAPPER_H 
#define GL_WRAPPER_H 

namespace GL_Wrapper
{
void DrawSphere(double r, int lats, int longs); 
void DrawWireBox(const double *const minBound, const double *const maxBound);
void DrawBox(const double *const minBound, const double *const maxBound);
};  // namespace GLWrapper

#endif 
