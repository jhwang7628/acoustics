#ifndef MACROS_H 
#define MACROS_H 

// return true if they are equal in the sense that they are closed enough
#ifndef EQUAL_FLOATS
#define EQUAL_FLOATS(x, y) \
    (fabs(x-y) < 1E-12)
#endif

#endif 
