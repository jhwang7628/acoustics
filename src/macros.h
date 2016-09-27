#ifndef MACROS_H 
#define MACROS_H 

//##############################################################################
// return true if they are equal in the sense that they are closed enough
//##############################################################################
#ifndef EQUAL_FLOATS
#define EQUAL_FLOATS(x, y) \
    (fabs((x)-(y)) < 1E-12)
#endif

//##############################################################################
// clear cin buffer
// (http://stackoverflow.com/questions/257091/how-do-i-flush-the-cin-buffer)
//##############################################################################
#ifndef PRE_CIN_CLEAR
#define PRE_CIN_CLEAR \
std::ios_base::sync_with_stdio(false); \
const int __bufferSize = (int)std::cin.rdbuf()->in_avail(); \
std::ios_base::sync_with_stdio(true); \
if (__bufferSize > 0) \
    std::cin.ignore(__bufferSize - 1);
#endif

#ifndef POST_CIN_CLEAR
#define POST_CIN_CLEAR \
std::cin.clear(); \
std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
#endif

//##############################################################################
#endif 
