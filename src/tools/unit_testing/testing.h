#ifndef _TESTING_H
#define _TESTING_H

#include <iostream> 
#include <utils/SimpleTimer.h>

#ifndef SDUMP
#define SDUMP(x)	" " << #x << "=[ " << x << " ] "
#endif 

#ifndef COUT_SDUMP
#define COUT_SDUMP(x) std::cout << SDUMP(x) << std::endl
#endif

#ifndef COUT_SDUMP_MAT
#define COUT_SDUMP_MAT(x) std::cout << #x << "=\n" << x << " " << std::endl
#endif

#endif 
