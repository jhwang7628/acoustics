#ifndef _TESTING_H
#define _TESTING_H

#ifndef SDUMP
#define SDUMP(x)	" " << #x << "=[ " << x << " ] "
#endif 

#ifndef COUT_SDUMP
#define COUT_SDUMP(x) std::cout << SDUMP(x) << std::endl
#endif


#endif 
