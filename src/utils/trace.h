#ifndef _TRACE_H_
#define _TRACE_H_

#include <cstdlib>
#include <stdio.h>
#ifdef WIN32
#include <conio.h>
#endif

#define DO_TRACE 1
#undef DO_TRACE

#ifdef DO_TRACE
//----------------------------------------
#define TRACE(priority,...) { \
	if (priority < 2) \
	{ \
		fprintf (stderr, "[ %s,%s,%d ] ", __FILE__, __FUNCTION__, __LINE__); \
		fprintf (stderr, __VA_ARGS__); \
		fprintf (stderr, "\n"); \
	} \
}

//----------------------------------------
#define TRACE_ERROR(...) { \
	fprintf (stderr, "[ ERROR @ %s,%s,%d ] ", __FILE__, __FUNCTION__, __LINE__); \
	fprintf (stderr, __VA_ARGS__); \
	fprintf (stderr, "\n"); \
	fprintf (stderr, "Quit? (1/0) >> "); \
	int quit; \
	scanf( "%d", &quit ); \
	if( quit != 0 ) exit (1); \
}

//----------------------------------------
#define TRACE_WARNING(...) { \
	fprintf (stderr, "[ WARNING ] %s: ", __FUNCTION__); \
	fprintf (stderr, __VA_ARGS__); \
	fprintf (stderr, "\n"); \
}

//----------------------------------------
#define TRACE_ASSERT(pred, ...)	\
	if( !(pred) )	\
	{	\
		TRACE_ERROR( "( " #pred " ) " __VA_ARGS__ );	\
	}
#else // DO_TRACE not defined
#define TRACE(priority,...) { \
}

//----------------------------------------
#define TRACE_ERROR(...) { \
}

//----------------------------------------
#define TRACE_WARNING(...) { \
}

//----------------------------------------
#define TRACE_ASSERT(pred, ...)	{ \
}
#endif // DO_TRACE

#endif

