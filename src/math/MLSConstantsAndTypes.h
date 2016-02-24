#ifndef _MLS_CONSTANTS_AND_TYPES_H
#define _MLS_CONSTANTS_AND_TYPES_H

enum MLS_lookup_t {FLEISHMAN, COMPACT};

//static const int MLS_ATTRIB_DIM = 3;
//const double MLS_FLEISHMAN_RAD = 0.005;
static const double MLS_FLEISHMAN_RAD = 0.018;
static const double MLS_COMPACT_RAD = 0.025;
static const MLS_lookup_t MLS_LOOKUP_TYPE = FLEISHMAN;
static const int MLS_LOOKUP_ORDER = 1; // for now, 0/1/2

#endif // _MLS_CONSTANTS_AND_TYPES_H

