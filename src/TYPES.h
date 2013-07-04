// Types.h: Project-wide type definitions
//
//////////////////////////////////////////////////////////////////////

#ifndef TYPES_H
#define TYPES_H

#include "config.h"

#include <set>
#include <vector>

#if 0
#include <linearalgebra/VEC3.h>
#endif

using namespace std;

#include <cmath>

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>

#include <map>

typedef std::vector<int>			        IntArray;
typedef std::vector<bool>			        BoolArray;
#if 0
typedef std::vector<VEC3F>                  Vector3Array;
#endif
typedef std::vector<REAL>			        FloatArray;

#if 0
typedef TVEC3<int>					        Tuple3i;
typedef std::vector<Tuple3i>	            Tuple3Array;
#endif

typedef std::pair<int,int>                  IndexRange;
typedef std::pair<int,int>                  IntPair;
typedef std::pair<int,REAL>                 KeyValuePair;
typedef std::vector<KeyValuePair>           KeyValueArray;
typedef std::pair<IndexRange, IndexRange>   RangePair;
typedef std::set<IntPair>                   IndexPairSet;

// 2D and 3D arrays for various data types
typedef boost::multi_array<REAL, 2>		    ScalarArray2D;
typedef boost::multi_array<REAL, 3>		    ScalarArray3D;
#if 0
typedef boost::multi_array<VEC3F, 2>	    Vec3Array2D;
typedef boost::multi_array<VEC3F, 3>	    Vec3Array3D;
#endif
typedef boost::multi_array<bool, 2>		    BoolArray2D;
typedef boost::multi_array<bool, 3>		    BoolArray3D;
typedef boost::multi_array<int, 2>		    IntArray2D;
typedef boost::multi_array<int, 3>		    IntArray3D;

template <class T, class R>
struct UnorderedMap {
    typedef boost::unordered_map<T, R> type;
};

static const int NUM_OCTANTS = 8;

enum AXIS {
    X_AXIS = 0,
    Y_AXIS,
    Z_AXIS
};

enum ACCELERATION_DIRECTION {
    TRANS_X = 0,
    TRANS_Y,
    TRANS_Z,
    ROT_X,
    ROT_Y,
    ROT_Z,
    NUM_ACCEL_DIRECTIONS
};

enum BoundaryCondition {
    BC_NONE = 0,
    BC_DIRICHLET,
    BC_NEUMANN
};

enum GridDirection {
    GRID_X = 0,
    GRID_Y,
    GRID_Z,
    NUM_DIRECTIONS
};

enum ObjectMeasure {
    SURFACE_AREA = 0,
    VOLUME,
    MATCH_RIGID
};

#endif
