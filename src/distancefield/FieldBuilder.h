//////////////////////////////////////////////////////////////////////
// FieldBuilder.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef FIELD_BUILDER_H
#define FIELD_BUILDER_H

#include <TYPES.h>

#include "distanceField.h"
#include "closestPointField.h"

#include <geometry/TriangleMesh.hpp>

#include <math.h>
#include <string>
#include <vector>

//////////////////////////////////////////////////////////////////////
// Distance field construction functions
//////////////////////////////////////////////////////////////////////
class DistanceFieldBuilder {
    public:
        // Builds a signed distance function out of the given .obj file
        // Optionally attempts to read/write from/to the given file name
        //
        //  File name is of the form %s_%04d.sdf where %s = filePrefix and %04d = resolution
        static DistanceField    *BuildSignedDistanceField( const std::string &objFileName,
                                                           int resolution,
                                                           const std::string &filePrefix );

        static ClosestPointField    *BuildSignedClosestPointField( const std::string &objFileName,
	                                                           int resolution,
	                                                           const std::string &filePrefix );
        
};

#endif
