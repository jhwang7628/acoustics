//////////////////////////////////////////////////////////////////////
// FieldBuilder.cpp
//////////////////////////////////////////////////////////////////////

#include "FieldBuilder.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Builds a signed distance function out of the given triangle mesh
// Optionally attempts to read/write from/to the given file name
//
//  File name is of the form %s_%04d.sdf where %s = filePrefix and %04d = resolution
//////////////////////////////////////////////////////////////////////
DistanceField *
DistanceFieldBuilder::BuildSignedDistanceField( const string &objFileName,
                                                int resolution,
                                                const string &filePrefix )
{
    char                     filename[ 1024 ];
    DistanceField           *field = NULL;
    int                      status = 1;

    field = new DistanceField();

    if ( !filePrefix.empty() ) {
        sprintf( filename, "%s_%04d.sdf", filePrefix.c_str(), resolution );

        status = field->load( string(filename) );
    }

    if ( status != 0 ) {
        // Couldn't read the field, so we need to build it explicitly
        field->computeSignedField( objFileName, resolution );

        if ( !filePrefix.empty() ) {
            field->save( string( filename ), false );
        }
    }

    return field;
}
