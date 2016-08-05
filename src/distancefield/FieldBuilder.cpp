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

ClosestPointField *
DistanceFieldBuilder::BuildSignedClosestPointField( const std::string &objFileName,
                              int resolution,
                              const std::string &filePrefix )
{
    char                     filename[ 1024 ];
    ClosestPointField        *field = NULL;
    int                      status = 1;

    field = new ClosestPointField();

    if ( !filePrefix.empty() ) {
        sprintf( filename, "%s_%04d.scpf", filePrefix.c_str(), resolution );

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
ClosestPointField *
DistanceFieldBuilder::BuildSignedClosestPointField( const std::string &objFileName,
                                                int resolution,
                                                const std::string &filePrefix,
                                                const Vec3d &boundingMin,
                                                const Vec3d &boundingMax)
{
    char                     filename[ 1024 ];
    ClosestPointField        *field = NULL;
    int                      status = 1;

    field = new ClosestPointField();

    if ( !filePrefix.empty() ) {
        sprintf( filename, "%s_%04d.scpf", filePrefix.c_str(), resolution );

        status = field->load( string(filename) );
    }

    if ( status != 0 ) {
        // Couldn't read the field, so we need to build it explicitly
        field->setBoundingBox(boundingMin, boundingMax); 
        field->computeSignedField( objFileName, resolution );

        if ( !filePrefix.empty() ) {
            field->save( string( filename ), false );
        }
    }

    return field;
}

AdaptiveDistanceField*
DistanceFieldBuilder::BuildAdaptiveDistanceField(const std::string& objFileName, const std::string &filePrefix,
                                                 double subdivideRadius, int maxAdfLevel, double maxError)
{
    char filename[1024];
    int status = 1;
    
    AdaptiveDistanceField* field = new AdaptiveDistanceField();
    if (!filePrefix.empty())
    {
        sprintf(filename, "%s_sr%g_ml%02d_me%g.adf", filePrefix.c_str(), subdivideRadius, maxAdfLevel, maxError);
        status = field->load(string(filename));
    }
    if (status != 0)
    {
        // build the field explicitly
        field->computeSignedField(objFileName, subdivideRadius, maxAdfLevel, maxError);
        if (!filePrefix.empty())
            field->save(string(filename));
    }
    return field;
}
