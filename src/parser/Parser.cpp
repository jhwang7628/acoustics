// Parser.cpp: Definition for theeParser class
//
//////////////////////////////////////////////////////////////////////

#include "Parser.h"

#include <io/TglMeshReader.hpp>

#include <utils/IO.h>

#if 0
#include <file/obj.h>
#endif

#include <transfer/CompressedMultiTermApproximation.h>
#include <transfer/MultiTermApproximation.h>
#include <transfer/PulseApproximation.h>

#include <iostream>
#include <vector>

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Parser *Parser::buildParser(string filename)
{
    // Construct a new tinyxml document and try to load the
    // given file
    TiXmlDocument *document = new TiXmlDocument();

    if (document->LoadFile(filename.c_str()))
    {
        Parser* parser = new Parser(document);
        parser->filename = filename;
        return parser;
    }
    else
    {
        return NULL;
    }
}

//////////////////////////////////////////////////////////////////////
// Construct a new parser
//////////////////////////////////////////////////////////////////////
Parser::Parser(TiXmlDocument *document)
    : document(document)
{
}

//////////////////////////////////////////////////////////////////////
// Clean up
//////////////////////////////////////////////////////////////////////
Parser::~Parser()
{
    delete document;
}

TriangleMesh<REAL> *Parser::getMesh()
{
    throw std::runtime_error("**ERROR** need to pass in root tag for xml parsing. for example: impact");
    return nullptr;
}
//////////////////////////////////////////////////////////////////////
// Construct a triangle mesh from data stored in the XML document
//////////////////////////////////////////////////////////////////////
TriangleMesh<REAL> *Parser::getMesh(const char * rootTag)
{
    const char              *fileName;
    REAL                     scale = 1.0;
    TriangleMesh<REAL>      *mesh = NULL;

    TiXmlElement *root = document->FirstChildElement(rootTag);
    if (!root)
    {
        cerr << "No " << rootTag << " tag in input file" << endl;
        return NULL;
    }

    // Read in the mesh element, which should provide the
    // mesh file name
    TiXmlElement *meshNode = root->FirstChildElement("mesh");
    if (!meshNode)
    {
        cerr << "No mesh tag in input file" << endl;
        return NULL;
    }

    // Try to get a file name from the mesh tag
    fileName = meshNode->Attribute("file");

    if (fileName == NULL)
    {
        cerr << "No mesh file name" << endl;
        return NULL;
    }

    // We also get a scale value from this tag (but don't
    // require it)
    if (meshNode->QueryRealAttribute("scale", &scale) != TIXML_SUCCESS)
    {
        // Assume no additional scaling
        scale = 1.0;
    }

    // Construct the mesh and attempt to read it from the given file
    mesh = new TriangleMesh<REAL>();

    if ( MeshObjReader::read( fileName, *mesh, false, false, scale ) == SUCC_RETURN ) {
        // Prepare the mesh for rendering, etc.
        mesh->generate_normals();

        return mesh;
    } else {
        // Unsuccessful
        delete mesh;

        return NULL;
    }
}

//////////////////////////////////////////////////////////////////////
// Mesh file name
//////////////////////////////////////////////////////////////////////
std::string Parser::getMeshFileName()
{
    const char            *fileName;

    TiXmlElement *root = document->FirstChildElement("impact");
    if (!root)
        root = document->FirstChildElement("impulse_response");
    if (!root)
    {
        cerr << "No impact tag in input file" << endl;
        return NULL;
    }

    // Read in the mesh element, which should provide the
    // mesh file name
    TiXmlElement *meshNode = root->FirstChildElement("mesh");
    if (!meshNode)
    {
        cerr << "No mesh tag in input file" << endl;
        return NULL;
    }

    // Try to get a file name from the mesh tag
    fileName = meshNode->Attribute("file");

    return string( fileName );
}

#if 0
//////////////////////////////////////////////////////////////////////
// Fetches equivalent source parametersf from the XML document
//////////////////////////////////////////////////////////////////////
Parser::EquivalentSourceParms Parser::getEquivalentSourceParms()
{
    string                   accelerationName;
    REAL                     accelerationPulseWidth;
    EquivalentSourceParms    sourceParms;

    accelerationName = queryRequiredAttr( "impact/acceleration", "pulsetype" );

    accelerationPulseWidth = queryRequiredReal( "impact/acceleration",
            "pulsewidth" );

    sourceParms._acceleration
        = InterpolationFunction::GetFunctionContainer(
                InterpolationFunction::GetFunctionName( accelerationName ),
                accelerationPulseWidth );

    sourceParms._multipoleTerms = queryRequiredInt( "impact/acceleration",
            "multipoleterms" );

    sourceParms._timeScale = queryRequiredReal( "impact/acceleration",
            "pulsetimescale" );

    sourceParms._pulseInterval = queryOptionalReal( "impact/acceleration",
            "pulseinterval", "-1.0" );

    sourceParms._timeStep = queryRequiredReal( "impact/acceleration",
            "integrationtimestep" );

    sourceParms._c = queryOptionalReal( "impact/acceleration", "c", "343.0" );
    sourceParms._density = queryOptionalReal( "impact/acceleration", "density",
            "1.0" );

    return sourceParms;
}
#endif

#if 0
//////////////////////////////////////////////////////////////////////
// Fetches parameters for source position training
//////////////////////////////////////////////////////////////////////
Parser::SourceTrainingParms Parser::getSourceTrainingParms()
{
    string                   accelerationName;
    REAL                     accelerationPulseWidth;
    SourceTrainingParms      trainingParms;

    accelerationName = queryRequiredAttr( "impact/training", "pulsetype" );

    accelerationPulseWidth = queryRequiredReal( "impact/training",
            "pulsewidth" );

    trainingParms._acceleration
        = InterpolationFunction::GetFunctionContainer(
                InterpolationFunction::GetFunctionName( accelerationName ),
                accelerationPulseWidth );

    trainingParms._timeScale = queryRequiredReal( "impact/training",
            "pulsetimescale" );

    trainingParms._timeStep = queryRequiredReal( "impact/training",
            "integrationtimestep" );

    trainingParms._c = queryOptionalReal( "impact/training", "c", "343.0" );
    trainingParms._density = queryOptionalReal( "impact/training", "density",
            "1.0" );

    // Get SDF parameters
    trainingParms._sdfResolution = queryRequiredInt( "impact/training",
            "fieldresolution" );

    trainingParms._sdfFilePrefix = queryOptionalAttr( "impact/training",
            "distancefield",
            "distancefield" );

    trainingParms._pointsPerVoxel = queryRequiredInt( "impact/training",
            "pointspervoxel" );

    // Get training parms
    trainingParms._sourcesPerIteration
        = queryRequiredInt( "impact/training",
                "sourcesperiteration" );

    trainingParms._candidatesPerIteration
        = queryRequiredInt( "impact/training",
                "candidatesperiteration" );

    trainingParms._maxSourcesPerBatch
        = queryRequiredInt( "impact/training",
                "maxsourcesperbatch" );

    trainingParms._errorTolerance = queryRequiredReal( "impact/training",
            "tolerance" );

    return trainingParms;
}
#endif

#if 0
//////////////////////////////////////////////////////////////////////
// Fetches parameters for faking a sphere-like collision (with an
// arbitrary mesh)
//////////////////////////////////////////////////////////////////////
Parser::SphericalImpactParms Parser::getSphereImpactParms()
{
    SphericalImpactParms       sphereParms;
    string                     listeningPosition;
    string                     outputFile;

    sphereParms._c = queryOptionalReal( "impact/sphereimpact", "c", "343.0" );
    sphereParms._density = queryOptionalReal( "impact/sphereimpact", "density",
            "1.0" );

    // Get SDF parameters
    sphereParms._sdfResolution = queryRequiredInt( "impact/sphereimpact",
            "fieldresolution" );

    sphereParms._sdfFilePrefix = queryOptionalAttr( "impact/sphereimpact",
            "distancefield",
            "distancefield" );

    // Get grid parameters
    sphereParms._gridResolution = queryRequiredInt( "impact/sphereimpact",
            "gridresolution" );

    sphereParms._gridScale = queryRequiredReal( "impact/sphereimpact",
            "gridscale" );

    sphereParms._timeStepFrequency = queryOptionalReal( "impact/sphereimpact",
            "timestepfrequency", "96000" );

    // Get boundary information
    sphereParms._sphereRadius = queryRequiredReal( "impact/sphereimpact",
            "sphereradius" );

    sphereParms._impactSpeed = queryRequiredReal( "impact/sphereimpact",
            "impactspeed" );

    // Material parameters
    sphereParms._boundary_E = queryRequiredReal( "impact/sphereimpact",
            "boundary_E" );
    sphereParms._boundary_nu = queryRequiredReal( "impact/sphereimpact",
            "boundary_nu" );
    sphereParms._boundaryDensity = queryRequiredReal( "impact/sphereimpact",
            "boundarydensity" );

    // Get sound output parameters
    listeningPosition = queryRequiredAttr( "impact/sphereimpact",
            "listeningposition" );
    sphereParms._listeningPosition = str2vec3d( listeningPosition );

    sphereParms._listeningPositions.clear();
    getListeningPositions( "sphereimpact", sphereParms._listeningPositions );

    outputFile = queryRequiredAttr( "impact/sphereimpact", "outputfile" );
    sphereParms._outputFile = outputFile;

    return sphereParms;
}
#endif

std::vector<Parser::ImpulseResponseParms::VolumetricSource> Parser::ImpulseResponseParms::QueryVolumetricSource(TiXmlNode *document, Parser *parser, const std::string &path, const REAL &soundSpeed)
{
    if (!document) throw std::runtime_error("**ERROR** document not established"); 
    
    std::vector<VolumetricSource> sources; 
    
    TiXmlNode *node = document; 
    node = parser->getNodeByPath(path); 
    while (node) 
    {
        Parser::ImpulseResponseParms::VolumetricSource source; 

        TiXmlElement *elm = node->ToElement(); 
        source.widthTime         = parser->queryRequiredReal(elm, "source_width_time"); 
        source.widthSpace        = parser->queryOptionalReal(elm, "source_width_space", soundSpeed*source.widthTime); 
        source.offsetTime        = parser->queryOptionalReal(elm, "source_offset_time", 0.0); 
        source.normalizeConstant = parser->queryOptionalReal(elm, "source_normalize_constant", 1.0/pow(sqrt(2.0*M_PI)*source.widthSpace,3)); 
        //source.widthSpace = c*widthTime; 
        //source.offsetTime = 0.0; // default 
        source.position.x = parser->queryRequiredReal(elm, "source_position_x"); 
        source.position.y = parser->queryRequiredReal(elm, "source_position_y"); 
        source.position.z = parser->queryRequiredReal(elm, "source_position_z"); 

        std::cout << SDUMP(source.position) << std::endl;

        source.flipSign = (parser->queryOptionalReal(elm, "source_sign_flip", 0.0) > 1E-10) ? true : false; 
        source.normalizeConstant *= (source.flipSign) ? -1.0 : 1.0;

        sources.push_back(source); 

        node = parser->getNextSiblingNode(node); 
    }

    if (sources.size()==0) 
        throw std::runtime_error("**ERROR** no sources found in the configuration file"); 

    return sources;
}

Parser::ImpulseResponseParms Parser::getImpulseResponseParms()
{
    ImpulseResponseParms       parms;

    parms._c = queryOptionalReal( "impulse_response/solver", "c", "343.0" );
    parms._density = queryOptionalReal( "impulse_response/solver", "density", "1.22521" );

    // Get SDF parameters
    parms._sdfResolution = queryRequiredInt( "impulse_response/solver", "fieldresolution" );
    parms._sdfFilePrefix = queryRequiredAttr( "impulse_response/solver", "distancefield" );

    // Get grid parameters
    parms._gridResolution = queryRequiredInt( "impulse_response/solver", "gridresolution" );
    //parms._gridScale = 1.0;
    parms._gridScale = queryOptionalReal( "impulse_response/solver", "gridscale", "1" );

    parms._cellSize = queryOptionalReal( "impulse_response/solver", "cellsize", "-1"); 
    //parms._fixedCellSize = ( queryOptionalInt( "impulse_response/solver", "fixedcellsize", "0" ) != 0 );
    //if ( parms._fixedCellSize )
    //    parms._cellSize = queryRequiredReal( "impulse_response/solver", "cellsize" );

    parms._timeStepFrequency = queryRequiredInt( "impulse_response/solver", "timestepfrequency" );
    parms._subSteps = queryRequiredInt( "impulse_response/solver", "substeps" );

    parms._stopTime = queryRequiredReal( "impulse_response/solver", "stop_time" ); 

    parms._outputPattern = queryRequiredAttr( "impulse_response/solver", "output_pattern" );

    parms._listeningFile = queryOptionalAttr( "impulse_response/solver", "listening_file", "none" );

    parms._useMesh = (queryOptionalInt("impulse_response/solver", "use_mesh", "1")==0) ? false : true; 
    parms._cornellBoxBoundaryCondition = (queryOptionalInt("impulse_response/solver", "cornell_box_boundary_condition", "0")==1) ? true : false; 

    // TODO very ineffective way to parse it but easier to implement for now,
    // can make it better

    parms._sources = ImpulseResponseParms::QueryVolumetricSource(document, this, "impulse_response/volumetric_source/source", parms._c); 


    parms._f = queryOptionalReal( "impulse_response/solver", "f", "500" );
    //std::vector<REAL> sourcesPosition_x = queryRequiredRealMultiple("impulse_response/volumetric_source/source", "source_position_x" ); 
    //std::vector<REAL> sourcesWidthTime  = queryRequiredRealMultiple("impulse_response/volumetric_source/source", "source_width_time" ); 
    //std::vector<REAL> sourcesPosition_y = queryRequiredRealMultiple("impulse_response/volumetric_source/source", "source_position_y" ); 
    //std::vector<REAL> sourcesPosition_z = queryRequiredRealMultiple("impulse_response/volumetric_source/source", "source_position_z" ); 

    //if (sourcesWidthTime.size() != sourcesPosition_x.size() || 
    //    sourcesWidthTime.size() != sourcesPosition_y.size() || 
    //    sourcesWidthTime.size() != sourcesPosition_z.size()) 
    //    throw std::runtime_error("**ERROR** source attribute misaligned. number of specified source_width_time is inconsistent with the number of source_position"); 

    //const int N_sources = sourcesWidthTime.size(); 

    //parms._sources.resize(N_sources); 


    //for (int ii=0; ii<N_sources; ii++) 
    //{
    //    Parser::ImpulseResponseParms::VolumetricSource source; 

    //    source.widthTime = sourcesWidthTime[ii]; 
    //    source.position_x = sourcesPosition_x[ii]; 
    //    source.position_y = sourcesPosition_y[ii]; 
    //    source.position_z = sourcesPosition_z[ii]; 

    //    parms._sources[ii] = source; 


    //}

    //for (int ii=0; ii<N_sources; ii++) 
    //    std::cout << parms._sources[ii] << std::endl; 


    //parms._sourceWidthTime = queryRequiredReal( "impulse_response/solver", "source_width_time" );
    //parms._sourcePosition_x = queryRequiredReal( "impulse_response/solver", "source_position_x" ); 
    //parms._sourcePosition_y = queryRequiredReal( "impulse_response/solver", "source_position_y" ); 
    //parms._sourcePosition_z = queryRequiredReal( "impulse_response/solver", "source_position_z" ); 

    return parms;
}

Parser::AcousticTransferParms Parser::getAcousticTransferParms()
{
    AcousticTransferParms      transferParams;
    string                     modeData;
    string                     outputFile;
    string                     multiTermOutputFile;

    transferParams._c = queryOptionalReal( "impact/acoustic_transfer",
            "c", "343.0" );
    transferParams._density = queryOptionalReal( "impact/acoustic_transfer",
            "density", "1.22521" );

    // Get SDF parameters
    transferParams._sdfResolution = queryRequiredInt( "impact/acoustic_transfer",
            "fieldresolution" );

    transferParams._sdfFilePrefix = queryRequiredAttr( "impact/acoustic_transfer",
            "distancefield" );

    // Get grid parameters
    transferParams._gridResolution = queryRequiredInt( "impact/acoustic_transfer",
            "gridresolution" );

    transferParams._gridScale = queryRequiredReal( "impact/acoustic_transfer",
            "gridscale" );

    transferParams._fixedCellSize
        = ( queryOptionalInt( "impact/acoustic_transfer", "fixedcellsize",
                    "0" ) != 0 );
    if ( transferParams._fixedCellSize )
    {
        transferParams._cellSize = queryRequiredReal( "impact/acoustic_transfer",
                "cellsize" );
    }

    transferParams._timeStepFrequency = queryRequiredInt(
            "impact/acoustic_transfer",
            "timestepfrequency" );
    transferParams._subSteps = queryRequiredInt( "impact/acoustic_transfer",
            "substeps" );
    transferParams._mode = queryOptionalInt( "impact/acoustic_transfer",
            "mode", "0");

    transferParams._nbar = queryRequiredReal( "impact/acoustic_transfer",
            "nbar" );

    transferParams._radiusMultipole = queryRequiredReal( "impact/acoustic_transfer",
            "radius_multipole" );

    modeData = queryRequiredAttr("impact/acoustic_transfer", "mode_data_file");
    transferParams._modeDataFile = modeData;

    transferParams._rigidPrefix = queryRequiredAttr( "impact/acoustic_transfer",
            "rigidprefix" );
    transferParams._rigidDensity = queryRequiredReal( "impact/acoustic_transfer",
            "rigiddensity" );


    outputFile = queryRequiredAttr( "impact/acoustic_transfer", "outputfile" );
    transferParams._outputFile = outputFile;

    transferParams._multipoleOutputFile
        = queryRequiredAttr( "impact/acoustic_transfer", "multipole_outputfile" );

    return transferParams;
}

//////////////////////////////////////////////////////////////////////
// Fetches parameters for acceleration pulse precomputation
//////////////////////////////////////////////////////////////////////
Parser::AccelerationPulseParms Parser::getAccelerationPulseParms()
{
    AccelerationPulseParms     pulseParms;
    string                     listeningPosition;
    string                     outputFile;
    string                     multiTermOutputFile;

    pulseParms._c = queryOptionalReal( "impact/acceleration_pulse",
            "c", "343.0" );
    pulseParms._density = queryOptionalReal( "impact/acceleration_pulse",
            "density", "1.22521" );

    // Get SDF parameters
    pulseParms._sdfResolution = queryRequiredInt( "impact/acceleration_pulse",
            "fieldresolution" );

    pulseParms._sdfFilePrefix = queryRequiredAttr( "impact/acceleration_pulse",
            "distancefield" );

    // Get grid parameters
    pulseParms._gridResolution = queryRequiredInt( "impact/acceleration_pulse",
            "gridresolution" );

    pulseParms._gridScale = queryRequiredReal( "impact/acceleration_pulse",
            "gridscale" );

    pulseParms._fixedCellSize
        = ( queryOptionalInt( "impact/acceleration_pulse", "fixedcellsize",
                    "0" ) != 0 );
    if ( pulseParms._fixedCellSize )
    {
        pulseParms._cellSize = queryRequiredReal( "impact/acceleration_pulse",
                "cellsize" );
    }

    pulseParms._timeStepFrequency = queryRequiredInt(
            "impact/acceleration_pulse",
            "timestepfrequency" );
    pulseParms._subSteps = queryRequiredInt( "impact/acceleration_pulse",
            "substeps" );

    // Get sound output parameters
    pulseParms._pulseTimeScale = queryRequiredReal( "impact/acceleration_pulse",
            "timescale" );

    pulseParms._listeningResolution = queryRequiredInt(
            "impact/acceleration_pulse",
            "listeningresolution" );

    pulseParms._listeningRadius = queryRequiredReal( "impact/acceleration_pulse",
            "listeningradius" );

    pulseParms._numShells = queryRequiredReal( "impact/acceleration_pulse",
            "numshells" );

    pulseParms._numTerms = queryRequiredReal( "impact/acceleration_pulse",
            "numterms" );

    pulseParms._radiusMultiplier = queryRequiredReal( "impact/acceleration_pulse",
            "radiusmultiplier" );

    pulseParms._rigidPrefix = queryRequiredAttr( "impact/acceleration_pulse",
            "rigidprefix" );
    pulseParms._rigidDensity = queryRequiredReal( "impact/acceleration_pulse",
            "rigiddensity" );

    pulseParms._compressionTolerance = queryOptionalReal(
            "impact/acceleration_pulse",
            "compressiontolerance",
            "0.01" );

    pulseParms._compressedOutputFile = queryOptionalAttr(
            "impact/acceleration_pulse",
            "compressedoutputfile",
            "pulsemapdata_comp" );

    pulseParms._symmetrizeField = queryOptionalInt( "impact/acceleration_pulse",
            "symmetrizecompressedfield",
            "0" ) == 1;

    if ( pulseParms._symmetrizeField ) {
        pulseParms._symmetricResolution = queryRequiredInt(
                "impact/acceleration_pulse",
                "symmetricresolution" );
    }

    pulseParms._useDirectionSet = queryOptionalInt( "impact/acceleration_pulse",
            "usedirectionset",
            "0" ) == 1;

    if ( pulseParms._useDirectionSet ) {
        pulseParms._directionSetFile = queryRequiredAttr(
                "impact/acceleration_pulse",
                "directionsetfile" );
    }


    outputFile = queryRequiredAttr( "impact/acceleration_pulse", "outputfile" );
    pulseParms._outputFile = outputFile;

    pulseParms._multiTermOutputFile
        = queryRequiredAttr( "impact/acceleration_pulse", "multitermoutputfile" );

    return pulseParms;
}


#if 0
//////////////////////////////////////////////////////////////////////
// Fetches parameters for setting up the input to an acoustic
// scattering calculation due to rigid body acceleration
//////////////////////////////////////////////////////////////////////
Parser::ScatteringTestParms Parser::getScatteringTestParms()
{
    ScatteringTestParms        scatterParms;

    string                     boundValue;

    scatterParms._c = queryOptionalReal( "impact/scattertest", "c", "343.0" );
    scatterParms._density = queryOptionalReal( "impact/scattertest", "density",
            "1.22521" );

    scatterParms._frequencyFile = queryRequiredAttr( "impact/scattertest",
            "frequencyfile" );

    getScatterObjects( "scattertest", scatterParms._scatterObjects );

    scatterParms._outputPrefix = queryRequiredAttr( "impact/scattertest",
            "outputprefix" );
    scatterParms._spectrumOutputPrefix = queryRequiredAttr( "impact/scattertest",
            "spectrumoutput" );

    scatterParms._responseFile = queryRequiredAttr( "impact/scattertest",
            "responsefile" );

    // Get finite difference grid information (optional)
    boundValue = queryOptionalAttr( "impact/scattertest", "FD_min_bound",
            "0.0 0.0 0.0" );
    scatterParms._FD_min_bound = str2vec3d( boundValue );

    boundValue = queryOptionalAttr( "impact/scattertest", "FD_max_bound",
            "1.0 1.0 1.0" );
    scatterParms._FD_max_bound = str2vec3d( boundValue );

    scatterParms._FD_grid_divisions = queryOptionalInt( "impact/scattertest",
            "FD_grid_divisions",
            "50" );

    scatterParms._FD_timeStepFrequency = queryOptionalInt( "impact/scattertest",
            "FD_timestepfrequency",
            "192000" );
    scatterParms._FD_subSteps = queryOptionalInt( "impact/scattertest",
            "FD_substeps", "5" );

    scatterParms._FD_outputFile = queryOptionalAttr( "impact/scattertest",
            "FD_outputfile",
            "FDTD_pressure" );

    scatterParms._FD_impactTimeScale = queryOptionalReal( "impact/scattertest",
            "FD_impacttimescale",
            "0.0001" );

    getListeningPositions( "scattertest", scatterParms._listeningPositions );

    return scatterParms;
}
#endif

//////////////////////////////////////////////////////////////////////
// Fetches scene objects
//////////////////////////////////////////////////////////////////////
Parser::SceneObjectParms Parser::getSceneObjectParms( bool loadPANData,
                                                      bool compressedFields,
                                                      bool symmetricFields )
{
    SceneObjectParms           parms;
    map<string, int>           meshIDMap;
    vector<string>             rigidFilePrefixes;
    IntArray                   sdfResolutions;
    vector<string>             sdfFileNames;
    vector<string>             pulseModelFileNames;
    FloatArray                 densities;

    bool                       smoothCurvatures = false;

    getMeshes( "scene", parms._meshes, rigidFilePrefixes,
               sdfFileNames, sdfResolutions, pulseModelFileNames,
               densities, meshIDMap, compressedFields );

    parms._matchWaveletSamplingRates = queryOptionalInt( "impact/scene",
                                                         "matchsamplingrate",
                                                         "0" ) == 1;

    smoothCurvatures = queryOptionalInt( "impact/scene",
                                         "smoothcurvatures",
                                         "0" ) == 1;

    // Build SDF meshes in parallel
    parms._distanceMeshes.resize( parms._meshes.size() );
#pragma omp parallel for schedule(static) default(shared)
    for ( size_t mesh_idx = 0; mesh_idx < parms._meshes.size(); mesh_idx++ ) {
        TriangleMesh<REAL>  *mesh = parms._meshes[ mesh_idx ];

#if 0
        parms._distanceMeshes[ mesh_idx ]
            = new ClosestPointMesh( *mesh, sdfResolutions[ mesh_idx ],
                    sdfFileNames[ mesh_idx ] );
#endif
        parms._distanceMeshes[ mesh_idx ] = new ClosestPointMesh( *mesh );
    }

    // Build rigid body, curvature and SDF information
    for ( size_t mesh_idx = 0; mesh_idx < parms._meshes.size(); mesh_idx++ )
    {
        TriangleMesh<REAL>  *mesh = parms._meshes[ mesh_idx ];

        parms._curvatureMeshes.push_back( new GTS_TriMesh( *mesh ) );
        //parms._curvatureMeshes.back()->precomputeMeanCurvatures();
        parms._curvatureMeshes.back()->precomputeMeanCurvatures( smoothCurvatures );

        parms._rigidMeshes.push_back(
                new RigidMesh( *mesh, rigidFilePrefixes[ mesh_idx ],
                    densities[ mesh_idx ] ) );

        if ( loadPANData && pulseModelFileNames[ mesh_idx ] != "null" )
        {
            parms._pulseModels.push_back( new RadialApproximation::AccelerationSet() );

            if ( compressedFields ) {
                CompressedMultiTermApproximation::ReadAccelerationSet(
                        pulseModelFileNames[ mesh_idx ],
                        *parms._pulseModels[ mesh_idx ],
                        symmetricFields,
                        parms._matchWaveletSamplingRates,
                        // No direction set suffix
                        NULL );
            }
            else {
                MultiTermApproximation::ReadAccelerationSet(
                        pulseModelFileNames[ mesh_idx ],
                        *parms._pulseModels[ mesh_idx ] );
            }

            parms._compactPulseModels.push_back(
                    parms._pulseModels[ mesh_idx ]->_allFields[ 0 ]
                    ->buildCompactApproximation(
                        *parms._pulseModels[ mesh_idx ] ) );
        }
        else
        {
            cout << "Setting null approximation in slot " << parms._pulseModels.size() << endl;
            parms._pulseModels.push_back( NULL );
            parms._compactPulseModels.push_back( NULL );
        }
    }

    getObjectMeshIDs( "scene", meshIDMap, parms._objectMeshIDs );

    parms._useProxies = queryOptionalInt( "impact/scene", "useproxies", "0" ) == 1;

    parms._randomizeImpulseTimes = queryOptionalInt( "impact/scene",
            "randomizeimpulsetimes",
            "0" ) == 1;

    if ( parms._randomizeImpulseTimes ) {
        parms._simTimeStep = queryRequiredReal( "impact/scene", "simtimestep" );
    }

    parms._collisionTimeScale = queryOptionalReal( "impact/scene", "collisiontimescale", "1.0" );

    if ( loadPANData && parms._useProxies ) {
        int                      measure = queryOptionalInt( "impact/scene", "proxymeasure", "0" );

        cout << SDUMP( measure ) << endl;

        if ( measure == SURFACE_AREA ) {
            parms._proxyMeasure = SURFACE_AREA;
        } else if ( measure == VOLUME ) {
            parms._proxyMeasure = VOLUME;
        } else {
            parms._proxyMeasure = MATCH_RIGID;
        }

        parms._proxyPathPrefix = queryRequiredAttr( "impact/scene", "proxypathprefix" );
        parms._proxyFilePrefix = queryRequiredAttr( "impact/scene", "proxyfileprefix" );
        parms._proxyMinScale = queryRequiredReal( "impact/scene", "proxyminscale" );
        parms._proxyIncrements = queryRequiredInt( "impact/scene", "proxyincrements" );

        parms._useSymmetricProxies
            = queryOptionalInt( "impact/scene", "symmetricproxies", "0" ) == 1;

        if ( parms._useSymmetricProxies ) {
            parms._symmetricProxyPrefix = queryRequiredAttr( "impact/scene",
                                                             "symmetricproxyprefix" );
        }

        parms._useDirectionSetProxies
            = queryOptionalInt( "impact/scene", "usedirectionsetproxies", "0" ) == 1;

        if ( parms._useDirectionSetProxies ) {
            parms._proxyDirectionSetFile = queryRequiredAttr(
                    "impact/scene", "proxydirectionsetfile" );
            parms._proxyDirectionSetSuffix = queryRequiredAttr(
                    "impact/scene", "proxydirectionsetsuffix" );
        }
    }

    return parms;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
string Parser::queryRequiredAttr( TiXmlElement *elem, std::string attr )
{
    const char* value = elem->Attribute(attr.c_str());

    if( value == NULL )
    {
        cerr << "[ERROR] Required attribute \"" << attr;
        cerr << "\" was not found in the XML elem: " << elem->Value() << endl;
        cerr << "Enter any character to continue, but ";
        cerr << "you should probably abort." << endl;
        char dummy;
        cin >> dummy;
        return string("ERROR");
    }
    else
    {
        return string( value );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
string Parser::queryOptionalAttr( TiXmlElement *node, std::string attr,
        std::string defaultValue )
{
    const char* value = node->Attribute(attr.c_str());

    if( value == NULL )
    {
        cout << "[STATUS] Optional attribute \"" << attr;
        cout << "\" not found. Using default: " << defaultValue << endl;
        return defaultValue;
    }
    else
    {
        return string( value );
    }
}

REAL Parser::queryOptionalReal( TiXmlElement *node, const std::string &attr, const REAL &defaultValue )
{
    const char* value = node->Attribute(attr.c_str());

    if( value == NULL )
    {
        cout << "[STATUS] Optional attribute \"" << attr;
        cout << "\" not found. Using default: " << defaultValue << endl;
        return defaultValue;
    }
    else
    {
        return atof( value );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
TiXmlNode* Parser::getNodeByPath( std::string pathStr )
{
    vector<string> path = UTILS_IO::split( pathStr, string("/") );
    TiXmlNode* node = document;
    for( size_t i = 0; i < path.size(); i++ )
    {
        node = node->FirstChildElement( path[i].c_str() );
        if( node == NULL )
        {
            cerr << "[ERROR] The required XML path did not exist: " << flush << endl;
            for( size_t j = 0; j < path.size(); j++ )
            {
                cerr << j << " " << path[j] << endl;
                return NULL;
            }
        }
    }

    return node;
}

TiXmlNode* Parser::getNextSiblingNode( TiXmlNode* node ) 
{
    return node->NextSibling(); 
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
string Parser::queryOptionalAttr( std::string path, std::string attr,
        std::string defaultVal )
{
    TiXmlNode* node = getNodeByPath( path );
    TiXmlElement* elm = node->ToElement();
    string val = queryOptionalAttr( elm, attr, defaultVal );
    cout << "Queried " << attr << "=\"" << val << "\"" << endl;
    return val;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
string Parser::queryRequiredAttr( std::string path, std::string attr )
{
    TiXmlNode* node = getNodeByPath( path );
    TiXmlElement* elm = node->ToElement();
    string val = queryRequiredAttr( elm, attr );
    cout << "Queried " << attr << "=\"" << val << "\"" << endl;
    return val;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
std::vector<std::string> Parser::queryRequiredAttrMultiple( std::string path, std::string attr )
{
    std::vector<std::string> allval; 
    TiXmlNode* node = getNodeByPath( path );
    if (node != NULL) 
    {
        TiXmlElement* elm = node->ToElement();
        string val = queryRequiredAttr( elm, attr );
        cout << "Queried " << attr << "=\"" << val << "\"" << endl;
        allval.push_back(val); 

        while(true) 
        {
            TiXmlNode* siblingNode = getNextSiblingNode(node); 
            if (siblingNode!=NULL) 
            {
                TiXmlElement* elm = siblingNode->ToElement();
                string val = queryRequiredAttr( elm, attr );
                cout << "Queried " << attr << "=\"" << val << "\"" << endl;
                allval.push_back(val); 

                delete node;  // clear data in node

                node = siblingNode; 
            }
            else
            {
                break; 
            }
        }
    }

    return allval;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
std::vector<REAL> Parser::queryRequiredRealMultiple( const std::string &path, const std::string &attr )
{

    std::vector<std::string> allval = queryRequiredAttrMultiple(path, attr); 
    std::vector<REAL> allval_r(allval.size()); 

    for (size_t ii=0; ii<allval.size(); ii++)
    {
        allval_r[ii] = atof(allval[ii].c_str()); 
    }

    return allval_r;
}



//////////////////////////////////////////////////////////////////////
// Gets a list of listening positions from the 
//////////////////////////////////////////////////////////////////////
void Parser::getListeningPositions( const char *inputElement,
                                    Vector3Array &positions )
{
    TiXmlElement              *root;
    TiXmlElement              *sphereRoot;
    TiXmlElement              *listNode;
    TiXmlElement              *positionNode;

    root = document->FirstChildElement( "impact" );
    if ( !root )
    {
        cerr << "Error: No impact node found" << endl;
        abort();
    }

    sphereRoot = root->FirstChildElement( inputElement );
    if ( !sphereRoot )
    {
        cerr << "Error: No sphereImpact node found" << endl;
        abort();
    }

    listNode = sphereRoot->FirstChildElement( "listeningpositions" );
    if ( !listNode )
    {
        cerr << "Error: No listening position list found" << endl;
        abort();
    }

    positionNode = listNode->FirstChildElement( "position" );

    while ( positionNode != NULL )
    {
        string                   positionValue;

        positionValue = queryRequiredAttr( positionNode, "value" );

        positions.push_back( str2vec3d( positionValue ) );

        positionNode = positionNode->NextSiblingElement( "position" );
    }

    TRACE_ASSERT( positions.size() > 0 );
}

//////////////////////////////////////////////////////////////////////
// Reads a mesh list
//////////////////////////////////////////////////////////////////////
void Parser::getMeshes( const char *inputElement,
                        vector<TriangleMesh<REAL> *> &meshes,
                        vector<string> &rigidFilePrefixes,
                        vector<string> &sdfFileNames,
                        IntArray &sdfResolutions,
                        vector<string> &pulseModelFileNames,
                        FloatArray &densities,
                        map<std::string, int> &meshIDMap,
                        bool compressedFields )
{
    TiXmlElement              *root;
    TiXmlElement              *inputRoot;
    TiXmlElement              *listNode;
    TiXmlElement              *meshNode;

    root = document->FirstChildElement( "impact" );

    if ( !root )
    {
        cerr << "Error: No impact node found" << endl;
        abort();
    }

    inputRoot = root->FirstChildElement( inputElement );
    if ( !inputRoot )
    {
        cerr << "Error: No " << inputElement << " node found" << endl;
        abort();
    }

    listNode = inputRoot->FirstChildElement( "meshList" );
    if ( !listNode )
    {
        cerr << "Error: No listening position list found" << endl;
        abort();
    }

    meshNode = listNode->FirstChildElement( "mesh" );

    while ( meshNode != NULL )
    {
        string                   meshFileName;
        string                   meshID;
        string                   rigidFilePrefix;
        string                   sdfFileName;
        string                   sdfResolutionValue;
        string                   pulseModelFileName;
        string                   densityValue;
        REAL                     scale;

        meshFileName = queryRequiredAttr( meshNode, "file" );
        meshID = queryRequiredAttr( meshNode, "id" );
        rigidFilePrefix = queryRequiredAttr( meshNode, "rigidprefix" );
        sdfFileName = queryRequiredAttr( meshNode, "distancefield" );
        sdfResolutionValue = queryRequiredAttr( meshNode, "fieldresolution" );

        if ( compressedFields ) {
            pulseModelFileName = queryRequiredAttr( meshNode, "compressedpulsefile" );
        } else {
            pulseModelFileName = queryRequiredAttr( meshNode, "pulsefile" );
        }

        densityValue = queryRequiredAttr( meshNode, "density" );

        // We also get a scale value from this tag (but don't
        // require it)
        if (meshNode->QueryRealAttribute("scale", &scale) != TIXML_SUCCESS)
        {
            // Assume no additional scaling
            scale = 1.0;
        }

        printf( "Assigning index %d to mesh %s\n", (int)meshes.size(),
                meshID.c_str() );
        meshIDMap[ meshID ] = (int)meshes.size();

#if 0
        if ( pulseModelFileName != "null" && compressedFields ) {
            pulseModelFileName += "_comp";
        }
#endif

        // Construct the mesh and attempt to read it from the given file
        TriangleMesh<REAL> *mesh = new TriangleMesh<REAL>();

        if ( MeshObjReader::read( meshFileName.c_str(), *mesh,
                                  false, false, scale ) == SUCC_RETURN )
        {
            // Prepare the mesh for rendering, etc.
            mesh->generate_normals();
        } else {
            cerr << "Error reading mesh from " << meshFileName << endl;

            // Unsuccessful
            delete mesh;

            abort();
        }

        meshes.push_back( mesh );
        rigidFilePrefixes.push_back( rigidFilePrefix );
        sdfFileNames.push_back( sdfFileName );
        sdfResolutions.push_back( atoi( sdfResolutionValue.c_str() ) );
        pulseModelFileNames.push_back( pulseModelFileName );
        densities.push_back( atof( densityValue.c_str() ) );

        meshNode = meshNode->NextSiblingElement( "mesh" );
    }
}

//////////////////////////////////////////////////////////////////////
// Gets the mesh ID for each object in a scene
//////////////////////////////////////////////////////////////////////
void Parser::getObjectMeshIDs( const char *inputElement,
                               map<string, int> &meshIDMap,
                               IntArray &objectMeshIDs )
{
    TiXmlElement              *root;
    TiXmlElement              *inputRoot;
    TiXmlElement              *listNode;
    TiXmlElement              *objectNode;

    int                        objectNumber = 0;

    root = document->FirstChildElement( "impact" );
    if ( !root )
    {
        cerr << "Error: No impact node found" << endl;
        abort();
    }

    inputRoot = root->FirstChildElement( inputElement );
    if ( !inputRoot )
    {
        cerr << "Error: No " << inputElement << " node found" << endl;
        abort();
    }

    listNode = inputRoot->FirstChildElement( "objectList" );
    if ( !listNode )
    {
        cerr << "Error: No listening position list found" << endl;
        abort();
    }

    objectNode = listNode->FirstChildElement( "object" );

    while ( objectNode != NULL )
    {
        string                   meshID;

        meshID = queryRequiredAttr( objectNode, "id" );

        TRACE_ASSERT( meshIDMap.find( meshID ) != meshIDMap.end() );

        objectMeshIDs.push_back( meshIDMap[ meshID ] );

        objectNode = objectNode->NextSiblingElement( "object" );

        printf( "Assigning object %d with ID %s to mesh %d\n", objectNumber,
                meshID.c_str(), objectMeshIDs.back() );

        objectNumber += 1;
    }
    printf( "\n" );
}

#if 0
//////////////////////////////////////////////////////////////////////
// Gets a list of scatter objects from the given input element
//////////////////////////////////////////////////////////////////////
void Parser::getScatterObjects( const char *inputElement,
        vector<ScatterObject> &objects )
{
    TiXmlElement              *root;
    TiXmlElement              *elementRoot;
    TiXmlElement              *listNode;
    TiXmlElement              *objectNode;

    root = document->FirstChildElement( "impact" );
    if ( !root )
    {
        cerr << "Error: No impact node found" << endl;
        abort();
    }

    elementRoot = root->FirstChildElement( inputElement );
    if ( !elementRoot )
    {
        cerr << "Error: No " << inputElement << " node found" << endl;
        abort();
    }

    listNode = elementRoot->FirstChildElement( "scatterobjects" );
    if ( !listNode )
    {
        cerr << "Error: No scatter object list found" << endl;
        abort();
    }

    objectNode = listNode->FirstChildElement( "object" );

    while ( objectNode != NULL )
    {
        string                   directionValue;
        string                   magnitudeValue;
        string                   scaleValue;
        string                   translationValue;
        string                   absorptionValue;

        string                   sdfResolutionValue;
        string                   gridResolutionValue;
        string                   gridScaleValue;
        string                   timeStepFrequencyValue;
        string                   subStepValue;
        string                   sphereRadiusValue;
        string                   sphereImpactSpeedValue;
        string                   boundary_E_value;
        string                   boundary_nu_value;
        string                   boundaryDensityValue;

        string                   forceSDF;

        ScatterObject            object;

        object._fileName = queryRequiredAttr( objectNode, "filename" );
        scaleValue = queryOptionalAttr( objectNode, "scale", "1.0" );
        translationValue = queryOptionalAttr( objectNode,
                "translation", "0.0 0.0 0.0" );
        directionValue = queryRequiredAttr( objectNode, "acceleration_direction" );
        magnitudeValue = queryRequiredAttr( objectNode, "acceleration_magnitude" );
        absorptionValue = queryRequiredAttr( objectNode, "absorption_coefficient" );

        forceSDF = queryOptionalAttr( objectNode, "force_sdf", "0" );

        object._scale = atof( scaleValue.c_str() );
        object._translation = str2vec3d( translationValue );
        object._accelerationDirection = str2vec3d( directionValue );
        object._accelerationMagnitude = atof( magnitudeValue.c_str() );
        object._absorptionCoefficient = atof( absorptionValue.c_str() );

        if ( abs( object._accelerationMagnitude ) > 0.0
                || atoi( forceSDF.c_str() ) != 0 )
        {
            sdfResolutionValue = queryRequiredAttr( objectNode, "fieldresolution" );
            gridResolutionValue = queryRequiredAttr( objectNode, "gridresolution" );
            gridScaleValue = queryRequiredAttr( objectNode, "gridscale" );
            timeStepFrequencyValue = queryRequiredAttr( objectNode,
                    "timestepfrequency" );
            subStepValue = queryRequiredAttr( objectNode, "substeps" );
            sphereRadiusValue = queryRequiredAttr( objectNode, "sphereradius" );
            sphereImpactSpeedValue = queryRequiredAttr( objectNode, "impactspeed" );
            boundary_E_value = queryRequiredAttr( objectNode, "boundary_E" );
            boundary_nu_value = queryRequiredAttr( objectNode, "boundary_nu" );
            boundaryDensityValue = queryRequiredAttr( objectNode, "boundarydensity" );

            object._sdfResolution = atoi( sdfResolutionValue.c_str() );
            object._gridResolution = atoi( gridResolutionValue.c_str() );
            object._gridScale = atof( gridScaleValue.c_str() );
            object._timeStepFrequency = atoi( timeStepFrequencyValue.c_str() );
            object._subSteps = atoi( subStepValue.c_str() );
            object._sphereRadius = atof( sphereRadiusValue.c_str() );
            object._impactSpeed = atof( sphereImpactSpeedValue.c_str() );
            object._boundary_E = atof( boundary_E_value.c_str() );
            object._boundary_nu = atof( boundary_nu_value.c_str() );
            object._boundaryDensity = atof( boundaryDensityValue.c_str() );
        }

        objects.push_back( object );

        objectNode = objectNode->NextSiblingElement( "object" );
    }

    TRACE_ASSERT( objects.size() > 0 );
}
#endif

