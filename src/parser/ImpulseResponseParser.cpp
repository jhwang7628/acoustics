#include <parser/ImpulseResponseParser.h> 

//##############################################################################
//##############################################################################
std::vector<VolumetricSource> ImpulseResponseParser::
QueryVolumetricSource(TiXmlNode *document, Parser *parser, const std::string &path, const REAL &soundSpeed)
{
    if (!document) throw std::runtime_error("**ERROR** document not established"); 
    
    std::vector<VolumetricSource> sources; 
    
    TiXmlNode *node = document; 
    node = parser->getNodeByPath(path); 
    while (node) 
    {
        VolumetricSource source; 

        TiXmlElement *elm        = node->ToElement(); 
        source.widthTime         = parser->queryRequiredReal(elm, "source_width_time"); 
        source.widthSpace        = parser->queryOptionalReal(elm, "source_width_space", soundSpeed*source.widthTime); 
        source.offsetTime        = parser->queryOptionalReal(elm, "source_offset_time", 0.0); 
        source.normalizeConstant = parser->queryOptionalReal(elm, "source_normalize_constant", 1.0/pow(sqrt(2.0*M_PI)*source.widthSpace,3)); 
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

//##############################################################################
//##############################################################################
ImpulseResponseParser::ImpulseResponseParms ImpulseResponseParser::
Parse()
{
    ImpulseResponseParms       parms;

    // Get SDF parameters
    parms._sdfResolution = queryRequiredInt(  "impulse_response/solver", "fieldresolution" );
    parms._sdfFilePrefix = queryRequiredAttr( "impulse_response/solver", "distancefield" );

    // Get grid parameters
    parms._gridResolution   = queryRequiredInt( "impulse_response/solver", "gridresolution" );
    parms._gridScale        = queryOptionalReal( "impulse_response/solver", "gridscale", "1" );

    // Solver settings 
    parms._timeStepFrequency = queryRequiredInt( "impulse_response/solver", "timestepfrequency" );
    parms._subSteps = queryRequiredInt( "impulse_response/solver", "substeps" );
    parms._stopTime = queryRequiredReal( "impulse_response/solver", "stop_time" ); 
    parms._cellSize = queryOptionalReal( "impulse_response/solver", "cellsize", "-1"); 

    parms._outputPattern = queryRequiredAttr( "impulse_response/solver", "output_pattern" );
    parms._listeningFile = queryOptionalAttr( "impulse_response/solver", "listening_file", "none" );

    // optional parameters 
    parms._c        = queryOptionalReal( "impulse_response/solver", "c", "343.0" );
    parms._density  = queryOptionalReal( "impulse_response/solver", "density", "1" );
    parms._useMesh = (queryOptionalInt("impulse_response/solver", "use_mesh", "1")==0) ? false : true; 
    parms._cornellBoxBoundaryCondition = (queryOptionalInt("impulse_response/solver", "cornell_box_boundary_condition", "0")==1) ? true : false; 
    parms._useGhostCellBoundary = (queryOptionalInt("impulse_response/solver", "use_ghost_cell", "1")==1) ? true : false; 
    parms._f = queryOptionalReal( "impulse_response/solver", "f", "500" );

    // set sources 
    parms._sources = QueryVolumetricSource(document, this, "impulse_response/volumetric_source/source", parms._c); 

    return parms;
}

//##############################################################################
// parse meshes from xml into objects 
//##############################################################################
void ImpulseResponseParser::
GetObjects(FDTD_Objects &objects) 
{
    // get the root node
    TiXmlElement *root, *inputRoot, *listNode;
    GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response"); 
    GET_FIRST_CHILD_ELEMENT_GUARD(inputRoot, root, "rigid_object"); 
    GET_FIRST_CHILD_ELEMENT_GUARD(listNode, inputRoot, "mesh_list"); 
    const std::string meshNodeName("mesh"); 

    TiXmlElement *meshNode;
    GET_FIRST_CHILD_ELEMENT_GUARD(meshNode, listNode, meshNodeName.c_str()); 
    while ( meshNode != NULL )
    {
        const std::string meshFileName = queryRequiredAttr(meshNode, "file");
        const std::string meshName = queryRequiredAttr(meshNode, "id");
        const std::string sdfFilePrefix = queryRequiredAttr(meshNode, "distancefield");
        const int sdfResolutionValue = queryRequiredInt(meshNode, "fieldresolution");
        const REAL scale = queryOptionalReal(meshNode, "scale", 1.0); 

        RigidObjectPtr object(new FDTD_RigidObject(meshFileName, sdfResolutionValue, sdfFilePrefix, meshName, scale));
        object->Initialize(); 

        objects.AddObject(meshName,object); 

        meshNode = meshNode->NextSiblingElement(meshNodeName.c_str());
    }
}








