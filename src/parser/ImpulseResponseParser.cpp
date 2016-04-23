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
// parse meshes from xml into objects 
//##############################################################################
void ImpulseResponseParser::
GetObjects(std::shared_ptr<FDTD_Objects> objects) 
{
    if (!objects) 
        throw std::runtime_error("**ERROR** passed in pointer null for GetObjects"); 

    // get the root node
    TiXmlDocument *document2 = &_document; 
    TiXmlElement *root, *inputRoot, *listNode;
    if (!document2)
        throw std::runtime_error("**ERROR** document null"); 

    GET_FIRST_CHILD_ELEMENT_GUARD(root, document2, "impulse_response"); 
    GET_FIRST_CHILD_ELEMENT_GUARD(inputRoot, root, "rigid_object"); 
    GET_FIRST_CHILD_ELEMENT_GUARD(listNode, inputRoot, "mesh_list"); 
    const std::string meshNodeName("mesh"); 

    TiXmlElement *meshNode;
    GET_FIRST_CHILD_ELEMENT_GUARD(meshNode, listNode, meshNodeName.c_str()); 
    while (meshNode != NULL)
    {
        const std::string meshFileName = queryRequiredAttr(meshNode, "file");
        const std::string meshName = queryRequiredAttr(meshNode, "id");
        const std::string sdfFilePrefix = queryRequiredAttr(meshNode, "distancefield");
        const int sdfResolutionValue = queryRequiredInt(meshNode, "fieldresolution");
        const REAL scale = queryOptionalReal(meshNode, "scale", 1.0); 

        RigidObjectPtr object = std::make_shared<FDTD_RigidObject>(meshFileName, sdfResolutionValue, sdfFilePrefix, meshName, scale);
        object->Initialize(); 

        objects->AddObject(meshName,object); 

        meshNode = meshNode->NextSiblingElement(meshNodeName.c_str());
    }
}

//##############################################################################
// parse solver settings from xml into struct
//##############################################################################
void ImpulseResponseParser::
GetSolverSettings(std::shared_ptr<PML_WaveSolver_Settings> settings) 
{
    if (!settings) 
        throw std::runtime_error("**ERROR** pointer null for GetSolverSettings"); 
    
    // get the element nodes required
    TiXmlElement *root, *solverNode;
    TiXmlDocument *document = &_document;
    GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response"); 
    GET_FIRST_CHILD_ELEMENT_GUARD(solverNode, root, "solver"); 

    // Physical paramters 
    settings->soundSpeed = queryOptionalReal(solverNode,"sound_speed", 343.0);
    settings->airDensity = queryOptionalReal(solverNode,"density", 1.0);

    // discretization settings 
    settings->cellSize           = queryRequiredReal(solverNode, "cellsize"); 
    settings->cellDivisions      = queryRequiredInt (solverNode, "gridresolution"); 
    settings->timeEnd            = queryRequiredReal(solverNode, "stop_time"); 
    settings->timeSavePerStep    = queryRequiredInt (solverNode, "substeps"); 
    const REAL timeStepFrequency= queryRequiredReal(solverNode, "timestepfrequency"); 
    settings->timeStepSize = 1.0/timeStepFrequency; 

    // IO settings
    settings->outputPattern      = queryRequiredAttr(solverNode, "output_pattern"); 

    // Boundary settings
    settings->PML_width          = queryRequiredReal(solverNode, "PML_width"); 
    settings->PML_strength       = queryRequiredReal(solverNode, "PML_strength"); 

    // Optional settings
    settings->listeningFile              = queryOptionalAttr(solverNode, "listening_file", "NOT_SPECIFIED" );
    settings->useMesh                    = (queryOptionalInt(solverNode, "use_mesh", "1")==0) ? false : true; 
    settings->useGhostCell               = (queryOptionalInt(solverNode, "use_ghost_cell", "1")==1) ? true : false; 
    settings->cornellBoxBoundaryCondition= (queryOptionalInt(solverNode, "cornell_box_boundary_condition", "0")==1) ? true : false; 

    // set sources 
    //parms._f = queryOptionalReal( "impulse_response/solver", "f", "500" );
    //parms._sources = QueryVolumetricSource(document, this, "impulse_response/volumetric_source/source", parms._c); 
}

//##############################################################################
// parse solver settings from xml into struct
//##############################################################################
void ImpulseResponseParser::
GetSources(std::shared_ptr<PML_WaveSolver_Settings> settings) 
{
    if (!settings) 
        throw std::runtime_error("**ERROR** pointer null for GetSources"); 
   
    // TODO {
      
    //// get the element nodes required
    //TiXmlElement *root, *solverNode;
    //TiXmlDocument *document = &_document;
    //GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response"); 
    //GET_FIRST_CHILD_ELEMENT_GUARD(solverNode, root, "solver"); 

    // } TODO
}
