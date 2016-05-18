#include <wavesolver/GaussianPressureSource.h>
#include <parser/ImpulseResponseParser.h> 

//##############################################################################
//##############################################################################
void ImpulseResponseParser::
GetPressureSources(const REAL &soundSpeed, std::vector<PressureSourcePtr> &pressureSources)
{
    // get the element nodes required
    TiXmlElement *root, *pressureSourceNode, *gaussianSourceNode;
    TiXmlDocument *document = &_document;
    GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response"); 
    GET_FIRST_CHILD_ELEMENT_GUARD(pressureSourceNode, root, "pressure_source"); 
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(gaussianSourceNode, pressureSourceNode, "gaussian_pressure_source"); 
        while (gaussianSourceNode) 
        {
            const REAL widthTime         = queryRequiredReal(gaussianSourceNode, "source_width_time"); 
            const REAL widthSpace        = queryOptionalReal(gaussianSourceNode, "source_width_space", soundSpeed*widthTime); 
            const REAL offsetTime        = queryOptionalReal(gaussianSourceNode, "source_offset_time", 0.0); 
            const bool flipSign = (queryOptionalReal(gaussianSourceNode, "source_sign_flip", 0.0) > 1E-10) ? true : false; 
            const REAL scaleSign = (flipSign) ? -1.0 : 1.0; 
            const REAL normalizeConstant = queryOptionalReal(gaussianSourceNode, "source_normalize_constant", 1.0/pow(sqrt(2.0*M_PI)*widthSpace,3))*scaleSign; 

            Vector3d sourcePosition; 
            sourcePosition.x = queryRequiredReal(gaussianSourceNode, "source_position_x"); 
            sourcePosition.y = queryRequiredReal(gaussianSourceNode, "source_position_y"); 
            sourcePosition.z = queryRequiredReal(gaussianSourceNode, "source_position_z"); 

            PressureSourcePtr sourcePtr(new GaussianPressureSource(sourcePosition, widthSpace, widthTime, offsetTime, normalizeConstant)); 
            pressureSources.push_back(std::move(sourcePtr)); 

            gaussianSourceNode = gaussianSourceNode->NextSiblingElement("gaussian_pressure_source"); 
        }
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No pressure sources found\n"; 
    }
}

//##############################################################################
// parse meshes from xml into objects 
//##############################################################################
void ImpulseResponseParser::
GetObjects(std::shared_ptr<FDTD_Objects> &objects) 
{
    if (!objects) 
    {
        std::cerr << "**WARNING** passed in pointer null for GetObjects. Initialize new object.\n"; 
        objects = std::make_shared<FDTD_Objects>();
    }

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
        const REAL initialPosition_x = queryOptionalReal(meshNode, "initial_position_x", 0.0); 
        const REAL initialPosition_y = queryOptionalReal(meshNode, "initial_position_y", 0.0); 
        const REAL initialPosition_z = queryOptionalReal(meshNode, "initial_position_z", 0.0); 

        RigidObjectPtr object = std::make_shared<FDTD_RigidObject>(meshFileName, sdfResolutionValue, sdfFilePrefix, meshName, scale);
        object->Initialize(); 
        object->ApplyTranslation(initialPosition_x, initialPosition_y, initialPosition_z); 
        objects->AddObject(meshName,object); 
        meshNode = meshNode->NextSiblingElement(meshNodeName.c_str());
    }
}

//##############################################################################
// parse solver settings from xml into struct
//##############################################################################
void ImpulseResponseParser::
GetSolverSettings(std::shared_ptr<PML_WaveSolver_Settings> &settings) 
{
    if (!settings) 
    {
        std::cerr << "**WARNING** pointer null for GetSolverSettings. Initialize new object.\n"; 
        settings = std::make_shared<PML_WaveSolver_Settings>(); 
    }
    
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
//void ImpulseResponseParser::
//GetSources(std::shared_ptr<PML_WaveSolver_Settings> settings) 
//{
//    if (!settings) 
//        throw std::runtime_error("**ERROR** pointer null for GetSources"); 
//   
//    // TODO {
//      
//    //// get the element nodes required
//    //TiXmlElement *root, *solverNode;
//    //TiXmlDocument *document = &_document;
//    //GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response"); 
//    //GET_FIRST_CHILD_ELEMENT_GUARD(solverNode, root, "solver"); 
//
//    // } TODO
//}
