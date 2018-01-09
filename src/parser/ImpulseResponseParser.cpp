#include <wavesolver/GaussianPressureSource.h>
#include <wavesolver/SpeakerPressureSource.h>
#include <parser/ImpulseResponseParser.h>
#include <io/ImpulseSeriesReader.h>
#include <wavesolver/VibrationalSource.h>
#include <wavesolver/ModalVibrationalSource.h>
#include <wavesolver/AccelerationNoiseVibrationalSource.h>
#include <wavesolver/ShellVibrationalSource.h>
#include <wavesolver/WaterVibrationalSource.h>
#include <wavesolver/FDTD_PlaneConstraint.h>
#include <wavesolver/WaterVibrationalSourceBubbles.h>
#include <wavesolver/SpeakerVibrationalSource.h>
#include <modal_model/SparseModalEncoder.h>

//##############################################################################
// parse meshes from xml into objects
//##############################################################################
void ImpulseResponseParser::
GetObjects(const std::shared_ptr<PML_WaveSolver_Settings> &solverSettings, std::shared_ptr<FDTD_Objects> &objects)
{
    if (!objects)
    {
        std::cerr << "**MESSAGE** passed in pointer null for GetObjects. Initialize new object.\n";
        objects = std::make_shared<FDTD_Objects>();
    }

    // get the root node
    TiXmlDocument *document2 = &_document;
    TiXmlElement *root, *inputRoot;
    if (!document2)
        throw std::runtime_error("**ERROR** document null");

    GET_FIRST_CHILD_ELEMENT_GUARD(root, document2, "impulse_response");
    GET_FIRST_CHILD_ELEMENT_GUARD(inputRoot, root, "scene");

    const std::string rigidSoundObjectNodeName("rigid_sound_object");
    const std::string rigidObjectNodeName("rigid_object");

    // parse and build rigid sound objects
    TiXmlElement *rigidSoundObjectNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(rigidSoundObjectNode, inputRoot, rigidSoundObjectNodeName.c_str());
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No rigid_sound_object found\n";
    }

    int meshCount = 0;
    std::vector<ImpulseSeriesReader> readers;
    while (rigidSoundObjectNode != NULL)
    {
        //const std::string meshFileName = queryRequiredAttr(rigidSoundObjectNode, "file");
        //const std::string sdfFilePrefix = queryRequiredAttr(rigidSoundObjectNode, "distancefield");
        //const int meshID = queryRequiredInt(rigidSoundObjectNode, "id");
        const std::string meshName = std::to_string(meshCount++);
        const std::string workingDirectory = queryRequiredAttr(rigidSoundObjectNode, "working_directory");
        const std::string objectPrefix = queryRequiredAttr(rigidSoundObjectNode, "object_prefix");
        const int sdfResolutionValue = queryRequiredInt(rigidSoundObjectNode, "fieldresolution");
        const REAL scale = queryOptionalReal(rigidSoundObjectNode, "scale", 1.0);
        const REAL initialPosition_x = queryOptionalReal(rigidSoundObjectNode, "initial_position_x", 0.0);
        const REAL initialPosition_y = queryOptionalReal(rigidSoundObjectNode, "initial_position_y", 0.0);
        const REAL initialPosition_z = queryOptionalReal(rigidSoundObjectNode, "initial_position_z", 0.0);
        FDTD_RigidObject::OptionalAttributes attr;
        attr.isFixed = (queryOptionalReal(rigidSoundObjectNode, "fixed", 0.0) > 1E-10) ? true : false;

        const bool buildFromTetMesh = true;
        RigidObjectPtr object_r = std::make_shared<FDTD_RigidSoundObject>(workingDirectory, sdfResolutionValue, objectPrefix, buildFromTetMesh, solverSettings, meshName, scale);
        RigidSoundObjectPtr object = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(object_r);
        object->SetOptionalAttributes(attr);
        object->SetAnimated(true);
        // load impulse from file
        const std::string impulseFile = queryRequiredAttr(rigidSoundObjectNode, "impulse_file");
        const std::string rigidsimConfigFile = queryRequiredAttr(rigidSoundObjectNode, "impulse_rigidsim_config_file");
        ImpulseSeriesReader reader(impulseFile, rigidsimConfigFile);
        std::shared_ptr<ImpulseSeriesObject> objectPtr = std::static_pointer_cast<ImpulseSeriesObject>(object);
        readers.push_back(reader);
        //reader.LoadImpulses(meshID, objectPtr, objects);
        //REAL impulseRangeStart, impulseRangeStop;
        //object->GetRangeOfImpulses(impulseRangeStart, impulseRangeStop);
        //std::cout << "Impulses Read for object " << meshName << ":\n"
        //          << " Number of impulses: " << object->Size() << "\n"
        //          << " Time step size for rigid sim: " << object->GetRigidsimTimeStepSize() << "\n"
        //          << " Time range of impulses: [" << impulseRangeStart << ", " << impulseRangeStop << "]\n"
        //          << "\n";

        // load modes from file
        const std::string modeFile = queryRequiredAttr(rigidSoundObjectNode, "mode_file");
        ModalMaterialList materials;
        GetModalMaterials(materials);
        const int materialID = queryRequiredInt(rigidSoundObjectNode, "material_id");
        auto materialPtr = materials.at(materialID);

        const REAL ODEStepSize = 1.0/queryRequiredReal(rigidSoundObjectNode, "modal_ODE_step_frequency");
        object->ModalAnalysisObject::Initialize(ODEStepSize, modeFile, materialPtr);
        object->FDTD_RigidSoundObject::Initialize();
        object->ApplyTranslation(initialPosition_x, initialPosition_y, initialPosition_z);

        // get source and attach to objects
        const bool has_modal_source     = queryRequiredBool(rigidSoundObjectNode, "has_modal_source"    );
        const bool has_acc_noise_source = queryRequiredBool(rigidSoundObjectNode, "has_acc_noise_source");
        if (has_modal_source)
        {
            std::cout << "has modal source\n";
            VibrationalSourcePtr sourcePtr(new ModalVibrationalSource(object));
            object->AddVibrationalSource(sourcePtr);
        }
        if (has_acc_noise_source)
        {
            std::cout << "has acc noise source\n";
            VibrationalSourcePtr sourcePtr(new AccelerationNoiseVibrationalSource(object));
            object->AddVibrationalSource(sourcePtr);
        }

        objects->AddObject(std::stoi(meshName), object_r);
        rigidSoundObjectNode = rigidSoundObjectNode->NextSiblingElement(rigidSoundObjectNodeName.c_str());
    }

    const int N_rigidSoundObject = objects->N();
    for (int o_idx=0; o_idx<N_rigidSoundObject; ++o_idx)
    {
        auto object = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(objects->GetPtr(o_idx));
        readers.at(o_idx).LoadImpulses(o_idx, object, objects);
        REAL impulseRangeStart, impulseRangeStop;
        object->GetRangeOfImpulses(impulseRangeStart, impulseRangeStop);
#ifdef DEBUG_PRINT
        std::cout << static_cast<ModalAnalysisObject>(*(object.get())) << std::endl;
        std::cout << "Impulses Read for object " << objects->GetMeshName(o_idx) << ":\n"
                  << " Number of impulses: " << object->N_Impulses() << "\n"
                  << " Time step size for rigid sim: " << object->GetRigidsimTimeStepSize() << "\n"
                  << " Time range of impulses: [" << impulseRangeStart << ", " << impulseRangeStop << "]\n"
                  << "\n";
        object->PrintAllImpulses();
#endif
    }

    // parse and build rigid objects. In the implementation, I still use the class
    // RigidSoundObject to store the non-sounding objects.
    TiXmlElement *rigidObjectNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(rigidObjectNode, inputRoot, rigidObjectNodeName.c_str());
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No rigid_object found\n";
    }
    while (rigidObjectNode != NULL)
    {
        const std::string meshName = std::to_string(meshCount++);
        const std::string workingDirectory = queryRequiredAttr(rigidObjectNode, "working_directory");
        const std::string objectPrefix = queryRequiredAttr(rigidObjectNode, "object_prefix");
        const int sdfResolutionValue = queryRequiredInt(rigidObjectNode, "fieldresolution");
        const REAL scale = queryOptionalReal(rigidObjectNode, "scale", 1.0);
        const REAL initialPosition_x = queryOptionalReal(rigidObjectNode, "initial_position_x", 0.0);
        const REAL initialPosition_y = queryOptionalReal(rigidObjectNode, "initial_position_y", 0.0);
        const REAL initialPosition_z = queryOptionalReal(rigidObjectNode, "initial_position_z", 0.0);
        FDTD_RigidObject::OptionalAttributes attr;
        attr.isFixed = (queryOptionalReal(rigidObjectNode, "fixed", 0.0) > 1E-10) ? true : false;
        attr.isThinStructure = queryOptionalBool(rigidObjectNode, "thin_structure", "0");

        const bool buildFromTetMesh = false;
        RigidObjectPtr object = std::make_shared<FDTD_RigidSoundObject>(workingDirectory, sdfResolutionValue, objectPrefix, buildFromTetMesh, solverSettings, meshName, scale, false);
        object->SetOptionalAttributes(attr);
        object->ApplyTranslation(initialPosition_x, initialPosition_y, initialPosition_z);

        // get speaker shader
        const bool hasSpeakerShader = queryOptionalBool(rigidObjectNode, "has_speaker_shader", "0");
        if (hasSpeakerShader)
        {
            const std::string speakerFile = queryRequiredAttr(rigidObjectNode, "speaker_file");
            IntArray vids;
            try
            {
                vids = queryRequiredIntListFromFile(rigidObjectNode, "speaker_vids_file");
            }
            catch (std::runtime_error &e)
            {
                vids = queryRequiredIntList(rigidObjectNode, "speaker_vids");
            }
            assert(vids.size()>0);

            VibrationalSourcePtr source(new SpeakerVibrationalSource(object));
            std::dynamic_pointer_cast<SpeakerVibrationalSource>(source)->Initialize(speakerFile, vids);
            object->AddVibrationalSource(source);
        }

        // get kinematics if any
        const std::string fileKinematics = queryOptionalAttr(rigidObjectNode, "file_kinematics", "N/A");
        if (fileKinematics != "N/A")
        {
            const REAL stepSize = 1.0/queryRequiredReal(rigidObjectNode, "kinematics_frame_rate");
            solverSettings->rigidsimDataRead = false; // this conflicts with rigid sim data, so disable it
            solverSettings->kinFileExists = true;

            KinematicsMetadata meta;
            meta.objId = meshName;
            meta.fileKinematics = fileKinematics;
            meta.stepSize = stepSize;
            solverSettings->objKinematicsMetadata[meshName] = meta;
            object->SetAnimated(true);
        }

        objects->AddObject(std::stoi(meshName), object);
        rigidObjectNode = rigidObjectNode->NextSiblingElement(rigidObjectNodeName.c_str());
    }

    // parse and build rigid object sequence. This class of object uses a sequence of (rigid) objs to
    // represent the deformation of an object.
    const std::string rigidObjectSeqNodeName("rigid_object_sequence");
    TiXmlElement *rigidObjectSeqNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(rigidObjectSeqNode, inputRoot, rigidObjectSeqNodeName.c_str());
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No rigid_object_seq found\n";
    }
    while (rigidObjectSeqNode != NULL)
    {
        const std::string meshName = std::to_string(meshCount++);
        const std::string workingDirectory = queryRequiredAttr(rigidObjectSeqNode, "working_directory");
        const std::string sequencePrefix = queryRequiredAttr(rigidObjectSeqNode, "obj_sequence_prefix");
        const REAL sequenceFramerate = 1./queryRequiredReal(rigidObjectSeqNode,  "obj_sequence_frequency");
        const int sdfResolutionValue = queryOptionalInt(rigidObjectSeqNode, "fieldresolution", "400");
        const REAL scale = queryOptionalReal(rigidObjectSeqNode, "scale", 1.0);
        const REAL initialPosition_x = queryOptionalReal(rigidObjectSeqNode, "initial_position_x", 0.0);
        const REAL initialPosition_y = queryOptionalReal(rigidObjectSeqNode, "initial_position_y", 0.0);
        const REAL initialPosition_z = queryOptionalReal(rigidObjectSeqNode, "initial_position_z", 0.0);
        const bool buildFromTetMesh = false;
        FDTD_RigidObject::OptionalAttributes attr;
        attr.isFixed = false;
        attr.isThinStructure = queryOptionalBool(rigidObjectSeqNode, "thin_structure", "0");

        // get speaker shader, use this shader to load rigid objects (similar to water bubbles shader)
        // NOTE: this is not a good abstraction because object and shader should be separate, however,
        //       I am doing it this way to keep it consistent with the water bubble shader.
        // for this reason, speaker shader needs to exist for the rigid object sequency type
        const bool hasSpeakerShader = queryOptionalBool(rigidObjectSeqNode, "has_speaker_shader", "0");
        if (!hasSpeakerShader)
            throw std::runtime_error("**ERROR** Does not support rigid object sequence without speaker shader source");

        const std::string speakerFile = queryRequiredAttr(rigidObjectSeqNode, "speaker_file");
        const std::string speakerVIdsDir = queryRequiredAttr(rigidObjectSeqNode, "speaker_vids_file_directory");
        const std::string speakerVIdsSuf = queryRequiredAttr(rigidObjectSeqNode, "speaker_vids_file_suffix");

        VibrationalSourcePtr source(new SpeakerVibrationalSource());
        auto speakerSrc = std::dynamic_pointer_cast<SpeakerVibrationalSource>(source);
        speakerSrc->SetSeqSampleRate(sequenceFramerate);
        SpeakerVibrationalSource::DataStep firstStep =
            speakerSrc->ReadObjSeqMetaData(workingDirectory, sequencePrefix, speakerVIdsDir, speakerVIdsSuf);
        speakerSrc->Initialize(speakerFile, firstStep.handles);

        // initialize object and connect object/source
        RigidObjectPtr object =
            std::make_shared<FDTD_RigidSoundObject>(
                    workingDirectory,
                    sdfResolutionValue,
                    firstStep.objFilePrefix,
                    buildFromTetMesh,
                    solverSettings,
                    meshName,
                    scale,
                    false
                    );
        object->SetOptionalAttributes(attr);
        object->ApplyTranslation(initialPosition_x, initialPosition_y, initialPosition_z);
        source->SetOwner(object);
        object->AddVibrationalSource(source);
        objects->AddObject(std::stoi(meshName), object);
        rigidObjectSeqNode = rigidObjectSeqNode->NextSiblingElement(rigidObjectSeqNodeName.c_str());
    }

    // parse and build water surface objects. Right now use rigid object to
    // test, but water surface might deform so its not suitable.  // TODO
    const std::string waterSurfaceObjectNodeName("water_surface_object");
    TiXmlElement *waterSurfaceObjectNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(waterSurfaceObjectNode, inputRoot, waterSurfaceObjectNodeName.c_str());
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No water_surface_object found\n";
    }
    while (waterSurfaceObjectNode != NULL)
    {
        const std::string meshName = std::to_string(meshCount++);
        const std::string workingDirectory = queryRequiredAttr(waterSurfaceObjectNode, "working_directory");
        const std::string objectPrefix = queryRequiredAttr(waterSurfaceObjectNode, "object_prefix");
        const int sdfResolutionValue = queryRequiredInt(waterSurfaceObjectNode, "fieldresolution");
        const REAL scale = queryOptionalReal(waterSurfaceObjectNode, "scale", 1.0);
        const REAL initialPosition_x = queryOptionalReal(waterSurfaceObjectNode, "initial_position_x", 0.0);
        const REAL initialPosition_y = queryOptionalReal(waterSurfaceObjectNode, "initial_position_y", 0.0);
        const REAL initialPosition_z = queryOptionalReal(waterSurfaceObjectNode, "initial_position_z", 0.0);
        const std::string inputRecordingFile = queryRequiredAttr(waterSurfaceObjectNode, "input_recording");
        const REAL decayRadius = queryOptionalReal(waterSurfaceObjectNode, "decay_radius", 5.0);

        const bool buildFromTetMesh = false;
        RigidObjectPtr object = std::make_shared<FDTD_RigidSoundObject>(workingDirectory, sdfResolutionValue, objectPrefix, buildFromTetMesh, solverSettings, meshName, scale);
        object->ApplyTranslation(initialPosition_x, initialPosition_y, initialPosition_z);
        VibrationalSourcePtr sourcePtr = std::make_shared<WaterVibrationalSource>(object, inputRecordingFile, decayRadius);
        object->AddVibrationalSource(sourcePtr);
        objects->AddObject(std::stoi(meshName), object);
        waterSurfaceObjectNode = waterSurfaceObjectNode->NextSiblingElement(waterSurfaceObjectNodeName.c_str());
    }

    // build shells
    std::string name;
    TiXmlElement *node;
    name = "shell_object";
    try {
        GET_FIRST_CHILD_ELEMENT_GUARD(node, inputRoot, name.c_str());
    } catch (const std::runtime_error &error) {
        std::cout << "No shell objects found\n";
    }
    while (node != NULL)
    {
        const std::string meshName           = std::to_string(meshCount++);
        const std::string workingDirectory   = queryRequiredAttr(node, "working_directory"        );
        const std::string objectPrefix       = queryRequiredAttr(node, "object_prefix"            );
        const std::string shellDataDirectory = queryRequiredAttr(node, "shell_data_directory"     );
        const bool has_shell_source          = queryOptionalBool(node, "has_shell_source"    , "1");
        RigidObjectPtr object = std::make_shared<FDTD_ShellObject>(workingDirectory,
                                                                   shellDataDirectory,
                                                                   -1,
                                                                   objectPrefix,
                                                                   solverSettings,
                                                                   meshName);
        auto shell_object = std::dynamic_pointer_cast<FDTD_ShellObject>(object);
        shell_object->Initialize();
        if (has_shell_source)
        {
            VibrationalSourcePtr sourcePtr(new ShellVibrationalSource(shell_object));
            object->AddVibrationalSource(sourcePtr);
        }
        objects->AddObject(std::stoi(meshName), object);
        node = node->NextSiblingElement(name.c_str());
    }

    // build constraints
    name = "constraint";
    try {
        GET_FIRST_CHILD_ELEMENT_GUARD(node, inputRoot, name.c_str());
    } catch (const std::runtime_error &error) {
        std::cout << "No constraint found\n";
    }
    while (node != NULL)
    {
        const std::string id = queryOptionalAttr(node, "id", std::to_string(objects->N_constraints()));
        const int direction = queryOptionalInt(node, "direction", "1");
        const int sign = queryOptionalInt(node, "sign", "1");
        const REAL height = queryOptionalReal(node, "height", 0.0);
        auto constraint = std::make_shared<FDTD_PlaneConstraint>(direction, sign, height);
        objects->AddConstraint(id, constraint);
        node = node->NextSiblingElement(name.c_str());
    }

    // Load data from the previous bubbles project
    // Currently using a static mesh
    const std::string bubblesNodeName("water_surface_bubbles_object");
    TiXmlElement *bubblesNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(bubblesNode, inputRoot, bubblesNodeName.c_str());
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No water_surface_bubbles_object found\n";
    }
    while (bubblesNode != NULL)
    {
        const std::string meshName = std::to_string(meshCount++);
        const std::string workingDirectory = queryRequiredAttr(bubblesNode, "working_directory");
        const std::string objectPrefix = queryRequiredAttr(bubblesNode, "object_prefix");
        const int sdfResolutionValue = queryRequiredInt(bubblesNode, "fieldresolution");
        const REAL scale = queryOptionalReal(bubblesNode, "scale", 1.0);
        const REAL initialPosition_x = queryOptionalReal(bubblesNode, "initial_position_x", 0.0);
        const REAL initialPosition_y = queryOptionalReal(bubblesNode, "initial_position_y", 0.0);
        const REAL initialPosition_z = queryOptionalReal(bubblesNode, "initial_position_z", 0.0);
        const std::string dataDir = queryRequiredAttr(bubblesNode, "data_dir");
        FDTD_RigidObject::OptionalAttributes attr;
        //attr.isFixed = (queryOptionalReal(rigidObjectNode, "fixed", 0.0) > 1E-10) ? true : false;
        attr.isFixed = false;

        const bool buildFromTetMesh = false;
        RigidObjectPtr object = std::make_shared<FDTD_RigidSoundObject>(workingDirectory, sdfResolutionValue, objectPrefix, buildFromTetMesh, solverSettings, meshName, scale, false);
        object->SetOptionalAttributes(attr);
        object->ApplyTranslation(initialPosition_x, initialPosition_y, initialPosition_z);
        VibrationalSourcePtr sourcePtr = std::make_shared<WaterVibrationalSourceBubbles>(object, dataDir, workingDirectory);
        object->AddVibrationalSource(sourcePtr);
        objects->AddObject(std::stoi(meshName), object);
        bubblesNode = bubblesNode->NextSiblingElement(bubblesNodeName.c_str());
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
        std::cerr << "**MESSAGE** pointer null for GetSolverSettings. Initialize new object.\n";
        settings = std::make_shared<PML_WaveSolver_Settings>();
    }

    // get the element nodes required
    TiXmlElement *root, *solverNode, *encoderNode;
    TiXmlDocument *document = &_document;
    GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response");
    GET_FIRST_CHILD_ELEMENT_GUARD(solverNode, root, "solver");

    // Physical paramters
    settings->soundSpeed = queryOptionalReal(solverNode,"sound_speed", 343.0);
    settings->airDensity = queryOptionalReal(solverNode,"density", 1.0);

    // discretization settings
    settings->cellSize           = queryRequiredReal(solverNode, "cellsize");
    settings->timeEnd            = queryRequiredReal(solverNode, "stop_time");
    settings->timeSavePerStep    = queryOptionalInt (solverNode, "substeps", "-1");
    settings->numberTimeSteps    = queryOptionalInt (solverNode, "number_of_timesteps", "0");
    const REAL timeStepFrequency= queryRequiredReal(solverNode, "timestepfrequency");
    settings->timeStepSize = 1.0/timeStepFrequency;
    settings->FV_boundarySubdivision = queryOptionalInt(solverNode, "fv_boundary_subdivision", "5");

    // IO settings
    settings->outputPattern = queryRequiredAttr(solverNode, "output_pattern");
    settings->writePressureFieldToDisk = (queryRequiredInt(solverNode, "write_pressure_field_to_disk")==0) ? false : true;

    // time parallelization settings
    settings->timeParallel        = queryOptionalInt(solverNode, "time_parallel", "0") != 0;
    if(settings->timeParallel)
    {
        settings->numTimeChunks   = queryRequiredInt(solverNode, "num_time_chunks") ;
        settings->modalParallelErrorTol = queryOptionalReal(solverNode, "modal_parallel_error_tol", 1e-7);
        settings->modalMinDampingPeriods = queryOptionalInt(solverNode, "modal_min_damping_periods", "2");
        settings->overlapTime = queryOptionalReal(solverNode, "overlap_time", 6e-3);
        settings->earlyTerminationRatio = queryOptionalReal(solverNode, "early_termination_ratio", 1E-3);
    }

    // Boundary settings
    settings->PML_width          = queryRequiredReal(solverNode, "PML_width");
    settings->PML_strength       = queryRequiredReal(solverNode, "PML_strength");

    // Optional settings
    settings->alpha                    = queryOptionalReal(solverNode, "air_viscosity_alpha"       , 0.0        );
    settings->useMesh                  = queryOptionalBool(solverNode, "use_mesh"                  , "1"        );
    settings->useGhostCell             = queryOptionalBool(solverNode, "use_ghost_cell"            , "1"        );
    settings->validateUsingFBem        = queryOptionalBool(solverNode, "validate_using_fbem"       , "0"        );
    settings->boundaryConditionPreset  = queryOptionalInt(solverNode, "boundary_condition_preset"  , "0"        );
    settings->fastForwardToEventTime   = queryOptionalReal(solverNode, "fast_forward_to_event_time", 0.0        );
    settings->adaptiveStartTime        = queryOptionalBool(solverNode, "adaptive_start_time"       , "0"        );
    const std::string boundaryHandling = queryOptionalAttr(solverNode, "boundary_handling"         , "rasterize");
    if (settings->alpha > 0.0)
        settings->useAirViscosity = true;
    if (boundaryHandling.compare("rasterize")==0)
        settings->boundaryHandlingType = PML_WaveSolver_Settings::BoundaryHandling::RASTERIZE;
    else if (boundaryHandling.compare("fully_coupled")==0)
        settings->boundaryHandlingType = PML_WaveSolver_Settings::BoundaryHandling::FULLY_COUPLED;
    else
        throw std::runtime_error("**ERROR** boundary handling type not understood: " + boundaryHandling);
    settings->onlyObjSequence = queryOptionalBool(solverNode, "only_obj_sequence", "0");

    // set sources
    //parms._f = queryOptionalReal( "impulse_response/solver", "f", "500" );
    //parms._sources = QueryVolumetricSource(document, this, "impulse_response/volumetric_source/source", parms._c);

    // parse object rigidsim results data if exists
    TiXmlElement *sceneNode, *rigidsimDataNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(sceneNode, root, "scene");
        GET_FIRST_CHILD_ELEMENT_GUARD(rigidsimDataNode, sceneNode, "object_rigidsim_data");
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No scene/object_rigidsim_data node found\n";
    }
    if (rigidsimDataNode)
    {
        settings->rigidsimDataRead = true;
        settings->fileDisplacement = queryRequiredAttr(rigidsimDataNode, "file_displacement");
        settings->fileVelocity = queryRequiredAttr(rigidsimDataNode, "file_velocity");
        settings->fileAcceleration = queryRequiredAttr(rigidsimDataNode, "file_acceleration");
    }

    // parse listening shell config
    TiXmlElement *listNode;
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(listNode, root, "listening_shell");
        if (listNode)
        {
            settings->useShell = true;
            settings->refShellFile = queryRequiredAttr(listNode, "reference_shell_file");
            settings->spacing      = queryRequiredReal(listNode, "finite_difference_spacing");
        }
    }
    catch (const std::runtime_error &e) { }

    // parse and construct modal encoder
    bool hasEncoderDefined;
    GET_FIRST_CHILD_ELEMENT_FLAG(encoderNode, root, "modal_encoder", hasEncoderDefined);
    if (hasEncoderDefined)
    {
        SparseModalEncoder::useEncoder = true;
        SparseModalEncoder::rank       = queryRequiredInt(encoderNode, "rank");
        SparseModalEncoder::epsilon    = queryRequiredReal(encoderNode, "epsilon");
    }

    // parse solver control policy
    TiXmlElement *node;
    GET_FIRST_CHILD_ELEMENT_GUARD(node, root, "solver_control_policy");
    std::string ptype = queryRequiredAttr(node, "type");
    if (ptype == "dynamic")
    {
        settings->solverControlPolicy = std::make_shared<Dynamic_Policy>();
        std::shared_ptr<Dynamic_Policy> p = std::dynamic_pointer_cast<Dynamic_Policy>(settings->solverControlPolicy);
        p->padding = queryOptionalInt(node, "padding", "20");
        p->type = "dynamic";
    }
    else if (ptype == "static")
    {
        settings->solverControlPolicy = std::make_shared<Static_Policy>();
        auto p = std::dynamic_pointer_cast<Static_Policy>(settings->solverControlPolicy);
        const REAL domainCenter_x    = queryOptionalReal(node, "domain_center_x", 0.0);
        const REAL domainCenter_y    = queryOptionalReal(node, "domain_center_y", 0.0);
        const REAL domainCenter_z    = queryOptionalReal(node, "domain_center_z", 0.0);
        p->cellDivisions = queryRequiredInt (node, "gridresolution");
        p->domainCenter.set(domainCenter_x, domainCenter_y, domainCenter_z);
        p->type = "static";
    }
    else if (ptype == "markers")
    {
        settings->solverControlPolicy = std::make_shared<Markers_Policy>();
        auto p = std::dynamic_pointer_cast<Markers_Policy>(settings->solverControlPolicy);
        p->boxSizeMarker = queryRequiredReal(node, "box_size_each_marker");
        p->markersDir    = queryRequiredAttr(node, "markers_dir"         );
        p->frameRate     = queryRequiredReal(node, "markers_frame_rate"  );
        p->startTime     = queryRequiredReal(node, "markers_start_time"  );
        p->type = "markers";
    }
    else
    {
        throw std::runtime_error("**ERRORR** solver control policy not understood");
    }
}

//##############################################################################
//##############################################################################
void ImpulseResponseParser::
GetPressureSources(const REAL &soundSpeed, std::vector<PressureSourcePtr> &pressureSources)
{
    TiXmlElement *root, *pressureSourceNode, *gaussianSourceNode;
    TiXmlDocument *document = &_document;
    GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response");
    pressureSourceNode = root->FirstChildElement("pressure_source");
    // parse gaussian src
    {
        gaussianSourceNode = pressureSourceNode->FirstChildElement("gaussian_pressure_source");
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

    // parse speaker src
    {
        TiXmlElement *node;
        node = pressureSourceNode->FirstChildElement("speaker_pressure_source");
        while (node)
        {
            const REAL t0 = queryRequiredReal(node, "start_time");
            const REAL w  = queryRequiredReal(node, "width_space");
            const std::string f = queryRequiredAttr(node, "speaker_file");
            PressureSourcePtr src(new SpeakerPressureSource());
            auto spksrc = std::dynamic_pointer_cast<SpeakerPressureSource>(src);
            spksrc->SetStartTime(t0);
            spksrc->SetWidthSpace(w);
            spksrc->SetSoundSpeed(soundSpeed);
            const std::string type = queryOptionalAttr(node, "positioning_type", "static");
            std::shared_ptr<SpeakerPressureSource::PositioningControl> control;
            if (type == "markers")
            {
                control = std::make_shared<SpeakerPressureSource::MarkersPositioningControl>();
                auto mControl =
                    std::dynamic_pointer_cast<SpeakerPressureSource::MarkersPositioningControl>(control);
                mControl->type      = type;
                mControl->dir       = queryRequiredAttr(node, "markers_dir");
                mControl->frameRate = queryRequiredReal(node, "markers_frame_rate");
                mControl->startTime = queryRequiredReal(node, "markers_start_time");
            }
            else if (type == "static")
            {
                control = std::make_shared<SpeakerPressureSource::PositioningControl>();
                Vector3d &pos = control->position;
                pos.x = queryRequiredReal(node, "position_x");
                pos.y = queryRequiredReal(node, "position_y");
                pos.z = queryRequiredReal(node, "position_z");
            }
            else
            {
                throw std::runtime_error("**ERROR** speaker position type not understood");
            }
            spksrc->SetPosControl(control);
            spksrc->Initialize(f);

            pressureSources.push_back(std::move(src));
            node = node->NextSiblingElement("speaker_pressure_source");
        }
    }
}

//##############################################################################
//##############################################################################
void ImpulseResponseParser::
GetListeningPoints(Vector3Array &listeningPoints)
{
    listeningPoints.clear();
    // get the element nodes required
    TiXmlElement *root, *listNode, *node;
    TiXmlDocument *document = &_document;
    GET_FIRST_CHILD_ELEMENT_GUARD(root, document, "impulse_response");
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(listNode, root, "listening_point_list");
        GET_FIRST_CHILD_ELEMENT_GUARD(node, listNode, "listening_point");
        while (node)
        {
            Vector3<REAL> position;
            position.x = queryRequiredReal(node, "x");
            position.y = queryRequiredReal(node, "y");
            position.z = queryRequiredReal(node, "z");
            listeningPoints.push_back(position);
            node = node->NextSiblingElement("listening_point");
        }
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No listening points found\n";
    }
}

//##############################################################################
//##############################################################################
void ImpulseResponseParser::
GetModalMaterials(ModalMaterialList &modalMaterials)
{
    // get the root node
    TiXmlDocument *document2 = &_document;
    TiXmlElement *root, *listNode;
    if (!document2)
        throw std::runtime_error("**ERROR** document null");

    GET_FIRST_CHILD_ELEMENT_GUARD(root, document2, "impulse_response");
    GET_FIRST_CHILD_ELEMENT_GUARD(listNode, root, "material_list");
    const std::string materialNodeName("modal_material");
    try
    {
        TiXmlElement *materialNode;
        GET_FIRST_CHILD_ELEMENT_GUARD(materialNode, listNode, materialNodeName.c_str());
        while (materialNode != NULL)
        {
            std::shared_ptr<ModalMaterial> material = std::make_shared<ModalMaterial>();
            material->alpha = queryRequiredReal(materialNode, "alpha");
            material->beta = queryRequiredReal(materialNode, "beta");
            material->density = queryRequiredReal(materialNode, "density");
            material->poissonRatio = queryRequiredReal(materialNode, "poisson_ratio");
            material->youngsModulus = queryRequiredReal(materialNode, "youngs_modulus");
            material->inverseDensity = 1./material->density;
            material->one_minus_nu2_over_E = (1.0 - pow(material->poissonRatio, 2)) / material->youngsModulus;
            material->id = modalMaterials.size();
            modalMaterials.push_back(material);
            materialNode = materialNode->NextSiblingElement(materialNodeName.c_str());
        }
    }
    catch (const std::runtime_error &error)
    {
        std::cout << "No materials found\n";
    }
}

//##############################################################################
//##############################################################################
ChunkPartitionParam_Ptr ImpulseResponseParser::
GetChunkPartitionParam()
{
    // get the root node
    TiXmlDocument *document2 = &_document;
    TiXmlElement *root, *node;
    if (!document2)
        throw std::runtime_error("**ERROR** document null");

    auto param = std::make_shared<ChunkPartitionParam>();
    try
    {
        GET_FIRST_CHILD_ELEMENT_GUARD(root, document2, "impulse_response");
        GET_FIRST_CHILD_ELEMENT_GUARD(node, root, "time_parallelization");
        GET_FIRST_CHILD_ELEMENT_GUARD(node, node, "chunk_partition");

        param->adaptive = queryRequiredBool(node, "adaptive");
        param->L_z      = queryRequiredReal(node, "L_z"     );
        param->N_0      = queryRequiredInt (node, "N_0"     );
        param->N_maxc   = queryRequiredInt (node, "N_maxc"  );
        param->outFile  = queryRequiredAttr(node, "out_file");
    }
    catch (std::runtime_error &e)
    {
        std::cout << "adaptive time-parallel = false\n";
    }
    return param;
}
