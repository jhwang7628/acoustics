#include "geometry/BoundingBox.h" 
#include "wavesolver/SimWorld.h"

//##############################################################################
// Static member definition
//##############################################################################
Vector3Array ListeningUnit::microphones; 
AudioOutput ListeningUnit::audioOutput; 
SimWorld::WorldRasterizer SimWorld::rasterizer; 

//##############################################################################
// Function GetSolverBBoxs
//##############################################################################
std::vector<std::pair<ActiveSimUnit_Ptr,BoundingBox>> SimWorld::
GetSolverBBoxs()
{
    std::vector<std::pair<ActiveSimUnit_Ptr,BoundingBox>> bboxs; 
    for (const auto &m : _simUnits)
        bboxs.push_back(
                std::pair<ActiveSimUnit_Ptr,BoundingBox>(
                    m,m->simulator->GetSolver()->GetSolverBBox())); 
    return bboxs; 
}

//##############################################################################
// Function UpdateSpeaker
//##############################################################################
void ActiveSimUnit::
UpdateSpeakers()
{
    if (!listen)
        return; 

    // delay line model
    const int N_mic = ListeningUnit::microphones.size(); 
    const REAL &cellSize = simulator->GetSolverSettings()->cellSize; 
    if (listen->speakers.size() != N_mic)
        listen->speakers.resize(N_mic); 
    for (int ii=0; ii<N_mic; ++ii)
    {
        const auto &mic = ListeningUnit::microphones.at(ii); 
        auto &spk = listen->speakers.at(ii); 
        spk = boxCenter + (mic - boxCenter)*(lowerRadiusBound-0.5*cellSize); 
    }
}

//##############################################################################
// Function Build
//##############################################################################
void SimWorld::
Build(ImpulseResponseParser_Ptr &parser)
{
    _objectCollections = std::make_shared<FDTD_Objects>();
    _simulatorSettings = std::make_shared<PML_WaveSolver_Settings>(); 

    // parse settings and construct objects in the scene
    parser->GetSolverSettings(_simulatorSettings); 
    parser->GetObjects(_simulatorSettings, _objectCollections); 
    SimWorld::rasterizer.cellSize = _simulatorSettings->cellSize; 
    SetWorldTime(_simulatorSettings->fastForwardToEventTime); 
    ListeningUnit::microphones = _simulatorSettings->listeningPoints; 

    // read and initialize animator if data exists
    if (_simulatorSettings->rigidsimDataRead)
    {
        _objectCollections->InitializeAnimator(_simulatorSettings->fileDisplacement, 
                                               _simulatorSettings->fileVelocity, 
                                               _simulatorSettings->fileAcceleration);
        UpdateObjectState(_state.time); 
    }

    // TEST: assign a box to each object
    for (auto &m : _objectCollections->_rigidObjects)
    {
        auto &obj = m.second; 
        std::cout << "Construct Sim unit for object: " << obj->GetMeshName()
                  << std::endl; 
        const auto &bbox = obj->GetBBox(); 
        ActiveSimUnit_Ptr simUnit = std::make_shared<ActiveSimUnit>(); 
        std::string *simulatorID = new std::string(std::to_string(m.first));
        simUnit->objects = std::make_shared<FDTD_Objects>(); 
        simUnit->objects->AddObject(std::stoi(obj->GetMeshName()), obj);  
        simUnit->simulator = std::make_shared<FDTD_AcousticSimulator>(simulatorID); 
        simUnit->simulator->SetParser(parser); 
        simUnit->simulator->SetSolverSettings(_simulatorSettings); 
        simUnit->simulator->SetSceneObjects(simUnit->objects); 
        auto meshPtr = obj->GetMeshPtr(); 
        const Vector3d meshCentroid_o = meshPtr->ComputeCentroid();
        const Vector3d meshCentroid_w = obj->ObjectToWorldPoint(meshCentroid_o);
        const Vector3d rastCentroid_w = SimWorld::rasterizer.cellCenter(
                                        SimWorld::rasterizer.rasterize(meshCentroid_w)); 
        const int divs = (int)std::ceil(
                meshPtr->boundingSphereRadius(meshCentroid_o)/_simulatorSettings->cellSize
                )*2 + 4;
        const BoundingBox simUnitBox(
                _simulatorSettings->cellSize, divs, rastCentroid_w); 
        simUnit->simulator->InitializeSolver(simUnitBox, _simulatorSettings); 
        simUnit->boxCenter = SimWorld::rasterizer.rasterize(rastCentroid_w); 
        simUnit->listen = std::make_unique<ListeningUnit>(); 
        simUnit->lowerRadiusBound = simUnitBox.minlength()/2.0; 
        simUnit->upperRadiusBound = simUnitBox.maxlength()/2.0; 
        _simUnits.insert(std::move(simUnit)); 
    }

    // setup filename for output
    char buffer[512];
    const std::string filename("all_audio.dat"); 
    snprintf(buffer, 512, _simulatorSettings->outputPattern.c_str(), 
             filename.c_str()); 
    ListeningUnit::audioOutput.OpenStream(std::string(buffer)); 
}

//##############################################################################
// Function UpdateObjectState
//##############################################################################
void SimWorld::
UpdateObjectState(const REAL &time)
{
    // logic for updating bbox and determine whether to move simbox
    _objectCollections->SetObjectStates(time); 
    for (auto &unit : _simUnits)
    {
        const auto &objects = unit->objects->GetRigidSoundObjects(); 
        Vector3d newCenter(0.,0.,0.); 
        for (const auto &objpair : objects) 
        {
            const auto &obj = objpair.second; 
            auto meshPtr = obj->GetMeshPtr(); 
            const Vector3d meshCentroid_o = meshPtr->ComputeCentroid();
            const Vector3d meshCentroid_w = obj->ObjectToWorldPoint(meshCentroid_o);
            newCenter += meshCentroid_w; 
        }
        newCenter /= (REAL)objects.size(); 
        unit->simulator->SetFieldCenter(newCenter); 
    }
}

//##############################################################################
// Function StepWorld
//##############################################################################
bool SimWorld::
StepWorld()
{
    bool continueStepping = true; 
    // update simulation
    for (auto &unit : _simUnits)
    {
        continueStepping = (unit->simulator->RunForSteps(1) || continueStepping); 
    }
    // update time and object states
    _state.time += _simulatorSettings->timeStepSize; 
    UpdateObjectState(_state.time); 

    return continueStepping;
}

//##############################################################################
// Function PreviewStepping
//##############################################################################
void SimWorld::
PreviewStepping(const int &previewSpeed) // TODO!
{
}

//##############################################################################
