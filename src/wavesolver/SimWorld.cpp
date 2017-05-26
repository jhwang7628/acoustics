#include "geometry/BoundingBox.h" 
#include "wavesolver/SimWorld.h"

//##############################################################################
// Function GetSolverBBoxs
//##############################################################################
std::vector<BoundingBox> SimWorld::
GetSolverBBoxs()
{
    std::vector<BoundingBox> bboxs; 
    for (const auto &m : _simUnits)
        bboxs.push_back(m->simulator->GetSolver()->GetSolverBBox()); 
    return bboxs; 
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

    // read and initialize animator if data exists
    if (_simulatorSettings->rigidsimDataRead)
    {
        _objectCollections->InitializeAnimator(_simulatorSettings->fileDisplacement, 
                                               _simulatorSettings->fileVelocity, 
                                               _simulatorSettings->fileAcceleration);
    }

    // TEST: assign a box to each object
    for (auto &m : _objectCollections->_rigidObjects)
    {
        auto &obj = m.second; 
        std::cout << "Construct Sim unit for object: " << obj->GetMeshName()
                  << std::endl; 
        const auto &bbox = obj->GetBBox(); 
        ActiveSimUnit_UPtr simUnit = std::make_unique<ActiveSimUnit>(); 
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
        const int divs = (int)std::ceil(
                meshPtr->boundingSphereRadius(meshCentroid_o)/_simulatorSettings->cellSize
                )*2 + 6;
        const BoundingBox simUnitBox(
                _simulatorSettings->cellSize, divs, meshCentroid_w); 
        simUnit->simulator->InitializeSolver(simUnitBox, _simulatorSettings); 
        _simUnits.insert(std::move(simUnit)); 
    }
}

//##############################################################################
// Function StepWorld
//##############################################################################
bool SimWorld::
StepWorld()
{
    bool continueStepping = true; 
    for (auto &u : _simUnits)
    {
        continueStepping = (u->simulator->RunForSteps(1) || continueStepping); 
    }
    _state.time += _simulatorSettings->timeStepSize; 
    _objectCollections->StepObjectStates(_state.time); 

    return continueStepping;
}
//##############################################################################
