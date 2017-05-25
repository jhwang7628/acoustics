#include "geometry/BoundingBox.h" 
#include "wavesolver/SimWorld.h"

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
        const Vector3d meshCentroid = meshPtr->ComputeCentroid(); 
        const int divs = (int)std::ceil(
                meshPtr->boundingSphereRadius(meshCentroid)/_simulatorSettings->cellSize
                )*2 + 6;
        const BoundingBox simUnitBox(
                _simulatorSettings->cellSize, divs, meshCentroid); 
        simUnit->simulator->InitializeSolver(simUnitBox, _simulatorSettings); 
        _simUnits.insert(std::move(simUnit)); 
    }
}

//##############################################################################
// Function Build
//##############################################################################
bool SimWorld::
StepWorld()
{
    bool continueStepping = true; 
    for (auto &u : _simUnits)
    {
        continueStepping = (u->simulator->RunForSteps(1) || continueStepping); 
    }
    _objectCollections->StepObjectStates(); 

    return continueStepping; // TODO 
}

//##############################################################################
