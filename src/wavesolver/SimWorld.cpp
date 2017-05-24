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
        LOGGING_INFO(
                "Construct Sim unit for object: %s\n", obj->GetMeshName().c_str()
        ); 
        const auto &bbox = obj->GetBBox(); 
        ActiveSimUnit_UPtr simUnit = std::make_unique<ActiveSimUnit>(); 
        simUnit->objects = std::make_shared<FDTD_Objects>(); 
        simUnit->objects->AddObject(std::stoi(obj->GetMeshName()), obj);  
        simUnit->simulator = std::make_shared<FDTD_AcousticSimulator>(); 
        simUnit->simulator->SetParser(parser); 
        simUnit->simulator->SetSolverSettings(_simulatorSettings); 
        simUnit->simulator->SetSceneObjects(simUnit->objects); 
        simUnit->simulator->InitializeSolver(_simulatorSettings); 
    }
}
//##############################################################################
