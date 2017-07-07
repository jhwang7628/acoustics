#include <unordered_map>
#include <Eigen/Dense>
#include "macros.h"
#include "geometry/BoundingBox.h" 
#include "wavesolver/SimWorld.h"

//##############################################################################
// Static member definition
//##############################################################################
SimWorld::WorldRasterizer SimWorld::rasterizer; 
Vector3Array ListeningUnit::microphones; 

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
// Function UpdateSpeakers
//##############################################################################
Vector3Array &ActiveSimUnit::
UpdateSpeakers()
{
    if (!listen)
        return listen->speakers; 

    // straight line connecting boxcenter and speaker
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
    return listen->speakers; 
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
                )*2 + 8;
        const BoundingBox simUnitBox(
                _simulatorSettings->cellSize, divs, rastCentroid_w); 
        simUnit->simulator->InitializeSolver(simUnitBox, _simulatorSettings); 
        simUnit->divisions = divs; 
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
    ListeningUnit::microphones = _simulatorSettings->listeningPoints; 
    AudioOutput::instance()->SetBufferSize(ListeningUnit::microphones.size()); 
    AudioOutput::instance()->OpenStream(std::string(buffer)); 
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
        // update speakers, interpolate pressure values and write to audio output
        const Vector3Array &spks = unit->UpdateSpeakers();
        FloatArray pressures(spks.size()); 
        Eigen::MatrixXd fetch; 
        unit->simulator->GetSolver()->FetchPressureData(spks, fetch, -1);
        // scale the pressure, NOTE: ignoring the phase for now
        for (int ii=0; ii<spks.size(); ++ii)
        {
            const auto &mic = ListeningUnit::microphones.at(ii); 
            const auto &spk = spks.at(ii); 
            pressures.at(ii) = 
                ListeningUnit::DelayLineScaling(spk, mic)* fetch(ii,0);
        }
        AudioOutput::instance()->AccumulateBuffer(pressures); 
    }
    AudioOutput::instance()->WriteAndResetBuffer(); 

    // update sim units topology
    CheckSimUnitBoundaries(); 

    // update time and object states
    _state.time += _simulatorSettings->timeStepSize; 
    UpdateObjectState(_state.time); 

    return continueStepping;
}

//##############################################################################
// Function CheckSimUnitBoundaries
//##############################################################################
bool SimWorld::
CheckSimUnitBoundaries()
{
    // establish interfaces
    _interfaces.clear();
    for (auto it_a=_simUnits.begin(); it_a!=_simUnits.end(); ++it_a)
    {
        auto it_b = it_a; 
        std::advance(it_b, 1); 
        for (; it_b!=_simUnits.end(); ++it_b)
        {
            const Vector3d centerDiff = 
                (*it_a)->boxCenter - (*it_b)->boxCenter; 
            const REAL maxDiff = std::max(
                    std::max(std::abs(centerDiff[0]), std::abs(centerDiff[1])),
                    std::abs(centerDiff[2])); 
            const REAL thres = (REAL)((*it_a)->divisions + (*it_b)->divisions)
                             / (REAL)2 * _simulatorSettings->cellSize; 
            if (maxDiff < thres || EQUAL_FLOATS(maxDiff, thres)) // <=
            {
                int interfaceDirection = -1; 
                for (int ii=0; ii<3; ++ii) 
                    if (EQUAL_FLOATS(std::abs(centerDiff[ii]), thres))
                        interfaceDirection = ii; 
                std::cout << "found interface at direction: " 
                          << interfaceDirection << std::endl;
                auto interface = std::make_shared<BoundaryInterface>(
                        (*it_a), (*it_b), GetWorldTime(), interfaceDirection); 
                _interfaces.push_back(interface); 
            }
        }
    }

    // establish cell neighbours
    for (auto &interface : _interfaces)
    {
        auto unit_a = interface->GetSimUnit_a(); 
        auto unit_b = interface->GetSimUnit_b(); 
        auto &grid_a = unit_a->simulator->GetGrid(); 
        auto &grid_b = unit_b->simulator->GetGrid(); 
        const int dir = interface->GetDirection(); 
        const bool a_on_top_of_b = 
            ((unit_a->boxCenter - unit_b->boxCenter)[dir] > 0);

        std::vector<int> bdIndices_a, bdIndices_b; 
        std::vector<Vector3d> bdPositions_a, bdPositions_b; 
        if (a_on_top_of_b) // grab pos b and neg a
        {
            grid_a.GetAllBoundaryCells(dir, -1, bdIndices_a, bdPositions_a); 
            grid_b.GetAllBoundaryCells(dir,  1, bdIndices_b, bdPositions_b); 
        }
        else  // grab neg b and pos a
        {
            grid_a.GetAllBoundaryCells(dir,  1, bdIndices_a, bdPositions_a); 
            grid_b.GetAllBoundaryCells(dir, -1, bdIndices_b, bdPositions_b); 
        }

        struct ProjCounter
        {
            int count = 0; 
            int otherUnitCellIdx = -1; 
        };
        const bool a_smaller_than_b = (unit_a->divisions < unit_b->divisions); 
        std::unordered_map<int,ProjCounter> counter; 
        const REAL offset = _simulatorSettings->cellSize;
        if (a_smaller_than_b)
        {
            for (const auto &ind : bdIndices_b)
                counter[ind].count = 1; 
            // offset a's boundary to b and see what cell its in
            for (int ii=0; ii<bdPositions_a.size(); ++ii)
            {
                const auto &pos = bdPositions_a.at(ii); 
                Vector3d posProj = pos;
                posProj[dir] += (a_on_top_of_b ? -1.0 : 1.0)*offset; 
                const int cellProj = grid_b.InPressureCell(posProj); 
                if (counter.find(cellProj) != counter.end())
                {
                    counter[cellProj].count += 1; 
                    counter[cellProj].otherUnitCellIdx = bdIndices_a.at(ii);
                }
                else 
                {
                    counter[cellProj].count = 1; 
                    counter[cellProj].otherUnitCellIdx = bdIndices_a.at(ii);
                }
            }
            for (const auto &m : counter) 
            {
                assert(m.second.count <= 2); 
                if (m.second.count == 2) 
                {
                    auto cellPair = std::make_pair(
                        m.second.otherUnitCellIdx, 
                        m.first); 
                    interface->AddCellPair(cellPair); 
                }
            }
        }
        else
        {
            for (const auto &ind : bdIndices_a)
                counter[ind].count = 1; 
            // offset b's boundary to a and see what cell its in
            for (int ii=0; ii<bdPositions_b.size(); ++ii)
            {
                const auto &pos = bdPositions_b.at(ii); 
                Vector3d posProj = pos;
                posProj[dir] += (a_on_top_of_b ? 1.0 : -1.0)*offset; 
                const int cellProj = grid_a.InPressureCell(posProj); 
                if (counter.find(cellProj) != counter.end())
                {
                    counter[cellProj].count += 1; 
                    counter[cellProj].otherUnitCellIdx = bdIndices_b.at(ii); 
                }
                else 
                {
                    counter[cellProj].count = 1; 
                    counter[cellProj].otherUnitCellIdx = bdIndices_b.at(ii); 
                }
            }
            for (const auto &m : counter) 
            {
                assert(m.second.count <= 2); 
                if (m.second.count == 2) 
                {
                    auto cellPair = std::make_pair(
                        m.first,
                        m.second.otherUnitCellIdx); 
                    interface->AddCellPair(cellPair); 
                }
            }
        }
    }
    return false; 
}

//##############################################################################
// Function PreviewStepping
//##############################################################################
void SimWorld::
PreviewStepping(const int &previewSpeed) // TODO!
{
}

//##############################################################################
