#include <queue>
#include <map>
#include <unordered_map>
#include <Eigen/Dense>
#include "macros.h"
#include "io/FDTD_ListenShell.hpp"
#include "wavesolver/SimWorld.h"
#include "igl/read_triangle_mesh.h"
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"

//##############################################################################
// Static member definition
//##############################################################################
SimWorld::WorldRasterizer SimWorld::rasterizer;
Vector3Array ListeningUnit::microphones;
std::shared_ptr<FDTD_ListenShell<REAL>> ListeningUnit::refShell;

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
    // straight line connecting boxcenter and speaker
    if (listen->mode == ListeningUnit::MODE::DELAY_LINE)
    {
        const int N_mic = ListeningUnit::microphones.size();
        const REAL &cellSize = simulator->GetSolverSettings()->cellSize;
        if (listen->speakers.size() != N_mic)
            listen->speakers.resize(N_mic);
        for (int ii=0; ii<N_mic; ++ii)
        {
            const auto &mic = ListeningUnit::microphones.at(ii);
            auto &spk = listen->speakers.at(ii);
            Vector3d v = mic - BoundingBoxCenter();
            v.normalize();
            const REAL r = lowerRadiusBound
                - (0.5 + simulator->GetSolverSettings()->PML_width)*cellSize;
            v *= r;
            spk = BoundingBoxCenter() + v;
        }
    }
    else if (listen->mode == ListeningUnit::MODE::DELAY_LINE_2ND)
    {
        const int N_mic = ListeningUnit::microphones.size();
        const int N_spk = N_mic*2;
        if (listen->speakers.size() != N_spk)
            listen->speakers.resize(N_spk);
        const BoundingBox bbox = GetBoundingBox();
        const REAL dt = simulator->GetSolverSettings()->timeStepSize;
        const REAL h = simulator->GetSolverSettings()->cellSize;
        // for each delay line, find x_b at the box boundary
        for (int ii=0; ii<N_mic; ++ii)
        {
            const Vector3d &x_l = ListeningUnit::microphones.at(ii);
            const Vector3d  x_c = bbox.center();
            const Vector3d v = (x_l - x_c).normalized();
            REAL t_min = std::numeric_limits<REAL>::max();
            for (int dd=0; dd<3; ++dd)
            {
                if (fabs(v[dd]) < 1E-16)
                    continue;
                const REAL pad = (0.5+simulator->GetSolverSettings()->PML_width);
                const REAL l = bbox.axislength(dd)/2.0
                             - sqrt(3.0)*h*pad;
                t_min = std::min(t_min, l/fabs(v[dd]));
            }

            const REAL dr = sqrt(3.0)*h;
            const REAL ct = simulator->GetSolverSettings()->soundSpeed*dt;
            const Vector3d xb = x_c + v*t_min;
            const Vector3d x0 = x_c + v*(ct*std::floor((xb-x_c).length()/ct));
            const Vector3d x1 = x0  - v*(ct*std::ceil(dr/ct));
            // first write the outer, then the inner
            listen->speakers.at(ii*2  ) = x0;
            listen->speakers.at(ii*2+1) = x1;
        }
    }
    else if (listen->mode == ListeningUnit::MODE::SHELL)
    {
        const int N1 = listen->outShell->N();
        const int N2 = listen->innShell->N();
        if (listen->speakers.size() != (N1+N2))
            listen->speakers.resize((N1+N2));
        auto &outPoints = listen->outShell->Points();
        auto &innPoints = listen->innShell->Points();
        auto it = std::copy(outPoints.begin(), outPoints.end(), listen->speakers.begin());
        it = std::copy(innPoints.begin(), innPoints.end(), it);
        assert(it == listen->speakers.end());
    }
    return listen->speakers;
}

//##############################################################################
// Function Build
//##############################################################################
void SimWorld::
Build(ImpulseResponseParser_Ptr &parser, const uint &indTimeChunks)
{
    _objectCollections = std::make_shared<FDTD_Objects>();
    _simulatorSettings = std::make_shared<PML_WaveSolver_Settings>();

    // parse settings and construct objects in the scene
    parser->GetSolverSettings(_simulatorSettings);
    _simulatorSettings->indTimeChunks = indTimeChunks;
    parser->GetObjects(_simulatorSettings, _objectCollections);
    SimWorld::rasterizer.cellSize = _simulatorSettings->cellSize;
    // logic for Time Parallelism
    // This is the old logic for calculating startTime:
    // REAL simTime = _simulatorSettings->timeStepSize * _simulatorSettings->numberTimeSteps;
    // REAL chunkTime = simTime / _simulatorSettings->numTimeChunks;
    // REAL startTime = _simulatorSettings->fastForwardToEventTime + chunkTime * _simulatorSettings->indTimeChunks;
    // It turns out, there is a lot of machine error built up from stepping time forwards, so we
    //    need to calculate this the same way as the integrator does: via repeated adding
    int timeStepsPerChunk = _simulatorSettings->numberTimeSteps / _simulatorSettings->numTimeChunks;
    REAL startTime = _simulatorSettings->fastForwardToEventTime;
    for( int i = 0; i < timeStepsPerChunk * _simulatorSettings->indTimeChunks; i++ )
        startTime += _simulatorSettings->timeStepSize;
    _simulatorSettings->numberTimeSteps /= _simulatorSettings->numTimeChunks;

    if( _simulatorSettings->timeParallel )
    {
        // This is the old logic for calculating the stop Boundary Time:
        // _simulatorSettings->stopBoundaryAccTime = startTime + chunkTime;
        // This is the new logic, taking into account machine error:
        _simulatorSettings->stopBoundaryAccTime = startTime;
        for( int i = 0; i < timeStepsPerChunk; i++)
            _simulatorSettings->stopBoundaryAccTime += _simulatorSettings->timeStepSize;
        _simulatorSettings->numberTimeSteps += _simulatorSettings->overlapTime / _simulatorSettings->timeStepSize;
        _terminationMonitor.r = _simulatorSettings->earlyTerminationRatio;
    }

    // setup adaptive chunks if necessary
    if (_simulatorSettings->timeParallel)
    {
        auto acParam = parser->GetChunkPartitionParam();
        std::ifstream stream(acParam->outFile.c_str());
        if (acParam->adaptive && stream)
        {
            std::string line;
            while(std::getline(stream, line))
                if (line.at(0) != '#')
                    break;
            std::istringstream iss(line);
            uint N_c;
            iss >> N_c;
            assert(N_c > 0 && indTimeChunks < N_c);
            REAL T[2];
            int cIndex;
            for (uint ii=0; ii<N_c; ++ii)
            {
                std::getline(stream, line);
                std::istringstream iss(line);
                iss >> cIndex >> T[0] >> T[1];
                if (ii == indTimeChunks)
                    break;
            }
            // reset start and end time and number of steps
            REAL T_r[2]; // rounded time range for this chunk, within this range shader is on
            T_r[0] = _simulatorSettings->fastForwardToEventTime;
            while (T_r[0] < T[0])
                T_r[0] += _simulatorSettings->timeStepSize;
            T_r[1] = T_r[0];
            int steps = 0;
            while (T_r[1] < T[1])
            {
                T_r[1] += _simulatorSettings->timeStepSize;
                ++steps;
            }
            startTime = T_r[0];
            _simulatorSettings->stopBoundaryAccTime = T_r[1];
            _simulatorSettings->numberTimeSteps = steps
                                                + (int)(_simulatorSettings->overlapTime / _simulatorSettings->timeStepSize);

            std::cout << "Shader ON time-range for chunk " << indTimeChunks
                      << ": [" << T_r[0] << ", " << T_r[1] << "]\n";
            std::cout << "Solve time-range: [" << T_r[0] << ", "
                      << T_r[0] + (REAL)_simulatorSettings->numberTimeSteps*_simulatorSettings->timeStepSize << "]\n";
        }
        else
        {
            if (acParam->adaptive)
                std::cerr << "**WARNING** Adaptive start time is used but no partition config file found. Run chunk analysis first.\n";
        }
    }

    SetWorldTime(startTime);
    _state.startTime = startTime;
    _simulatorSettings->fastForwardToEventTime = startTime;

    // read and initialize animator if data exists
    if (_simulatorSettings->kinFileExists)
    {
        _objectCollections->InitializeAnimator(_simulatorSettings->objKinematicsMetadata);
    }
    else if (_simulatorSettings->rigidsimDataRead)
    {
        _objectCollections->InitializeAnimator(_simulatorSettings->fileDisplacement,
                                               _simulatorSettings->fileVelocity,
                                               _simulatorSettings->fileAcceleration);
    }
    UpdateObjectState(_state.time);

    std::map<ActiveSimUnit_Ptr, BoundingBox> candidateUnit;
    assert(_simulatorSettings->solverControlPolicy);
    if (_simulatorSettings->solverControlPolicy->type == "dynamic")
    {
        auto policy = std::dynamic_pointer_cast<Dynamic_Policy>(_simulatorSettings->solverControlPolicy);
        // assign a box to each object
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
            simUnit->objects->AddConstraints(_objectCollections);
            simUnit->simulator = std::make_shared<FDTD_AcousticSimulator>(simulatorID);
            simUnit->simulator->SetParser(parser);
            simUnit->simulator->SetSolverSettings(_simulatorSettings);
            simUnit->simulator->SetSceneObjects(simUnit->objects);
            simUnit->simulator->SetOwner(simUnit);

            auto meshPtr = obj->GetMeshPtr();
            const Vector3d meshCentroid_o = meshPtr->ComputeCentroid();
            const Vector3d meshCentroid_w = obj->ObjectToWorldPoint(meshCentroid_o);
            const Vector3d rastCentroid_w = SimWorld::rasterizer.cellCenter(
                                            SimWorld::rasterizer.rasterize(meshCentroid_w));
            const int divs = (int)std::ceil(
                    meshPtr->boundingSphereRadius(meshCentroid_o)/_simulatorSettings->cellSize
                    )*2 + policy->padding + (int)(_simulatorSettings->PML_width);
            const BoundingBox simUnitBox(
                    _simulatorSettings->cellSize, divs, rastCentroid_w);
            simUnit->divisions = divs;
            simUnit->listen = std::make_unique<ListeningUnit>();
            simUnit->lowerRadiusBound = simUnitBox.minlength()/2.0;
            simUnit->upperRadiusBound = simUnitBox.maxlength()/2.0;
            simUnit->unitID = simulatorID;

            // make sure this object is not inside some other existed unit
            bool createUnit = true;
            for (auto &pair : candidateUnit)
            {
                if (pair.second.isInside(meshCentroid_w))
                {
                    createUnit = false;
                    pair.first->objects->AddObject(std::stoi(obj->GetMeshName()), obj);
                    break;
                }
            }
            if (createUnit)
            {
                candidateUnit[simUnit] = simUnitBox;
                _simUnits.insert(std::move(simUnit));
            }
        }
    }
    else if (_simulatorSettings->solverControlPolicy->type == "markers")
    {
        auto policy =
            std::dynamic_pointer_cast<Markers_Policy>(_simulatorSettings->solverControlPolicy);
        ActiveSimUnit_Ptr simUnit = std::make_shared<ActiveSimUnit>();
        std::string *simulatorID = new std::string("0");
        simUnit->objects = _objectCollections;
        simUnit->simulator = std::make_shared<FDTD_AcousticSimulator>(simulatorID);
        simUnit->simulator->SetParser(parser);
        simUnit->simulator->SetSolverSettings(_simulatorSettings);
        simUnit->simulator->SetSceneObjects(simUnit->objects);
        simUnit->simulator->SetOwner(simUnit);

        const BoundingBox simUnitBox = FindSolverBBoxWithMarkers(policy, simUnit->divisions);
        simUnit->listen = std::make_unique<ListeningUnit>();
        simUnit->lowerRadiusBound = simUnitBox.minlength()/2.0;
        simUnit->upperRadiusBound = simUnitBox.maxlength()/2.0;
        simUnit->unitID = simulatorID;
        candidateUnit[simUnit] = simUnitBox;
        _simUnits.insert(std::move(simUnit));
    }
    else
    {
        auto policy = std::dynamic_pointer_cast<Static_Policy>(_simulatorSettings->solverControlPolicy);
        ActiveSimUnit_Ptr simUnit = std::make_shared<ActiveSimUnit>();
        std::string *simulatorID = new std::string("0");
        simUnit->objects = _objectCollections;
        simUnit->simulator = std::make_shared<FDTD_AcousticSimulator>(simulatorID);
        simUnit->simulator->SetParser(parser);
        simUnit->simulator->SetSolverSettings(_simulatorSettings);
        simUnit->simulator->SetSceneObjects(simUnit->objects);
        simUnit->simulator->SetOwner(simUnit);

        const BoundingBox simUnitBox(
                _simulatorSettings->cellSize,
                policy->cellDivisions,
                policy->domainCenter);
        simUnit->divisions = policy->cellDivisions;
        simUnit->listen = std::make_unique<ListeningUnit>();
        simUnit->lowerRadiusBound = simUnitBox.minlength()/2.0;
        simUnit->upperRadiusBound = simUnitBox.maxlength()/2.0;
        simUnit->unitID = simulatorID;
        candidateUnit[simUnit] = simUnitBox;
        _simUnits.insert(std::move(simUnit));
    }

    // initialize solver
    for (auto &pair : candidateUnit)
    {
        pair.first->simulator->InitializeSolver(pair.second, _simulatorSettings);
        pair.first->simulator->GetGrid().grid_id = pair.first->unitID;
    }

    // build listening shell

    // setup listening models
    char buffer[512];
    std::string filename("all_audio.dat");
    if(_simulatorSettings->timeParallel)
    {
        snprintf(buffer, 512, "%05d_all_audio.dat", _simulatorSettings->indTimeChunks);
        filename = std::string(buffer);
    }
    snprintf(buffer, 512, _simulatorSettings->outputPattern.c_str(),
             filename.c_str());
    int N_listen;
    switch ((*_simUnits.begin())->listen->mode)
    {
        case ListeningUnit::SHELL:
        {
            N_listen = (*_simUnits.begin())->listen->outShell->N() +
                       (*_simUnits.begin())->listen->innShell->N();
            break;
        }
        case ListeningUnit::DELAY_LINE:
        {
            N_listen = _simulatorSettings->listeningPoints.size();
            break;
        }
        case ListeningUnit::DELAY_LINE_2ND:
        {
            N_listen = _simulatorSettings->listeningPoints.size()*2;
            break;
        }
    }
    ListeningUnit::microphones = _simulatorSettings->listeningPoints;
    for (auto &p : _simUnits)
        p->UpdateSpeakers();
    AudioOutput::instance()->SetBufferSize(N_listen);
    AudioOutput::instance()->OpenStream(std::string(buffer));

    // write start time
    {
        char buffer[512];
        std::string filename = "start_time";
        if (_simulatorSettings->timeParallel)
        {
            snprintf(buffer, 512, "%05d_%s", _simulatorSettings->indTimeChunks,
                     filename.c_str());
            filename = std::string(buffer);
        }
        snprintf(buffer, 512, _simulatorSettings->outputPattern.c_str(),
                 filename.c_str());
        std::ofstream stream(buffer);
        if (stream)
        {
            stream << std::setprecision(18) << std::fixed
                   << _state.startTime << std::endl;
            stream.close();
        }
    }

    // write listening models metadata
    {
        for (auto &unit : _simUnits)
        {
            const std::string filename =
                unit->simulator->CompositeFilename("listening_model_metadata");
            std::ofstream stream(filename.c_str());
            if (stream)
            {
                std::string mode;
                switch (unit->listen->mode)
                {
                    case ListeningUnit::MODE::DELAY_LINE:
                    {
                        mode = "delay_line";
                        break;
                    }
                    case ListeningUnit::MODE::DELAY_LINE_2ND:
                    {
                        mode = "delay_line_2nd";
                        break;
                    }
                    case ListeningUnit::MODE::SHELL:
                    {
                        mode = "shell";
                        break;
                    }
                }

                stream << "# Listening model metadata\n"
                       << "# Listening model: " << mode << "\n"
                       << "# Number of far-field listening points: "
                       << ListeningUnit::microphones.size() << "\n"
                       << "# Number of near-field (in-box) sample points: "
                       << unit->listen->speakers.size() << "\n"
                       << "# Number of box center points: 1\n";
                stream << std::setprecision(16) << std::fixed;
                // microphone pos (far-field)
                stream << ListeningUnit::microphones.size() << "\n";
                for (const auto &v : ListeningUnit::microphones)
                    stream << v[0] << " " << v[1] << " " << v[2] << "\n";
                // speakers pos (in-box)
                stream << unit->listen->speakers.size() << "\n";
                for (const auto &v : unit->listen->speakers)
                    stream << v[0] << " " << v[1] << " " << v[2] << "\n";
                // sim box center
                stream << "1\n"
                       << unit->BoundingBoxCenter().x << " "
                       << unit->BoundingBoxCenter().y << " "
                       << unit->BoundingBoxCenter().z << "\n";
                stream.close();
            }
        }
    }
}

//##############################################################################
// Function UpdateObjectState
//##############################################################################
void SimWorld::
UpdateObjectState(const REAL &time)
{
    // update objects
    _objectCollections->SetObjectStates(time);

    if (_simulatorSettings->timeParallel && (time >= _simulatorSettings->stopBoundaryAccTime ) )
        if( !_objectCollections->AreEvalsDisabled() )
        {
            std::cout << "Disabling Boundary Evaluations!\n";
            _objectCollections->DisableAllEvals();
        }

    if (_simulatorSettings->solverControlPolicy->type != "dynamic")
        return;

    // logic for updating bbox and determine whether to move simbox
    for (auto &unit : _simUnits)
    {
        const auto &objects = unit->objects->GetRigidObjects();
        if (objects.size() == 0) continue;
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
        unit->unitCenter = newCenter;
    }
}

//##############################################################################
// Function StepWorld
//##############################################################################
bool SimWorld::
StepWorld()
{
    bool continueStepping = true;
    bool belowThres = false;
    // update simulation
    std::cout << "================ Step START ================\n";
    bool changed = _objectCollections->UpdateSourceTimes(_state.time);

    int count=0;
    for (auto &unit : _simUnits)
    {
        std::cout << "-------- unit " << count
                  << " START -------- \n";
        if (changed)
        {
            std::cout << "Mesh changed" << std::endl;
            unit->simulator->GetGrid().setMeshChanged();
        }

        continueStepping = (unit->simulator->RunForSteps(1) || continueStepping);
        // update speakers, interpolate pressure values and write to audio output
        const Vector3Array &spks = unit->UpdateSpeakers();
        FloatArray pressures(spks.size());
        Eigen::MatrixXd fetch;
        unit->simulator->GetSolver()->FetchPressureData(spks, fetch, -1);
        // scale the pressure, NOTE: ignoring the phase for now
        REAL pAbsPeak = 0.0;
        for (int ii=0; ii<spks.size(); ++ii)
        {
            const auto &spk = spks.at(ii);
            if (unit->listen->mode == ListeningUnit::DELAY_LINE)
            {
                const auto &mic = ListeningUnit::microphones.at(ii);
                pressures.at(ii) =
                    ListeningUnit::DelayLineScaling(spk, mic)* fetch(ii,0);
            }
            else if (unit->listen->mode == ListeningUnit::SHELL)
            {
                pressures.at(ii) = fetch(ii,0);
            }
            else if (unit->listen->mode == ListeningUnit::DELAY_LINE_2ND)
            {
                pressures.at(ii) = fetch(ii,0);
            }
            else
            {
                assert(false); // should not happen
            }
            pAbsPeak = std::max(pressures.at(ii), pAbsPeak);
        }

        // update pressure peak for the termination criterion
        const REAL simTime = unit->simulator->GetSimulationTime();
        _terminationMonitor.p_chunk = std::max(_terminationMonitor.p_chunk,
                                               pAbsPeak);
        if (_simulatorSettings->timeParallel &&
            simTime > _simulatorSettings->stopBoundaryAccTime)
        {
            if (pAbsPeak > _terminationMonitor.r * _terminationMonitor.p_chunk &&
                pAbsPeak > _terminationMonitor.absLowerBound)
            {
                _terminationMonitor.lastTimeAboveThres = simTime;
            }
            // pressure has fallen below threshold for the last L steps, stop stepping
            if ((simTime - _terminationMonitor.lastTimeAboveThres) > _terminationMonitor.windowSize)
            {
                belowThres = true;
                std::cout << "Chunk termination condition triggered, will terminate after this step";
                COUT_SDUMP(_terminationMonitor.lastTimeAboveThres);
                COUT_SDUMP(simTime);
                COUT_SDUMP(_terminationMonitor.windowSize);
                COUT_SDUMP(_terminationMonitor.absLowerBound);
                COUT_SDUMP(_terminationMonitor.r);
                COUT_SDUMP(_terminationMonitor.p_chunk);
                COUT_SDUMP(pAbsPeak);
            }
        }

        // update audio output buffer
        AudioOutput::instance()->AccumulateBuffer(pressures);
        std::cout << "-------- unit " << count
                  << " STOP -------- \n";
        count++;
    }
    AudioOutput::instance()->WriteAndResetBuffer();

    // update sim units topology
    CheckSimUnitBoundaries();

    // update time and object states
    _state.time += _simulatorSettings->timeStepSize;
#if 1
    UpdateObjectState(_state.time);
#else
    // simple translation
    //auto objects = _objectCollections->GetRigidSoundObjects();
    //for (auto &m : objects)
    //{
    //    m.second->ApplyTranslation(0.0, -2.0E-5, 0.0);
    //    std::cout << "center = " << m.second->GetBBox().Center() <<std::endl;
    //}
    // monopole vibration
    //const REAL omega = 2.0*M_PI*1500.0;
    //const REAL r0 = 0.05;
    //const REAL dr = 0.01;
    //auto objects = _objectCollections->GetRigidObjects();
    //for (auto &m : objects)
    //{
    //    const REAL sdot = -dr*omega/r0*cos(omega*_state.time);
    //    std::cout << "scaling = " << 1.0+sdot*_simulatorSettings->timeStepSize << std::endl;
    //    std::cout << "center = " << m.second->GetBBox().Center() <<std::endl;
    //    m.second->ApplyScale(1.0+sdot*_simulatorSettings->timeStepSize);
    //}
    // dipole vibration
    const REAL omega = 2.0*M_PI*1500.0;
    const REAL scale = omega*omega*0.01;
    auto objects = _objectCollections->GetRigidObjects();
    for (auto &m : objects)
    {
        m.second->ApplyTranslation(0., -1./omega*cos(omega*_state.time)*_simulatorSettings->timeStepSize*scale, 0.);
        std::cout << "center = " << m.second->GetBBox().Center() <<std::endl;
    }
#endif
    std::cout << "================ Step STOP ================\n";

    return continueStepping && !belowThres;
}

//##############################################################################
// Function CheckSimUnitBoundaries
//##############################################################################
bool SimWorld::
CheckSimUnitBoundaries()
{
    // clear all stored interfaces
    for (auto &unit : _simUnits) // stored in Mac_Grid
        unit->simulator->GetGrid().ClearBoundaryInterface();

    // establish interfaces
    for (auto it_a=_simUnits.begin(); it_a!=_simUnits.end(); ++it_a)
    {
        auto it_b = it_a;
        std::advance(it_b, 1);
        for (; it_b!=_simUnits.end(); ++it_b)
        {
            const Vector3d centerDiff =
                (*it_a)->BoundingBoxCenter() - (*it_b)->BoundingBoxCenter();
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
                std::cout << "locking sim box move\n";
                (*it_a)->simulator->GetSolver()->SetSimBoxForceFix(true);
                (*it_b)->simulator->GetSolver()->SetSimBoxForceFix(true);
                auto interface = std::make_shared<BoundaryInterface>(
                        (*it_a), (*it_b), GetWorldTime(), interfaceDirection);
                bool exist = false;
                for (auto &exist_int : _interfaces)
                    if (exist_int->Equal(*interface))
                        exist = true;
                if (!exist) // first time discover this interface
                {
                    (*it_a)->objects->Join((*it_b)->objects);
                    (*it_b)->objects->Join((*it_a)->objects);
#ifdef USE_MERGE_SIM_BOXES
                    _interfaces.push_back(interface);
#endif
                }
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
        grid_a.AddBoundaryInterface(interface);
        grid_b.AddBoundaryInterface(interface);
        if (interface->initialized)
            continue;

        const bool a_on_top_of_b =
            ((unit_a->BoundingBoxCenter()
             -unit_b->BoundingBoxCenter())[dir] > 0);
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
        interface->initialized = true;
    }
    return (_interfaces.size() > 0);
}

//##############################################################################
// Function PreviewStepping
//##############################################################################
void SimWorld::
PreviewStepping(const uint &previewSpeed)
{
    // update simulation
    int count=0;
    for (auto &unit : _simUnits)
    {
        unit->simulator->PreviewStepping(previewSpeed);
        count++;
    }

    // update sim units topology
    CheckSimUnitBoundaries();

    // update time and object states
    _state.time += _simulatorSettings->timeStepSize*previewSpeed;
    UpdateObjectState(_state.time);
    _objectCollections->UpdateSourceTimes(_state.time);
    std::cout << "================ Step STOP ================\n";
}

//##############################################################################
// Function ResetStartTime
//##############################################################################
void SimWorld::
ResetStartTime(const REAL &startTime)
{
    for (auto &unit : _simUnits)
        unit->simulator->ResetStartTime(startTime);
}

//##############################################################################
// Function ClearAllSources
//##############################################################################
void SimWorld::
ClearAllSources()
{
    // pressure source
    _objectCollections->GetPressureSources().clear();
    // object sources
    auto &objects = _objectCollections->GetRigidObjects();
    for (auto &obj : objects)
        obj.second->ClearVibrationalSources();
}

//##############################################################################
// Function ClearAllSources
//##############################################################################
void SimWorld::
PostRunWrite()
{
    // write chunk time range
    {
        char buffer[512];
        std::string filename = "time_range";
        if (_simulatorSettings->timeParallel)
        {
            snprintf(buffer, 512, "%05d_%s", _simulatorSettings->indTimeChunks,
                     filename.c_str());
            filename = std::string(buffer);
        }
        snprintf(buffer, 512, _simulatorSettings->outputPattern.c_str(),
                 filename.c_str());
        std::ofstream stream(buffer);
        if (stream)
        {
            stream << std::setprecision(18) << std::fixed
                   << _state.startTime << std::endl
                   << _simulatorSettings->stopBoundaryAccTime << std::endl
                   << _state.time << std::endl;
            stream.close();
        }
    }

    // write a txt file to label it finish
    {
        char buffer[512];
        std::string filename = "sim-done";
        if (_simulatorSettings->timeParallel)
        {
            snprintf(buffer, 512, "%05d_%s", _simulatorSettings->indTimeChunks,
                     filename.c_str());
            filename = std::string(buffer);
        }
        snprintf(buffer, 512, _simulatorSettings->outputPattern.c_str(),
                 filename.c_str());
        std::ofstream stream(buffer);
        stream.close();
    }
}

//##############################################################################
// Function RunChunksAnalysis
//##############################################################################
void SimWorld::
RunChunksAnalysis(const ChunkPartitionParam_Ptr &param, const int N_t)
{
    std::cout << "\n\n==========================================================\n";
    std::cout << "Running chunks analysis ...\n";
    std::cout << "==========================================================\n";
    // first print some info
    std::cout << " World start time     : " << _state.time             << "\n";
    std::cout << " Number of total steps: " << N_t                     << "\n";
    std::cout << " Number of objects    : " << _objectCollections->N() << "\n";
    const auto &objs = _objectCollections->GetRigidObjects();;
    for (const auto &m : objs)
    {
        const auto &obj = m.second;
        std::cout << "  Object " << obj->GetMeshName()
                  << " has " << obj->NumberVibrationalSources()
                  << " vibrational shader sources\n";
    }
    std::cout << " Chunk Partitioning Parameters: \n"
              << "  adaptive = " << param->adaptive << "\n"
              << "  L_z      = " << param->L_z      << "\n"
              << "  N_0      = " << param->N_0      << "\n"
              << "  N_maxc   = " << param->N_maxc   << "\n"
              << "  out_file = " << param->outFile  << "\n";
    // compute zero masks
    std::cout << "\nComputing zero masks ...\n";
    auto &ws = _simulatorSettings;
    const REAL dt = ws->timeStepSize;
    std::vector<REAL> q;
    std::vector<Eigen::Vector2d> Z, Z_bar;
    const REAL T[2] = {_state.time,
                       _state.time + (REAL)N_t *dt};
    for (int ii=0; ii<N_t; ++ii)
    {
        const REAL t = T[0] + (REAL)ii*dt;
        if (_objectCollections->ShaderIsZero(t))
        {
            q.push_back(t);
        }
        else
        {
            if (q.size() != 0 && (t - q[0]) >= param->L_z)
            {
                Z.push_back({q[0], t});
            }
            q.clear();
        }
    }
    if (q.size() != 0)
    {
        Z.push_back({q[0], T[1]});
    }

    std::cout << " Z = \n";
    for (const auto &z : Z)
        std::cout << "    " << z.transpose() << std::endl;

    // compute complement
    std::cout << "\nComputing complement masks ...\n";
    if (Z.size() == 0)
    {
        Z_bar.push_back({T[0], T[1]});
    }
    else
    {
        if (!EQUAL_FLOATS(Z[0][0], T[0]))
            Z_bar.push_back({T[0], Z[0][0]});
        for (int ii=0; ii<Z.size()-1; ++ii)
            Z_bar.push_back({Z[ii][1], Z[ii+1][0]});
        if (!EQUAL_FLOATS(Z[Z.size()-1][1], T[1]))
            Z_bar.push_back({Z[Z.size()-1][1], T[1]});
    }

    std::cout << " Z_bar = \n";
    for (const auto &z : Z_bar)
        std::cout << "    " << z.transpose() << std::endl;

    // helper function to compute intersections
    auto Intersect = [](const Eigen::Vector2d &I,
                        const std::vector<Eigen::Vector2d> &Z,
                        std::vector<Eigen::Vector2d> &out)
    {
        out.clear();
        for (const Eigen::Vector2d &z : Z)
        {
            if ((z[0] < I[0] && z[1] < I[0]) ||
                (z[0] > I[1] && z[1] > I[1]))
                continue;
            out.push_back({std::max(I[0], z[0]),
                           std::min(I[1], z[1])});
        }
    };
    auto PrintRange = [](const Eigen::Vector2d &s)
    {
        char buf[512];
        snprintf(buf, 512, "%.6f, %.6f", s[0], s[1]);
        return std::string(buf);
    };

    // initialize uniform chunks and shrink them
    std::cout << "\nInitializing uniform chunks and shrink the bounds ...\n";
    std::map<int, Eigen::Vector2d> S;
    const REAL h_0 = (T[1] - T[0])/(REAL)(param->N_0);
    for (int ii=0; ii<param->N_0; ++ii)
    {
        const Eigen::Vector2d I_i = {T[0] + (REAL)(ii+0)*h_0,
                                     T[0] + (REAL)(ii+1)*h_0};
        std::vector<Eigen::Vector2d> z_bar;
        Intersect(I_i, Z_bar, z_bar);
        if (z_bar.size() > 0)
        {
            S[ii] = {std::max(I_i[0], z_bar[0             ][0]),
                     std::min(I_i[1], z_bar[z_bar.size()-1][1])};
        }
    }

    std::cout << " (vanilla) S = \n";
    for (const auto &s : S)
        std::cout << "    " << s.first << ": "
                  << PrintRange(s.second) << std::endl;

    // greedily pick the largest zero-mask and remove them from chunks until
    // reaches limits
    std::cout << "\nSubdivide chunks starting from the largest zero part ...\n";
    std::sort(Z.begin(), Z.end(),
              [](const Eigen::Vector2d &a, const Eigen::Vector2d &b)
              {return (a[1]-a[0]) > (b[1]-b[0]);});

    std::cout << " sorted Z = \n";
    for (const auto &z : Z)
        std::cout << "    " << z[1]-z[0] << ": " << PrintRange(z) << std::endl;

    int i = S.size();
    int largest_idx = 0;
    while (i < param->N_maxc && largest_idx < Z.size())
    {
        const auto &z_q = Z.at(largest_idx);
        ++largest_idx;
        if (z_q[1]-z_q[0] < param->L_z)
            break;
        std::vector<Eigen::Vector2d> inter;
        std::vector<Eigen::Vector2d> z_qs = {z_q};
        std::map<int, Eigen::Vector2d> S_freeze = S;
        for (const auto &s : S_freeze)
        {
            Intersect(s.second, z_qs, inter);
            for (const Eigen::Vector2d z : inter)
            {
                if (z[1]-z[0] > param->L_z)
                {
                    assert(S.find(i) == S.end());
                    S[i] = {S_freeze[s.first][0], z[0]};
                    ++i;
                    assert(S.find(i) == S.end());
                    S[i] = {z[1], S_freeze[s.first][1]};
                    ++i;
                    S.erase(s.first);
                }
            }
        }
    }

    std::cout << " S = \n";
    for (const auto &s : S)
        std::cout << "    " << s.first << ": "
                  << PrintRange(s.second) << std::endl;

    // writing it to output file
    std::cout << "\nWriting chunks to file: " << param->outFile << std::endl;
    {
        std::ofstream stream(param->outFile.c_str());
        if (stream)
        {
            stream << "#Chunk partitioning config file for wavesolver time-parallelization\n";
            stream << "#Parameters:\n"
                   << "#  adaptive = " << param->adaptive << "\n"
                   << "#  L_z      = " << param->L_z      << "\n"
                   << "#  N_0      = " << param->N_0      << "\n"
                   << "#  N_maxc   = " << param->N_maxc   << "\n"
                   << "#  out_file = " << param->outFile  << "\n"
                   << "#\n"
                   << "#Partitions:\n";
            stream << S.size() << "\n";
            int count=0;
            for (const auto &s : S)
            {
                stream << std::setprecision(16) << std::fixed
                       << count << " " << s.second[0] << " " << s.second[1] << "\n";
                ++count;
            }
        }
    }
}

//##############################################################################
// Function FindSolverBBoxWithMarkers
//##############################################################################
BoundingBox SimWorld::
FindSolverBBoxWithMarkers(std::shared_ptr<Markers_Policy> policy, int &divs)
{
    std::cout << "Find Solver bbox with markers ...\n";
    using namespace boost::filesystem;
    // iterate and find all marker files
    path p(policy->markersDir.c_str());
    std::vector<path> filenames;
    for (directory_iterator it(p); it!=directory_iterator(); ++it)
        if (it->path().extension().string() == ".obj")
            filenames.push_back(it->path());
    auto ParseFileID = [](const path &a)
    {
        std::vector<std::string> tokens;
        boost::split(tokens, a.stem().string(), [](char c){return c == '_';});
        return std::stoi(tokens.at(1));
    };
    std::sort(filenames.begin(), filenames.end(),
            [&](const path &a, const path &b){
                return ParseFileID(a) < ParseFileID(b);});

    // bound the solver box
    const REAL &h = _simulatorSettings->cellSize;
    const REAL &H = policy->boxSizeMarker;
    const REAL H_over_2 = H/2.0;
    Vector3d minBound(F_MAX, F_MAX, F_MAX);
    Vector3d maxBound(F_LOW, F_LOW, F_LOW);
    for (const auto &f : filenames)
    {
        const int idx = ParseFileID(f);
        const REAL t = policy->startTime + (REAL)idx/policy->frameRate;
        if (t > _state.startTime &&
            t < _simulatorSettings->stopBoundaryAccTime)
        {
            // compute centroid of marker
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            igl::read_triangle_mesh(f.string().c_str(), V, F);
            const Eigen::Vector3d C = V.colwise().sum() / (REAL)V.rows();

            // min/max the box
            for (int dd=0; dd<3; ++dd)
            {
                minBound[dd] = std::min(minBound[dd], C[dd]-H_over_2);
                maxBound[dd] = std::max(maxBound[dd], C[dd]+H_over_2);
            }
        }
    }
    const Vector3d center = (maxBound + minBound) / 2.0;
    const Vector3d diff = (maxBound - minBound);
    REAL l = 0.0;
    for (int dd=0; dd<3; ++dd)
        l = std::max(l, diff[dd]);
    divs = (int)(l/h);
    minBound = center - l/2.0;
    maxBound = center + l/2.0;
    BoundingBox bbox;
    bbox.setMinBound(minBound);
    bbox.setMaxBound(maxBound);
    std::cout << " ... Done. Box = " << minBound << " " << maxBound << std::endl;
    return bbox;
}

//##############################################################################
