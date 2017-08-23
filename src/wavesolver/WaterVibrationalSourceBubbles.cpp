#include <wavesolver/WaterVibrationalSourceBubbles.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <sndgen/WavReader.hpp>
#include <utils/STL_Wrapper.h>
#include <fstream>

#include "bubbles/bubbleForcing.hpp"
#include "bubbles/FileInput.hpp"
#include "bubbles/ODEInt.hpp"

//##############################################################################
//##############################################################################
WaterVibrationalSourceBubbles::
WaterVibrationalSourceBubbles(RigidObjectPtr owner, const std::string &dataDir)
    : VibrationalSource(owner), _surfaceMesh(owner->GetMeshPtr())
{
    Initialize(dataDir);
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    // TODO: make this explicit
    // step if necessary
#pragma omp critical
    while (std::fabs(time - _curTime) > 1e-12)
    {
        step(time);
    }

    // transform the sample point to object frame
    //const Eigen::Vector3d samplePointObject_e = _modelingTransformInverse * Eigen::Vector3d(position.x, position.y, position.z);
    //const Vector3d samplePointObject(samplePointObject_e[0], samplePointObject_e[1], samplePointObject_e[2]);

    const Vector3d samplePointObject = _owner->WorldToObjectPoint(position);

    // use spatial partitioning
    int closestTriangle;
    if (hintTriangle < 0)
    {
#ifdef USE_OPENMP
#pragma omp critical
#endif
        _surfaceMesh->FindNearestTriangle(samplePointObject, closestTriangle);
    }
    else
    {
        std::vector<int> neighbours;
        std::dynamic_pointer_cast<TriangleMeshGraph<REAL>>(_surfaceMesh)->FindKNearestTrianglesGraph(1, samplePointObject, 1, hintTriangle, neighbours);
        closestTriangle = neighbours.at(0);
    }

    return _projectedAccel(closestTriangle);
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    return Evaluate(_surfaceMesh->vertex(vertexID), vertexNormal, time, -1);

    //throw std::runtime_error("**ERROR** Cannot sample water vibrational source using vertexID");
    //return 0.0;
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
Initialize(const std::string &dataDir)
{
    std::cout << "Initialize WaterVibrationalSourceBubbles from directory: " << dataDir << std::endl;
    _fileInfo = parseFileNames(dataDir);

    // DEBUGGING
    for (auto iter = _fileInfo.begin(); iter != _fileInfo.end(); ++iter)
    {
        std::cout << iter->first << ":\n"
                  << "   " << iter->second.meshFile << "\n"
                  << "   " << iter->second.datFile << "\n"
                  << "   " << iter->second.freqFile << std::endl;
    }

    _curTime = -1;
    _t1 = _t2 = -1;
    _dt = 1.0 / 192000.;

    FreqType ft = CAPACITANCE;
    std::string infoFile = dataDir + std::string("/bemOutput/oscillators/trackedBubInfo.txt");
    parseConfigFile(infoFile, ft);
    makeOscillators(_bubbles);
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
step(REAL time)
{
    // First load new solution data if necessary
    if (time >= _t2)
    {
        // Advance
        _t1 = _t2;
        _m1 = _m2;
        _v1 = _v2;
        _kd1 = _kd2;
        _fullKd1 = _fullKd2;
        _b1 = _b2;

        // New t2
        auto iter = _fileInfo.upper_bound(time);
        if (iter == _fileInfo.end())
        {
            // Past the last solution data time
            // TODO: is this the best solution?
            _t2 = 50000;
            _m2 = _m1;
            _v2 = _v1;
            _kd2 = _kd1;
            _fullKd2 = _fullKd1;
            _b2 = _b1;
        }
        else
        {
            _t2 = iter->first;
            _m2.loadGmsh(iter->second.meshFile);
            _b2 = parseFreqFile(iter->second.freqFile);
            _v2 = loadSurfaceDatFile(_b2,
                                     iter->second.datFile,
                                     _m2);

            _kd2.reset(new PointKDTree(_m2.m_surfTriCenters.data(), _m2.m_surfTriCenters.size(), false));
            _fullKd2.reset(new PointKDTree(_m2.m_allTriCenters.data(), _m2.m_allTriCenters.size(), false));
        }

        if (_t1 < 0)
        {
            _t1 = -1;
            _m1 = _m2;
            _v1 = _v2;
            _kd1 = _kd2;
            _fullKd1 = _fullKd2;
            _b1 = _b2;
        }
    }

    // Update all oscillators to the correct time (one timestep after this time)
    updateOscillators(time);

    // Compute velocity before and after current time step
    computeVelocities(time);

    // Compute acceleration
    _accel = (_velT2 - _velT1) / (2 * _dt);

    // Project to surface
    // This step is only necessary until the wavesolver handles deforming geometry
    projectToSurface();

    _curTime = time;
    //_curTime += _dt;

    std::cout << "Bubbles time: " << _curTime << ", t1: " << _t1 << ", t2: " << _t2 << std::endl;
    std::cout << "vel1 max: " << _velT1.cwiseAbs().maxCoeff() << ", vel2 max: " << _velT2.cwiseAbs().maxCoeff() << std::endl;
    std::cout << "Accel max: " << _accel.cwiseAbs().maxCoeff() << ", accel l2: " << _accel.norm() << std::endl;
    std::cout << "Projected accel max: " << _projectedAccel.cwiseAbs().maxCoeff() << std::endl;
    std::cout << "Accel inf/nan: " << ! _projectedAccel.array().isFinite().all() << std::endl;
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
updateOscillators(REAL time)
{
    using namespace MathUtils;
    using namespace Eigen;

    // Temp debugging
    using namespace std;
    static bool first = true;
    static vector<ofstream> oscOutputs;

    if (first)
    {
        for (int i = 0; i < _oscillators.size(); ++i)
        {
            ostringstream os;
            os << "oscillator_" << i << ".txt";
            oscOutputs.emplace_back(os.str().c_str());
        }

        first = false;
    }

    RK4<Vector2d> integrator;
    typedef BubbleOscillator<Vector2d> DS;

    for (int i = 0; i < _oscillators.size(); ++i)
    {
        Oscillator &osc = _oscillators[i];

        // Skip if there are not enough frequency solves
        if (osc.m_frequencies.getData().rows() < 1) continue;

        // If this oscillator is not active yet, skip it
        if (!osc.isActive(time)) continue;

        DS bubOsc(osc, 0);

        // Can step now
        if (osc.m_currentTime < 0)
        {
            // Initializing this oscillator
            osc.m_currentTime = time;

            osc.m_state = integrator.step(osc.m_state,
                                          time,
                                          _dt,
                                          bubOsc);

            osc.m_lastVals.col(2) = osc.m_state;

            oscOutputs[i] << osc.m_state(1) << endl;
        }

        while (osc.m_currentTime < time)
        {
            osc.m_state = integrator.step(osc.m_state,
                                          time,
                                          _dt,
                                          bubOsc);

            oscOutputs[i] << osc.m_state(1) << endl;

            osc.m_lastVals.col(0) = osc.m_lastVals.col(1);
            osc.m_lastVals.col(1) = osc.m_lastVals.col(2);
            osc.m_lastVals.col(2) = osc.m_state;

            osc.m_currentTime += _dt;
        }
    }
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
parseConfigFile(const std::string &infoFile, FreqType fType)
{
    _bubbles.clear();

    std::ifstream in(infoFile.c_str());

    while (in.good())
    {
        std::pair<int, Bubble> bub = parseBubbleInfo(in, fType);

        if (bub.first >= 0)
        {
            _bubbles.insert(bub);
        }

        //std::cout << "\rLoaded " << bub.first << "    ";
        //std::cout.flush();
    }

    std::cout << "Loaded " << _bubbles.size() << " bubbles" << std::endl;
}

// Returns the (angular) Minnaert frequency
static double
minnaertFreq(double radius)
{
    return std::sqrt(3 * GAMMA * ATM - 2 * SIGMA / radius) / (radius * std::sqrt(RHO_WATER));
}

//##############################################################################
//##############################################################################
std::pair<int, Bubble> WaterVibrationalSourceBubbles::
parseBubbleInfo (std::ifstream &in,
                 FreqType fType)
{
    std::pair<int, Bubble> bubAndNum;
    Bubble &bub = bubAndNum.second;

    std::vector<double> times, freqs, transfs, pressures;
    std::vector<int> bubNums;

    double t, freq;
    std::complex<double> trans;
    double pressure;

    std::string line;

    // First line is Bub <number>
    std::getline(in, line);

    if (line.empty())
    {
        bubAndNum.first = -1;
        return bubAndNum;
    }

    {
        std::istringstream is(line.substr(4));
        is >> bubAndNum.first;
    }

    // Second line is the start info
    std::getline(in, line);

    char startType = line.at(7);
    std::istringstream is(line.substr(8));
    is >> bub.m_startTime;
    int bubNum;
    while (is >> bubNum)
    {
        bub.m_prevBubbles.push_back(bubNum);
    }

    bub.m_startType = Bubble::parseEventType(startType);

    // Next line is radius info
    std::getline(in, line);
    line = line.substr(6);
    bub.m_radius = std::stod(line);

    while ( in.peek() != 'E')
    {
        std::getline(in, line);
        std::istringstream newis(line);
        newis >> t >> freq >> trans >> pressure >> bubNum;

        times.push_back(t);

        if (fType == CAPACITANCE)
        {
            freqs.push_back(2 * M_PI * freq);
        }
        else
        {
            freqs.push_back(minnaertFreq(bub.m_radius));
        }

        transfs.push_back(std::abs(trans));
        pressures.push_back(pressure);
        bubNums.push_back(bubNum);
    };

    // Last line is end data
    std::getline(in, line);
    char endType = line[5];
    line = line.substr(6);
    std::istringstream newis(line);
    newis >> bub.m_endTime;
    while (newis >> bubNum)
    {
        bub.m_nextBubbles.push_back(bubNum);
    }

    bub.m_endType = Bubble::parseEventType(endType);

    // TODO: set a Minnaert freq for the bubble if no solve data present?
    // How to handle transfer?
    if (!times.empty())
    {
        // Update freqs for start and end
        if (bub.m_startTime < times.at(0))
        {
            //times.insert(times.begin(), bub.m_startTime);
            //freqs.insert(freqs.begin(), freqs.at(0));
            //transfs.insert(transfs.begin(), transfs.at(0));
            //pressures.insert(pressures.begin(), pressures.at(0));
        }

        if (bub.m_endTime > times.back())
        {
            //times.push_back(bub.m_endTime);
            //freqs.push_back(freqs.back());
            //transfs.push_back(transfs.back());
            //pressures.push_back(pressures.back());
        }

        bub.m_times = times;
        bub.m_frequencies = freqs;
        bub.m_transfer = transfs;
        bub.m_pressure = pressures;
        bub.m_nums = bubNums;
    }

    return bubAndNum;

    //std::cout << "Radius: " << radius << std::endl
              //<< "freqs: \n" << frequencies.getData() << std::endl
              //<< "transfer: \n" << transfer.getData() << std::endl;
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
makeOscillators(const std::map<int, Bubble> &singleBubbles)
{
#define CHAIN_BUBBLES
    _oscillators.clear();
    std::set<int> used;

    int numFiltered = 0;

    for (const auto& bubPair : singleBubbles)
    {
        //std::cout << "\rTracking " << bubPair.first << " of " << singleBubbles.size() << "       ";
        //std::cout.flush();

        if (!used.count(bubPair.first))
        {
            const Bubble* curBub = &bubPair.second;
            int curBubNum = bubPair.first;

            Oscillator osc;
            osc.m_startTime = curBub->m_startTime;
            Bubble::EventType oscStartType = curBub->m_startType;
            double firstLength = curBub->m_endTime - curBub->m_startTime;

            //if (osc.m_startTime < 2.6 ||
                //osc.m_startTime > 3.0)
            //{
                //continue;
            //}

            std::vector<double> times;
            std::vector<double> radii;
            std::vector<double> freqs;
            std::vector<double> transfer;
            std::vector<double> pressure;

            // Keep track of the last forcing event time, to filter them
            double prevStartTime = -std::numeric_limits<double>::infinity();

            bool goodRadiiSeq = true;

            do
            {
                //std::cout << curBubNum << std::endl;

                used.insert(curBubNum);
                //if (times.empty())
                {
                    if (!curBub->m_times.empty() && curBub->m_startTime < curBub->m_times.at(0))
                    {
                        times.push_back(curBub->m_startTime);
                        radii.push_back(curBub->m_radius);
                        freqs.push_back(curBub->m_frequencies.at(0));
                        transfer.push_back(curBub->m_transfer.at(0));
                        pressure.push_back(curBub->m_pressure.at(0));
                    }
                }

                times.insert(times.end(),
                             curBub->m_times.begin(),
                             curBub->m_times.end());

                radii.insert(radii.end(),
                             curBub->m_times.size(),
                             curBub->m_radius);

                freqs.insert(freqs.end(),
                             curBub->m_frequencies.begin(),
                             curBub->m_frequencies.end());

                transfer.insert(transfer.end(),
                                curBub->m_transfer.begin(),
                                curBub->m_transfer.end());

                pressure.insert(pressure.end(),
                                curBub->m_pressure.begin(),
                                curBub->m_pressure.end());


                for (int i = 0; i < curBub->m_times.size(); ++i)
                {
                    osc.m_trackedBubbleNumbers[curBub->m_times[i]] = curBub->m_nums.at(i);
                }

                // Check to filter out bad low freq bubble in pouringGlass17
                bool forceThisBub = true;
                //if (curBub->m_frequencies.size() >= 1 && curBub->m_startTime >= .9 && curBub->m_startTime < 1.1 &&
                    //curBub->m_frequencies.at(0) <= 500 * 2 * M_PI)
                if (false && curBub->m_frequencies.size() > 0 && curBub->m_frequencies.at(0) <= 500 * 2 * M_PI)
                {
                    forceThisBub = false;
                    //std::cout << "Filtered bub at time: " << curBub->m_startTime << std::endl;
                    ++numFiltered;
                }

                // Filter forcing events, only use this one if it is as least 1ms
                // since the last one
                //if (curBub->m_startTime - prevStartTime > 0.001)
                //if (curBub->m_frequencies.size() > 0 && forceThisBub)
                if (forceThisBub)
                {
                    std::shared_ptr<ForcingFunction> force = makeForcingFunc(*curBub,
                                                                             curBubNum,
                                                                             singleBubbles);

                    osc.m_forcing.insert(std::make_pair(curBub->m_startTime, force));
                    prevStartTime = curBub->m_startTime;
                }

                osc.m_bubbleNumbers.push_back(curBubNum);

                bool lastBub = true;

                if (curBub->m_endType == Bubble::MERGE)
                {
                    // If this is the largest parent bubble, continue,
                    // else end this bubble

                    if (curBub->m_nextBubbles.size() != 1)
                    {
                        throw std::runtime_error(std::string("did not merge to one bubble ") + std::to_string(curBubNum));
                    }

                    const Bubble& nextBub = singleBubbles.at(curBub->m_nextBubbles.at(0));

                    int largestParent = largestBubble(nextBub.m_prevBubbles,
                                                      singleBubbles);

                    if (largestParent == curBubNum)
                    {
#ifdef CHAIN_BUBBLES
                        lastBub = false;
                        curBubNum = curBub->m_nextBubbles.at(0);
                        curBub = &singleBubbles.at(curBubNum);
#endif
                    }
                }
                else if (curBub->m_endType == Bubble::SPLIT)
                {
                    // Continue this bubble to the largest child bubble where this is the largest parent

                    if (curBub->m_nextBubbles.size() < 2)
                    {
                        std::cerr << "split to less than 2 bubbles" << " " << curBub->m_startTime << " " << curBub->m_endTime << " " << curBubNum << std::endl;
                        throw std::runtime_error("split to less than 2 bubbles");
                    }

                    // Sort children in order of size
                    std::multimap<double, int, std::greater<double>> children;
                    for (int childNum : curBub->m_nextBubbles)
                    {
                        children.insert(std::make_pair(singleBubbles.at(childNum).m_radius,
                                                       childNum));
                    }

                    bool found = false;

                    for (auto &child : children)
                    {
                        if (!used.count(child.second))
                        {
                            int largestParent = largestBubble(singleBubbles.at(child.second).m_prevBubbles,
                                                              singleBubbles);

                            if (largestParent == curBubNum)
                            {
#ifdef CHAIN_BUBBLES
                                found = true;
                                curBubNum = child.second;
                                curBub = &singleBubbles.at(curBubNum);
                                lastBub = false;
#endif
                                break;
                            }
                        }
                    }
                }

                if (lastBub)
                {
                    osc.m_endType = curBub->m_endType;

                    if (!curBub->m_times.empty() && curBub->m_endTime > times.back())
                    {
                        //times.push_back(curBub->m_endTime);
                        //radii.push_back(curBub->m_radius);
                        //freqs.push_back(curBub->m_frequencies.back());
                        //transfer.push_back(curBub->m_transfer.back());
                        //pressure.push_back(curBub->m_pressure.back());
                    }

                    // Extend with an exponential
                    double extendTime = 0;
                    // TODO: figure out what to do here with the wavesolver
                    if (0 && times.size() >= 3 && curBub->m_endType == Bubble::COLLAPSE)
                    {
                        extendTime = 0.030; // TODO: set this based on bubble size?
                        double factor = 1;

                        double st = 0.001 + times.back() - times.front();
                        Eigen::VectorXd extraTimes = Eigen::VectorXd::LinSpaced(31, st, st + extendTime);

                        int ind = times.size() - 1;

                        if (0)
                        {
                            for (int i = 0; i < extraTimes.size(); ++i)
                            {
                                times.push_back(extraTimes(i) + times.at(0));
                                freqs.push_back(freqs.back());
                                transfer.push_back(transfer.back());
                                radii.push_back(radii.back());
                                pressure.push_back(pressure.back());
                            }
                        }
                        else
                        {
                            double y2 = freqs.at(ind);
                            double y1 = freqs.at(ind-1);
                            double t2 = times.at(ind) - times.at(0);
                            double t1 = times.at(ind-1) - times.at(0);
                            double m = (y2 - y1) / (t2 - t1) * factor;

                            if (m < 1000*1000)
                            {
                                double freqDiff = y2 - y1;

                                double bf = m / y2;
                                double cf = y2 * std::exp(-m * t2 / y2);

                                y2 = transfer.at(ind);
                                y1 = transfer.at(ind-1);
                                m = (y2 - y1) / (t2 - t1) * factor;

                                double xfrDiff = y2 - y1;

                                double bt = m / y2;
                                double ct = y2 * std::exp(-m * t2 / y2);

                                //if (freqDiff > 0 && freqDiff < 100 * 2 * M_PI
                                    //&& fabs(radii.at(ind) - radii.at(ind-1)) < 1e-10)
                                if (fabs(radii.at(ind) - radii.at(ind-1)) < 1e-12 && freqs.at(ind) > freqs.at(ind-1))
                                {
                                    for (int i = 0; i < extraTimes.size(); ++i)
                                    {
                                        double freq = cf * std::exp(bf * extraTimes(i));
                                        double trans = ct * std::exp(bt * extraTimes(i));

                                        //if (freq > 20000 * 2 * M_PI || freq >
                                            //freqs.back() * 1.2 ||
                                            //fabs(trans) > transfer.back() * 1.1)
                                        //{
                                            //break;
                                        //}

                                        times.push_back(extraTimes(i) + times.at(0));
                                        freqs.push_back(freq);
                                        //transfer.push_back(trans);
                                        transfer.push_back(transfer.back());
                                        radii.push_back(radii.back());
                                        pressure.push_back(pressure.back());
                                    }
                                }
                            }
                        }
                    }

                    // Create the pressure forcing function
                    if (times.size() >= 1)
                    {
                        const double precTime = 0.001;
                        times.insert(times.begin(), times.front() - precTime);
                        //osc.m_startTime -= precTime;

                        if (oscStartType == Bubble::ENTRAIN)// && firstLength >= 0.0005)
                        {
                            //pressure.insert(pressure.begin(), pressure.front());
                            pressure.insert(pressure.begin(), 101325.0);
                        }
                        else
                        {
                            pressure.insert(pressure.begin(), pressure.front());
                        }

                        std::shared_ptr<ForcingFunction> pForce(new PressureForcing(times, pressure, 0.0003));
                        //std::shared_ptr<ForcingFunction> pForce(new PressureForcing(times, pressure, 0.25 * freqs.front() / 2 / M_PI));
                        //osc.m_forcing.insert(std::make_pair(osc.m_startTime, pForce));

                        freqs.insert(freqs.begin(), freqs.front());
                        transfer.insert(transfer.begin(), transfer.front());
                        radii.insert(radii.begin(), radii.front());

                        // Add extra times at the end
                        while (times.size() < 4)
                        {
                            times.push_back(times.back() + precTime);
                            pressure.push_back(pressure.back());
                            freqs.push_back(freqs.back());
                            transfer.push_back(transfer.back());
                            radii.push_back(radii.back());
                        }
                    }

                    osc.m_frequencies.setData (Eigen::Map<Eigen::VectorXd> (times.data(), times.size()),
                                               Eigen::Map<Eigen::VectorXd> (freqs.data(), freqs.size()));

                    osc.m_transfer.setData (Eigen::Map<Eigen::VectorXd> (times.data(), times.size()),
                                            Eigen::Map<Eigen::VectorXd> (transfer.data(), transfer.size()));

                    osc.m_radii.setData (Eigen::Map<Eigen::VectorXd> (times.data(), times.size()),
                                         Eigen::Map<Eigen::VectorXd> (radii.data(), radii.size()));

                    osc.m_pressure.setData (Eigen::Map<Eigen::VectorXd> (times.data(), times.size()),
                                            Eigen::Map<Eigen::VectorXd> (pressure.data(), pressure.size()));

                    // Debugging
                    //double maxR = 1, minR = 1;
                    //if (osc.m_radii.getData().rows() > 0)
                    //{
                        //maxR = osc.m_radii.getData().col(1).maxCoeff();
                        //minR = osc.m_radii.getData().col(1).minCoeff();
                    //}

                    //if ((maxR - minR) / maxR > 0.4)
                    //{
                        //std::cerr << "Bad radii sequence: " << _oscillators.size() << " " << osc.m_radii.getData().col(1).maxCoeff()
                                  //<< " " << osc.m_radii.getData().col(1).minCoeff() << std::endl;

                        //goodRadiiSeq = false;
                    //}

                    // End this bubble
                    osc.m_endTime = curBub->m_endTime + extendTime;
                    curBub = NULL;

                    // TODO: tmp debugging
                    if (freqs.size() == 0 || freqs.at(0) > 700 * 2 * M_PI) goodRadiiSeq = false;

                    // Check to filter out bad low freq bubble in pouringGlass17
                    //if (freqs.size() >= 1 && osc.m_startTime >= .97 && osc.m_startTime < 1.05 &&
                        //freqs.at(0) <= 500 * 2 * M_PI)
                    //{
                        //goodRadiiSeq = false;
                        //++numFiltered;
                    //}
                }
            } while (curBub);

            if (goodRadiiSeq)
            {
                _oscillators.push_back(osc);
            }
        }
    }

    std::cout << std::endl;
    std::cout << "numFiltered: " << numFiltered << std::endl;
}


//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
computeVelocities(REAL time)
{
    using namespace std;
    using namespace Eigen;

    static bool first = true;
    static vector<ofstream> oscOutputs;

    if (first)
    {
        for (int i = 0; i < _oscillators.size(); ++i)
        {
            ostringstream os;
            os << "velValues_" << i << ".txt";
            oscOutputs.emplace_back(os.str().c_str());
        }

        first = false;
    }

    bool useT1 = std::fabs(time - _t1) < std::fabs(_t2 - time);

    // Interpolate onto the mesh that is closest
    _m = useT1 ? &_m1 : &_m2;
    //_kd = useT1 ? _fullKd1 : _fullKd2;
    _kd = useT1 ? _kd1 : _kd2;

    double localT1 = time - _dt;
    double localT2 = time + _dt;

    _velT1 = VectorXd::Zero(_m->m_surfTris.size());
    _velT2 = VectorXd::Zero(_m->m_surfTris.size());

    std::vector<int> closest1, closest2;

	// Loop through the oscillators
    for (int i = 0; i < _oscillators.size(); ++i)
    {
        Oscillator &osc = _oscillators[i];

        if (!osc.isActive(time) || osc.m_trackedBubbleNumbers.empty()) continue;

        bool existT1 = osc.m_startTime <= _t1 && osc.m_endTime > _t1;
        bool existT2 = osc.m_startTime <= _t2 && osc.m_endTime > _t2;

        if (!existT1 && !existT2) continue;

        std::shared_ptr<PointKDTree> tree1(existT1 ? _kd1 : _kd2);
        std::shared_ptr<PointKDTree> tree2(existT2 ? _kd2 : _kd1);

        SurfaceVelocityData &vel1 = existT1 ? _v1 : _v2;
        SurfaceVelocityData &vel2 = existT2 ? _v2 : _v1;

        Mesh &localM1 = existT1 ? _m1 : _m2;
        Mesh &localM2 = existT2 ? _m2 : _m1;

        double lookupT = existT1 ? _t1 : _t2;
        auto iter = osc.m_trackedBubbleNumbers.lower_bound(lookupT - 1e-15);
        if (iter == osc.m_trackedBubbleNumbers.end())
        {
            cout << "t1: " << _t1 << ", t2: " << _t2 << endl;
            for (auto b : osc.m_trackedBubbleNumbers)
            {
                cout << b.first << " " << b.second << endl;
            }

            throw runtime_error("bad tracked bubble numbers lookup t1");
        }

        int bubbleNumber1 = iter->second;

        lookupT = existT2 ? _t2 : _t1;
        iter = osc.m_trackedBubbleNumbers.lower_bound(lookupT - 1e-15);
        if (existT2 && iter == osc.m_trackedBubbleNumbers.end())
        {
            cout << "t1: " << _t1 << ", t2: " << _t2 << endl;
            for (auto b : osc.m_trackedBubbleNumbers)
            {
                cout << b.first << " " << b.second << endl;
            }

            cout << "end time: " << osc.m_endTime << endl;
            cout << "times: " << osc.m_frequencies.getData().col(0) << endl;

            throw runtime_error("bad tracked bubble numbers lookup t2");
        }

        int bubbleNumber2 = iter->second;

        oscOutputs[i] << time << " " << bubbleNumber1 << " " << bubbleNumber2 << endl;

        //if (!existT1)
        //{
        //    bubbleNumber1 = bubbleNumber2;
        //}
        //else if (!existT2)
        //{
        //    bubbleNumber2 = bubbleNumber1;
        //}

        // DEBUGGING
        bool useFirstOnly = true;
        struct DebugData
        {
            std::shared_ptr<PointKDTree> tree;
            SurfaceVelocityData vel;
            int bubNum;
            Mesh localM;
        };

        if (useFirstOnly)
        {
            // Use the first velocity without any interpolation
            if (!osc.debugData)
            {
                DebugData * ddata = new DebugData();
                ddata->tree = tree1;
                ddata->vel = vel1;
                ddata->bubNum = bubbleNumber1;
                ddata->localM = localM1;

                osc.debugData = static_cast<void*>(ddata);
            }
        }
        // DEBUGGING

        for (int j = 0; j < _m->m_surfTris.size(); ++j)
        {
            if (_m->m_triType.at(_m->m_surfTris.at(j)) != Mesh::FLUID_AIR)
            {
                continue;
            }

            // Now interpolate to correct times
            MLSVal val1, val2;
            MLSPoint p = _m->m_surfTriCenters[j];

            tree1->find_nearest(p,
                                5,
                                closest1);

            tree2->find_nearest(p,
                                5,
                                closest2);

            // Values at t1 and t2
            try
            {
                // TODO: vel1 has only fluid surface data, surfTriCenters also has rigid triangles
                val1 = _mls.lookup(p,
                                   localM1.m_surfTriCenters,
                                   vel1.at(bubbleNumber1),
                                   -1,
                                   NULL,
                                   &closest1);

                val2 = _mls.lookup(p,
                                   localM2.m_surfTriCenters,
                                   vel2.at(bubbleNumber2),
                                   -1,
                                   NULL,
                                   &closest2);

                if (std::fabs(val1(0)) > 1e6)
				{
					std::cout << "bad val1: " << val1 << std::endl;
					exit(1);
				}

                if (std::fabs(val2(0)) > 1e6)
				{
					std::cout << "bad val2: " << val2 << std::endl;
					exit(1);
				}
            }
            catch (...)
            {
                std::cout << "time: " << time << ", i: " << i << ", j: " << j << std::endl;
                std::cout << "num1: " << bubbleNumber1 << ", num2: " << bubbleNumber2 << std::endl;

                for (auto iter = vel1.begin(); iter != vel1.end(); ++iter)
                {
                    std::cout << iter->first << std::endl;
                }

                exit(1);
            }

            if (useFirstOnly)
            {
                DebugData *ddata = static_cast<DebugData*>(osc.debugData);


                ddata->tree->find_nearest(p,
                                          5,
                                          closest1);

                val1 = _mls.lookup(p,
                                   ddata->localM.m_surfTriCenters,
                                   ddata->vel.at(ddata->bubNum),
                                   -1,
                                   NULL,
                                   &closest1);

                val2 = val1;
            }

            // DEBUGGING, set surface velocity to constant 1 for now
            val1.setConstant(1);
            val2.setConstant(1);

            // Now interpolate
            double pct = (localT1 - _t1) / (_t2 - _t1);
            //if (pct < 0 || pct > 1)
            //{
            //    std::cout << "localT1: " << localT1 << ", t1: " << _t1 << ", t2: " << _t2 << std::endl;
            //    throw std::runtime_error("bad pct");
            //}
            if (pct < 0)
            {
                _velT1(j) += osc.m_lastVals.col(0)(0) * val1(0);
            }
            else if (pct > 1)
            {
                _velT1(j) += osc.m_lastVals.col(0)(0) * val2(0);
            }
            else
            {
                _velT1(j) += osc.m_lastVals.col(0)(0) * ( pct * val2(0) + (1 - pct) * val1(0) );
            }

            pct = (localT2 - _t1) / (_t2 - _t1);
            //if (pct < 0 || pct > 1)
            //{
            //    std::cout << "localT2: " << localT2 << ", t1: " << _t1 << ", t2: " << _t2 << std::endl;
            //    throw std::runtime_error("bad pct");
            //}
            if (pct < 0)
            {
                _velT2(j) += osc.m_lastVals.col(2)(0) * val1(0);
            }
            else if (pct > 1)
            {
                _velT2(j) += osc.m_lastVals.col(2)(0) * val2(0);
            }
            else
            {
                _velT2(j) += osc.m_lastVals.col(2)(0) * ( pct * val2(0) + (1 - pct) * val1(0) );
            }
        }
	}
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
projectToSurface()
{
    _projectedAccel.resize(_surfaceMesh->num_triangles());
    _projectedAccel.setZero();

    // For each triangle, find closest one on parent mesh and project
    //for (int i = 0; i < _accel.size(); ++i)
    for (int i = 0; i < _surfaceMesh->num_triangles(); ++i)
    {
        Point3<REAL> p = _surfaceMesh->triangle_centroid(i);

        Eigen::Vector3d ep(p.x, p.y, p.z);

        int nearestTri = _kd->find_nearest(ep);

        if (_m->m_triType.at(_m->m_surfTris.at(nearestTri)) == Mesh::FLUID_AIR)
        {
            int index = nearestTri;
            if (std::fabs(_accel(index)) > std::fabs(_projectedAccel(i)))
            {
                _projectedAccel(i) = _accel(index);
            }

            //Vector3<double> boxN = _surfaceMesh->triangle_normal(i).normalized();
            //Eigen::Vector3d surfN = _m->triangleNormal(_m->m_surfTris.at(nearestTri));

            //Eigen::Vector3d boxNe;
            //boxNe << boxN.x, boxN.y, boxN.z;

            //if (surfN.dot(boxNe) < 0.7)
            //{
            //    std::cout << "surfN: " << surfN.transpose() << std::endl;
            //    std::cout << "boxN: " << boxNe.transpose() << std::endl;
            //    std::cout << "nearestTri: " << nearestTri << std::endl;

            //    Tuple3ui tri = _surfaceMesh->triangles().at(i);
            //    Point3<double> v1= _surfaceMesh->vertices().at(tri.x);
            //    std::cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << std::endl;
            //    Point3<double> v2= _surfaceMesh->vertices().at(tri.y);
            //    std::cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << std::endl;
            //    Point3<double> v3= _surfaceMesh->vertices().at(tri.z);
            //    std::cout << "v3: " << v3.x << " " << v3.y << " " << v3.z << std::endl;

            //    Vector3<double> n = ((v2 - v1).crossProduct(v3-v1)).normalized();
            //    std::cout << "n: " << n.x << " " << n.y << " " << n.z << std::endl;
            //}
        }
    }

    // DEBUGGING, write to file
    //using namespace std;
    //ofstream of("surfAccel.ply");

    //of << "ply" << endl;
    //of << "format ascii 1.0" << endl;
    //of << "element vertex " << _surfaceMesh->num_vertices() << endl;
    //of << "property float x" << endl;
    //of << "property float y" << endl;
    //of << "property float z" << endl;
    //of << "element face " << _surfaceMesh->num_triangles() << endl;
    //of << "property list uchar uint vertex_index" << endl;
    //of << "property float accel" << endl;
    //of << "end_header" << endl;

    //for (int i = 0; i < _surfaceMesh->num_vertices(); ++i)
    //{
    //    of << _surfaceMesh->vertex(i).x << " "
    //       << _surfaceMesh->vertex(i).y << " "
    //       << _surfaceMesh->vertex(i).z << endl;
    //}

    //for (int i = 0; i < _surfaceMesh->num_triangles(); ++i)
    //{
    //    of << "3 "
    //       << _surfaceMesh->triangles()[i].x << " "
    //       << _surfaceMesh->triangles()[i].y << " "
    //       << _surfaceMesh->triangles()[i].z << " "
    //       << _projectedAccel(i) << endl;
    //}
}

