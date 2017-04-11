/**
 * Test the ODE integrators on a damped harmonic oscillator.
 */

#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <cmath>
#include <sndfile.h>
#include <Eigen/Dense>
#include <MathUtils/ODEInt.h>
#include <bubbleForcing.hpp>
#include <timeSeries.hpp>
#include <constants.h>
#include <bubblePop.hpp>

#include <vtkXMLUnstructuredGridReader.h>

using namespace Eigen;
using namespace MathUtils;

static double sampleFreq;
static double dt;
static double timeLength;

// Use row vectors for sound buffers
// then data will be in correct order for sndfile writes
// Matters with more than one channel
static RowVectorXd buffer;

static BubblePop pop;

static std::default_random_engine s_rnd;
static std::uniform_real_distribution<double> s_popAmp(.001, .03); // Relative to bubble sound

// Class to represent bubble events
class Bubble
{
public:
    enum EventType
    {
        ENTRAIN,
        MERGE,
        SPLIT,
        COLLAPSE
    };

    EventType m_startType;
    EventType m_endType;

    // Bubbles that this one came from
    std::vector<int> m_prevBubbles;

    // Bubbles that this one merged/split into
    std::vector<int> m_nextBubbles;

    double m_startTime;
    double m_endTime;

    std::vector<double> m_times;
    std::vector<double> m_frequencies;
    std::vector<double> m_transfer;
    std::vector<double> m_pressure;
    std::vector<int> m_nums;

    // Equivalent radius
    double m_radius;
};

// Class used to represent an oscillator
// can be multiple Bubble objects chained together
class Oscillator
{
public:
    double m_startTime;
    double m_endTime;

    Bubble::EventType m_endType;

    TimeSeries<double> m_frequencies;
    TimeSeries<double> m_transfer;
    TimeSeries<double> m_radii;
    TimeSeries<double> m_pressure;

    std::multimap<double, std::shared_ptr<ForcingFunction>> m_forcing;
    std::vector<int> m_bubbleNumbers;
    std::map<double, int> m_trackedBubbleNumbers;

    double m_currentTime;
    Eigen::Vector3d m_lastVals;
    Eigen::Vector2d m_state;

    Oscillator()
        : m_currentTime(-1)
    {
        m_lastVals.setZero();
        m_state.setZero();
    }

    bool
    isActive(REAL time)
	{
		if (m_startTime > time) return false;

        double endTime = m_endTime > 1000 ? m_startTime + .001 : m_endTime;

        if (m_currentTime > 0 && m_currentTime > endTime && m_state.norm() <= 1e-15)
        {
            return false;
        }

        return true;
	}
};

static int
largestBubble(const std::vector<int> &bubbles,
              const std::map<int, Bubble> &allBubbles)
{
    double r = 0;
    int maxB = 0;

    for (int i : bubbles)
    {
        if (allBubbles.at(i).m_radius > r)
        {
            maxB = i;
            r = allBubbles.at(i).m_radius;
        }
    }

    return maxB;
}

/**
 * Bubble harmonic oscillator
 * q'' + 2 * beta * q' + omega * q = F/m
 * with m = 4 pi r^3 / rho
 */
template<typename State>
class BubbleOscillator
{
public:
    typedef typename State::Scalar T;
    typedef BubbleOscillator<State> Evaluator;
    typedef TimeSeries<T> TsType;

    double curBeta;

    static double
    undampedNaturalFreq (T radius)
    {
        using namespace std;

        return sqrt(3 * GAMMA * ATM - 2 * SIGMA/radius) / (radius * sqrt(RHO_WATER));
    }

    static double
    calcBeta (T radius,
              T w0)
    {
        using namespace std;

        T dr = w0 * radius / CF;
        T dvis = 4 * MU / (RHO_WATER * w0 * radius*radius);
        T phi = 16. * GTH * G / (9 * (GAMMA-1)*(GAMMA-1) * w0);
        T dth = 2*(sqrt(phi - 3) - (3 * GAMMA - 1) / (3 * (GAMMA - 1))) / (phi - 4);

        T dtotal = dr + dvis + dth;

        return w0 * dtotal / sqrt(dtotal * dtotal + 4);
        //return 0.5 * w0 * dtotal;
    }

    BubbleOscillator (const Oscillator &osc,
                      int bubNum,
                      bool write = false)
        : m_osc(osc),
          //m_r(radius),
          //m_m(RHO_WATER / (4. * M_PI * radius)),
          //m_forceFunc(forcing),
          //m_frequencies(frequencies),
          writeFrequencies(write)
    {
        if (writeFrequencies)
        {
            of.open("frequencies-" + std::to_string(bubNum) + ".txt");
        }
    }

    State
    operator() (const State &currState,
                T currTime)
    {
        // currState = [x'; x]

        double startTime = m_osc.m_frequencies.getData()(0,0);

        if (currTime < startTime)
        {
            return currState;
        }

        // The forcing here is in units of pressure
        T force = 0;

        auto forceIter = m_osc.m_forcing.begin();
        while (forceIter != m_osc.m_forcing.end() && forceIter->first <= currTime)
        {
            bool allForcing = false;
            if (allForcing)
            {
                force += forceIter->second->value(currTime - forceIter->first);
            }
            else
            {
                auto nextIter = forceIter;
                ++nextIter;
                if (nextIter == m_osc.m_forcing.end() ||
                    nextIter->first > currTime)
                {
                    force += forceIter->second->value(currTime - forceIter->first);
                    break;
                }
            }
            ++forceIter;
        }

        T w0 = m_osc.m_frequencies.interp(currTime);

        double r = m_osc.m_radii.interp(currTime);
        T beta = calcBeta(r,
                          w0);

        if (writeFrequencies)
        {
            of << std::setprecision(12) << currTime << " " << w0 / 2 / M_PI << " "
               << force << " " << r << " " << beta << std::endl;
            of.flush();
        }

        curBeta = beta;

        // Working in the pressure-volume frame
        //double mvp = RHO_WATER / (4 * M_PI * r);

        // Use the actual calculated mass, not the leighton spherical one
        double k = GAMMA * m_osc.m_pressure.interp(currTime) / (4./3. * M_PI * r*r*r);
        double mvp = k / w0 / w0;

        T acc = force / (mvp) - 2 * beta * currState(0) - w0 * w0 * currState(1);

        State retState;
        retState.resize(2, 1);
        retState << acc, currState(0);

        return retState;
    }

private:
    const Oscillator &m_osc;

    bool writeFrequencies;

    std::ofstream of;
};

Bubble::EventType
eventType(char type)
{
    switch (type)
    {
        case 'N':
            return Bubble::ENTRAIN;
            break;

        case 'M':
            return Bubble::MERGE;
            break;

        case 'S':
            return Bubble::SPLIT;
            break;

        case 'C':
            return Bubble::COLLAPSE;
            break;

        default:
            throw std::runtime_error("Invalid event type");
    }
}

enum FreqType
{
    CAPACITANCE,
    MINNAERT
};

// Returns the (angular) Minnaert frequency
static double
minnaertFreq(double radius)
{
    return std::sqrt(3 * GAMMA * ATM - 2 * SIGMA / radius) / (radius * std::sqrt(RHO_WATER));
}

Bubble
parseDataFile (const char *fileName,
               FreqType fType)
{
    Bubble bub;

    std::ifstream in(fileName);

    std::vector<double> times, freqs, transfs;

    double t, freq;
    std::complex<double> trans;

    // First line is the start info
    std::string line;
    std::getline(in, line);

    char startType = line.at(7);
    std::istringstream is(line.substr(8));
    is >> bub.m_startTime;
    int bubNum;
    while (is >> bubNum)
    {
        bub.m_prevBubbles.push_back(bubNum);
    }

    bub.m_startType = eventType(startType);

    // Next line is radius info
    std::getline(in, line);
    line = line.substr(6);
    bub.m_radius = std::stod(line);

    while ( in.peek() != 'E')
    {
        std::getline(in, line);
        std::istringstream newis(line);
        newis >> t >> freq >> trans;

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

    bub.m_endType = eventType(endType);

    // TODO: set a Minnaert freq for the bubble if no solve data present?
    // How to handle transfer?
    if (!times.empty())
    {
        // Update freqs for start and end
        if (bub.m_startTime < times.at(0))
        {
            times.insert(times.begin(), bub.m_startTime);
            freqs.insert(freqs.begin(), freqs.at(0));
            transfs.insert(transfs.begin(), transfs.at(0));
        }

        if (bub.m_endTime > times.back())
        {
            times.push_back(bub.m_endTime);
            freqs.push_back(freqs.back());
            transfs.push_back(transfs.back());
        }

        bub.m_times = times;
        bub.m_frequencies = freqs;
        bub.m_transfer = transfs;
    }

    return bub;

    //std::cout << "Radius: " << radius << std::endl
              //<< "freqs: \n" << frequencies.getData() << std::endl
              //<< "transfer: \n" << transfer.getData() << std::endl;
}

std::pair<int, Bubble>
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

    bub.m_startType = eventType(startType);

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

    bub.m_endType = eventType(endType);

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

int
bubbleNumber(char *fileName)
{
    std::string fName(fileName);
    int pos = fName.find_first_of("-") + 1;

    return std::stoi(fName.substr(pos));
}

std::map<int, Bubble>
parseConfigFile(char* infoFile, FreqType fType)
{
    std::map<int, Bubble> bubbles;

    std::ifstream in(infoFile);

    while (in.good())
    {
        std::pair<int, Bubble> bub = parseBubbleInfo(in, fType);

        if (bub.first >= 0)
        {
            bubbles.insert(bub);
        }

        std::cout << "\rLoaded " << bub.first << "    ";
        std::cout.flush();
    }

    std::cout << std::endl;

    return bubbles;
}

static std::shared_ptr<ForcingFunction>
makeForcingFunc(const Bubble& bub,
                int curBubNum,
                const std::map<int, Bubble>& bubbles)
{
    bool noSplit = false;
    bool noMerge = false;
    bool noEntrain = false;

    std::shared_ptr<ForcingFunction> forcing;

    //if (bub.m_endTime - bub.m_startTime < 0.01)
    //{
        //forcing.reset(new ZeroForcing());
    //}

    // If this bubble came from another bubble, set the state correctly
    if (bub.m_startType == Bubble::SPLIT)
    {
        if (bub.m_prevBubbles.size() > 1)
        {
            throw std::runtime_error("split from more than one bubble");
        }

        int prevBub = bub.m_prevBubbles.at(0);
        double minR = std::numeric_limits<double>::infinity();
        int minBub = -1;

        if (bubbles.count(prevBub))
        {
            for ( auto b : bubbles.at(prevBub).m_nextBubbles )
            {
                minR = std::min(minR, bubbles.at(b).m_radius);
                minBub = b;
            }
        }

        //if (prevBub == 270)
        //{
            //std::cout << "r: " << bub.m_radius << " minR: " << minR << " start: " << bub.m_startTime
                      //<< " end: " << bub.m_endTime << std::endl;
        //}

        //std::cout << "minR: " << minR << std::endl;

        //forcing.reset(new CzerskiJetForcing(bub.m_radius));

        if (bubbles.at(prevBub).m_radius < bub.m_radius)
        {
            forcing.reset(new ZeroForcing());
        }
        else
        {
            if (curBubNum == minBub)
            {
                forcing.reset(new CzerskiJetForcing(bub.m_radius,
                                                    5000,
                                                    //bub.m_endTime - bub.m_startTime,
                                                    //std::min(bub.m_endTime - bub.m_startTime, 0.5 / (3.0 / minR)),
                                                    //ETA_SPLIT,
                                                    ETA,
                                                    false,
                                                    false));
                                                    //true));
            }
            else
            {
                forcing.reset(new CzerskiJetForcing(bub.m_radius,
                                                    5000,
                                                    //bub.m_endTime - bub.m_startTime,
                                                    //std::min(bub.m_endTime - bub.m_startTime, 0.5 / (3.0 / minR)),
                                                    //ETA_SPLIT,
                                                    ETA,
                                                    false,
                                                    false));
                                                    //true));
            }
        }
        //forcing.reset(new CzerskiJetForcing(minR,
                                            //std::min(bub.m_endTime - bub.m_startTime, 0.5 / (3.0 / minR))));
        if (noSplit)
        {
            forcing.reset(new ZeroForcing());
        }
    }
    else if (bub.m_startType == Bubble::MERGE)
    {
        if (bub.m_prevBubbles.size() < 2)
        {
            forcing.reset(new ZeroForcing());
            return forcing;
            //throw std::runtime_error("merged from less than two bubbles");
        }

        if (bub.m_prevBubbles.size() > 2)
        {
            forcing.reset(new ZeroForcing());
        }
        else
        {
            // Make sure all parent bubbles merged completely with this one
            // Sometimes a bubble can split into two, and the smaller fragment
            // can merge with another one. This is probably just resolution and/or
            // bubble tracking issues
            bool allMerge = true;
            for (auto p : bub.m_prevBubbles)
            {
                allMerge = allMerge && bubbles.at(p).m_endType == Bubble::MERGE;
            }

            if (allMerge)
            {
                int p1 = bub.m_prevBubbles.at(0);
                int p2 = bub.m_prevBubbles.at(1);
                double r1 = 0, r2 = 0;

                if (bubbles.count(p1) && bubbles.count(p2))
                {
                    r1 = bubbles.at(p1).m_radius;
                    r2 = bubbles.at(p2).m_radius;
                }
                else
                {
                    throw std::runtime_error("Missing parent");
                }

                if (r1 + r2 > bub.m_radius)
                {
                    double v1 = 4./3. * M_PI * r1 * r1 * r1;
                    double v2 = 4./3. * M_PI * r2 * r2 * r2;
                    double vn = 4./3. * M_PI * bub.m_radius * bub.m_radius * bub.m_radius;

                    double diff = v1 + v2 - vn;

                    if (diff > std::max<double>(v1, v2))
                    {
                        // In this case both parent bubbles must have split off...
                        forcing.reset(new ZeroForcing());
                    }
                    else
                    {
                        if (v1 > v2)
                        {
                            v1 -= diff;
                        }
                        else
                        {
                            v2 -= diff;
                        }

                        r1 = std::pow(3./4. / M_PI * v1, 1./3.);
                        r2 = std::pow(3./4. / M_PI * v2, 1./3.);

                        forcing.reset(new MergeForcing(bub.m_radius,
                                                       r1,
                                                       r2,
                                                       5000));
                                                       //bub.m_endTime - bub.m_startTime));
                    }
                }
                else
                {
                    forcing.reset(new MergeForcing(bub.m_radius,
                                                   r1,
                                                   r2,
                                                   5000));
                                                   //bub.m_endTime - bub.m_startTime));
                }
            }
            else
            {
                forcing.reset(new ZeroForcing());
            }

            if (noMerge)
            {
                forcing.reset(new ZeroForcing());
            }
        }
    }
    else if (bub.m_startType == Bubble::ENTRAIN)
    {
        forcing.reset(new CzerskiJetForcing(bub.m_radius,
                                            5000,
                                            //bub.m_endTime - bub.m_startTime,
                                            //ETA_ENTRAIN,
                                            ETA,
                                            true,
                                            false,
                                            1));

        if (noEntrain)
        {
            forcing.reset(new ZeroForcing());
        }
    }
    else
    {
        throw std::runtime_error("bad start event");
    }

    return forcing;
}

// Attach the bubbles together to get the full oscillator lifetime
std::vector<Oscillator>
makeOscillators(const std::map<int, Bubble> &singleBubbles)
{
#define CHAIN_BUBBLES
    std::vector<Oscillator> oscillators;
    std::set<int> used;

    int numFiltered = 0;

    for (const auto& bubPair : singleBubbles)
    {
        std::cout << "\rTracking " << bubPair.first << " of " << singleBubbles.size() << "       ";
        std::cout.flush();

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
                    std::cout << "Filtered bub at time: " << curBub->m_startTime << std::endl;
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
                    if (1 && times.size() >= 3 && curBub->m_endType == Bubble::COLLAPSE)
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
                        //std::cerr << "Bad radii sequence: " << oscillators.size() << " " << osc.m_radii.getData().col(1).maxCoeff()
                                  //<< " " << osc.m_radii.getData().col(1).minCoeff() << std::endl;

                        //goodRadiiSeq = false;
                    //}

                    // End this bubble
                    osc.m_endTime = curBub->m_endTime + extendTime;
                    curBub = NULL;

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
                oscillators.push_back(osc);
            }
        }
    }

    std::cout << std::endl;
    std::cout << "numFiltered: " << numFiltered << std::endl;

    return oscillators;
}

void
runOscillator(int oscNum,
              Oscillator& osc)
{
    typedef BubbleOscillator<Vector2d> DS;

    // Really checking if it is less than 2, but we add
    // an extra for the start and end ones
    if (osc.m_frequencies.getData().rows() < 1
        )

        //bub.m_radius < 0.001)
        //bub.m_startType == Bubble::MERGE)
        //bub.m_startType == Bubble::SPLIT)
    //if (false)
    {
        return;
    }

    //MidpointMethod<Vector2d> integrator;
    RK4<Vector2d> integrator;

    Vector2d state(0, 0);

    DS bubOsc(osc,
              oscNum);

    //std::cout << "Radius: " << bub.m_radius << std::endl;
    //std::cout << "Freqs: " << bub.m_frequencies.getData() << std::endl;

    // TODO: check indices here, are we double counting the end state
    // of one bubble and the start state of another?
    double t = osc.m_startTime;
    int counter = 0;
    double extraTime = 0.0;
    int i = std::max(0, static_cast<int>(osc.m_startTime * sampleFreq + 0.5));
    double endTime = osc.m_endTime > 1000 ? osc.m_startTime + .001 : osc.m_endTime;
    for (;
         i < std::min<int>(buffer.cols(), static_cast<int>((endTime + extraTime) * sampleFreq + 0.5)) ||
         (state.norm() > 1e-15 && i < buffer.cols());
         //i < buffer.cols();
         ++i, t += dt, ++counter)
    {
        state = integrator.step(state,
                                t,
                                dt,
                                bubOsc);

        //double tmp = std::exp(-bubOsc.curBeta * (t - osc.m_startTime));
        //double blend = tmp >= 0.85 ? std::exp( -std::pow(tmp - 0.85, 2) / 0.0028125) : 1;

        buffer(i) = state(1) * osc.m_transfer.interp(t);
        //buffer(i) = state(1);


        //buffer(i) = state(1);
        //buffer(i) = state(1) * osc.m_transfer.interp(t) * blend;
        //buffer(i) = state(1) * bub.m_transfer.linearInterp(t);
    }

    // Add a pop sound at the end of collapsed oscillators

    if (1 && osc.m_endType == Bubble::COLLAPSE)
    {
        const TimeSeries<double>::DataVector &radiiData = osc.m_radii.getData();
        //Eigen::ArrayXd popSound = pop.makePop(radiiData(radiiData.rows()-1, 1), sampleFreq, 0.004);
        Eigen::ArrayXd popSound = pop.makePop(radiiData(radiiData.rows()-1, 1), sampleFreq);

        if (buffer.cols() - i > popSound.rows())
        {
            double bubMag = buffer.array().abs().maxCoeff();

            // Add the scaled pop sound
            buffer.array().segment(i, popSound.rows()) += popSound * bubMag * s_popAmp(s_rnd);
            //buffer.array().segment(i, popSound.rows()) += popSound * bubMag * .03;
        }
    }
}

template<int ROWS, int COLS>
void
writeSoundFile(const std::string &name,
               const Matrix<double, ROWS, COLS>& data,
               int sampleRate)
{
    SF_INFO info;
    info.frames = data.cols();
    info.samplerate = sampleRate;
    info.channels = data.rows();
    info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    //info.sections = 1;
    //info.seekable = 0;

    //std::cout << "frames: " << info.frames << " channels: " << info.channels << std::endl;

    SNDFILE* sf = sf_open(name.c_str(),
                          SFM_WRITE,
                          &info);

    if (!sf)
    {
        std::cerr << "Error: " << sf_strerror(NULL) << std::endl;
        throw std::runtime_error("error opening sound file for writing");
    }

    sf_count_t nWritten = sf_writef_double(sf, data.data(), data.cols());

    if (nWritten != data.cols())
    {
        throw std::runtime_error("incorrect number of frames written");
    }

    // TODO: check return code
    (void)sf_close(sf);
}

static
void
outputOscillatorInfo(const std::vector<Oscillator>& oscillators)
{
    std::ofstream of("oscillatorBubInfo.txt");

    for (size_t i = 0; i < oscillators.size(); ++i)
    {
        of << "Bub " << std::setw(8) << std::setfill('0') << i << std::endl;

        const Oscillator &osc = oscillators.at(i);

        of << osc.m_startTime << std::endl;
        of << osc.m_endTime << std::endl;

        for (auto j : osc.m_bubbleNumbers)
        {
            of << j << std::endl;
        }

        for (auto j : osc.m_trackedBubbleNumbers)
        {
            of << j.second << " " << j.first << std::endl;
        }

        of << std::fixed << std::setprecision(6);

        for (int i = 0; i < osc.m_radii.getData().rows(); ++i)
        {
            of << osc.m_radii.getData()(i, 0) << " " << osc.m_radii.getData()(i,1) << " "
               << osc.m_frequencies.getData()(i,1) / 2 / M_PI << " "
               << std::setw(12) << std::setfill(' ') << osc.m_transfer.getData()(i, 1) << " " << osc.m_pressure.getData()(i,1) << std::endl;
        }
    }

    of.close();
}

// Class to hold a triangle mesh
class Mesh
{
public:
    enum
    {
        FLUID_AIR = 1,
        SOLID = 2
    }

	std::vector<Eigen::Vector3d> m_vertices;
	std::vector<Eigen::Vector3i> m_triangles;
    std::vector<int> m_triType;

    std::vector<int> m_surfTris; // list of surface triangles (where the velocity solution data is)
    std::vector<Eigen::Vector3d> m_surfTriCenters;

	loadGmsh(const std::string &fileName)
	{
        std::ifstream in(fileName.c_str());

		std::string line;

		// Skip first four lines
		for (int i = 0; i < 4; ++i)
		{
			std::getline(in, line);
		}

		// Next line is # of vertices
		int numVerts;
		in >> numVerts;

		m_vertices.resize(numVerts);

		// Read the vertices
		for (int i = 0; i < numVerts; ++i)
		{
		    int index;
		    double x, y, z;

		    in >> index >> x >> y >> z;

            m_vertices[i] << x, y, z;
		}

		// Skip two lines
		in >> line >> line;

		// Read # of triangles
		int numFaces;
		in >> numFaces;

		m_triangles.resize(numFaces);
		m_triType.resize(numFaces);

		// Read triangles
		for (int i = 0; i < numFaces; ++i)
        {
            int ignore, type, index, v1, v2, v3;
            in >> ignore >> ignore >> ignore;
            in >> type >> index >> v1 >> v2 >> v3;

            m_triangles[i] << v1 - 1, v2 - 1, v3 - 1;
            m_triType[i] = type;

            // TODO: confirm whether solution data is full mesh
            // or just fluid surface data
            if (type == FLUID_AIR || type == RIGID)
            {
                m_surfTris.push_back(i);

                m_surfTriCenters.push_back( 1./3. * (m_vertices[v1] + m_vertices[v2] + m_vertices[v3]) );
            }
        }
	}

	// TODO: implement filter?
};

typedef struct
{
    int bubNum;
    double freq; // Negative for bad/dead bubbles
    std::complex<double> xfrIn; // TODO: unnecessary here
} BubbleInputInfo;

static std::vector<BubbleInputInfo>
parseFreqFile(const std::string &fName)
{
    using namespace std;
    vector<BubbleInputInfo> bubInfo;

    ifstream in(fName.c_str());

    string line;

    bubInfo.clear();

    // First line is time
    getline(in, line);

    // Read first real line
    getline(in, line);

    vector<string> data;
    while (!line.empty())
    {
        boost::split(data, line, boost::is_any_of(" \n"), boost::token_compress_on);

        if (data.size() < 2)
        {
            break;
        }

        BubbleInputInfo curInfo;
        curInfo.bubNum = atoi(data[0].c_str());

        if (data.size() == 2)
        {
            // Bad bubble
            curInfo.freq = -1;
        }
        else
        {
            // Good bubble

            curInfo.freq = atof(data[1].c_str());

            vector<string> xfrData;
            boost::split(xfrData, data[2], boost::is_any_of("(),"), boost::token_compress_on);

            curInfo.xfrIn = complex<double>(atof(xfrData[1].c_str()), atof(xfrData[2].c_str()));
        }

        bubInfo.push_back(curInfo);

        getline(in, line);
    }

    return bubInfo;
}

// indexed by bubble id, then vector of triangle velocities
typedef std::map<int, std::vector<Eigen::Vector1d>> SurfaceVelocityData;
typedef MLSModeInterpolator<double, 3, 1> MLSInterp; // TODO: should this be 1d or 3d? (interpolate normal velocities or full velocity vectors?)
typedef KDTree<3, Eigen::Vector1d, Dist> PointKDTree;

Mesh m1, m2;
SurfaceVelocityData v1, v2;
MLSInterp mls;
std::shared_ptr<PointKDTree> kd1, kd2; // TODO: add copy/move semantics to the kd tree class so shared pointers aren't necessary
std::vector<BubbleInputInfo> b1, b2;
double t1, t2; // surrounding times for surface data
Eigen::VectorXd velT1, velT2;


struct FileNames
{
    std::string meshFile;
    std::string datFile; // helmholtz solution file
    std::string freqFile; // frequency info file
};

std::map<double, FileNames> fileInfo; // indexed by time

SurfaceVelocityData
loadSurfaceDatFile(const std::vector<BubbleInputInfo> &bubInfo,
                   const std::string &fileName,
                   const Mesh &mesh)
{
    std::ifstream in(fileName.c_str());
    SurfaceVelocityData output;

    // Count number of fluid surface triangles
    int airTris = 0;

    for (int type : mesh.m_triType)
    {
        if (type == Mesh::FLUID_AIR)
            ++airTris;
    }

    double rho = 1.184; // Density of air

    for (int i = 0; i < bubInfo.size(); ++i)
    {
        // skip bad bubbles
        if (bubInfo[i].freq < 0) continue;

        double omega = 2 * M_PI * bubInfo[i].freq;

        std::complex<double> factor(0, -omega * rho);
        std::complex<double> val;

        // Read the dirichlet data and discard it
        for (int j = 0; j < mesh.m_vertices.size(); ++j)
        {
            in.read((char*)&val.real(), sizeof(val.real()));
            in.read((char*)&val.imag(), sizeof(val.imag()));
        }

        output[bubInfo[i].bubNum].resize(airTris);

        std::vector<Eigen::Vector1d> &curOutput = output[bubInfo[i].bubNum];

        // Read the velocity data and store it
        for (int j = 0; j < airTris; ++j)
        {
            in.read((char*)&val.real(), sizeof(val.real()));
            in.read((char*)&val.imag(), sizeof(val.imag()));

            // Convert from pressure gradient to velocity here
            curOutput.at(j) << (val / factor).real();
        }
    }

    return output;
}

void
updateOscillators(REAL time)
{
    RK4<Vector2d> integrator;
    typedef BubbleOscillator<Vector2d> DS;

    for (int i = 0; i < m_oscillators.size(); ++i)
    {
        Oscillator &osc = m_oscillators[i];

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
                                          dt,
                                          bubOsc);

            osc.m_lastVals(2) = osc.m_state(1);
        }

        while (osc.m_currentTime < time)
        {
            osc.m_state = integrator.step(osc.m_state,
                                          time,
                                          dt,
                                          bubOsc);

            osc.m_lastVals(0) = osc.m_lastVals(1);
            osc.m_lastVals(1) = osc.m_lastVals(2);
            osc.m_lastVals(2) = osc.m_state(1);

            osc.m_currentTime += dt;
        }
    }
}

void
computeVelocities(REAL time)
{
    using namespace std;
    using namespace Eigen;

    bool useT1 = std::fabs(time - t1) < std::fabs(t2 - time);

    // Interpolate onto the mesh that is closest
    Mesh &m = useT1 ? m1 : m2;

    double localT1 = time - dt;
    double localT2 = time + dt;

    velT1 = VectorXd::Zero(m.m_surfTris.size());
    velT2 = VectorXd::Zero(m.m_surfTris.size());

    std::vector<int> closest1, closest2;

	// Loop through the oscillators
    for (int i = 0; i < m_oscillators.size(); ++i)
    {
        Oscillator &osc = m_oscillators[i];

        if (!osc.isActive(time)) continue;

        bool existT1 = osc.m_startTime <= t1;
        bool existT2 = osc.m_startTime <= t2;

        if (!existT1 && !existT2) continue;

        std::shared_ptr<PointKDTree> tree1(existT1 ? kd1 : kd2);
        std::shared_ptr<PointKDTree> tree2(existT2 ? kd2 : kd1);

        SurfaceVelocityData &vel1 = existT1 ? v1 : v2;
        SurfaceVelocityData &vel2 = existT2 ? v2 : v1;

        Mesh &localM1 = existT1 ? m1 : m2;
        Mesh &localM2 = existT2 ? m2 : m1;

        auto iter = osc.m_trackedBubbleNumbers.lower_bound(t1 - 1e-15);
        if (iter == osc.m_trackedBubbleNumbers.end())
        {
            throw runtime_error("bad tracked bubble numbers lookup t1");
        }

        int bubbleNumber1 = iter->second;

        iter = osc.m_trackedBubbleNumbers.lower_bound(t2 - 1e-15);
        if (existT2 && iter == osc.m_trackedBubbleNumbers.end())
        {
            throw runtime_error("bad tracked bubble numbers lookup t2");
        }

        int bubbleNumber2 = iter->second;

        if (!existT1)
        {
            bubbleNumber1 = bubbleNumber2;
        }
        else if (!existT2)
        {
            bubbleNumber2 = bubbleNumber1;
        }

        for (int j = 0; j < m.m_surfTris.size(); ++j)
        {
            // Now interpolate to correct times
            MLSVal val1, val2;
            MLSPoint p = m.m_surfTriCenters[j];

            tree1.find_nearest(p,
                               6,
                               closest1);

            tree2.find_nearest(p,
                               6,
                               closest2);

            // Values at t1 and t2
            val1 = mls.lookup(p,
                              localM1.m_surfTriCenters,
                              vel1.at(bubbleNumber1),
                              -1,
                              NULL,
                              closest1);

            val2 = mls.lookup(p,
                              localM2.m_surfTriCenters,
                              vel2.at(bubbleNumber2),
                              -1,
                              NULL,
                              closest2);

            // Now interpolate
            double pct = (localT1 - t1) / (t2 - t1);
            velT1(j) += osc.m_lastVals(0) * ( pct * val2 + (1 - pct) * val1 );

            pct = (localT2 - t1) / (t2 - t1);
            velT2(j) += osc.m_lastVals(2) * ( pct * val2 + (1 - pct) * val1 );
        }
	}
}

std::shared_ptr<PointKDTree>
createKDTree(const Mesh &m)
{
    std::shared_ptr<PointKDTree> tree(new PointKDTree(m.m_surfTriCenters.data(), m.m_surfTriCenters.size(), false));

    return tree;
}

void
updateSurfaceVelocityData(REAL time)
{
    if (time <= t2)
    {
        return;
    }

    auto iter = fileInfo.lower_bound(time);

    if (iter == fileInfo.end())
    {
        // arghhh
        throw std::runtime_error("handle this");
    }

    t1 = t2;
    t2 = iter->first;

    m1 = m2;
    m2.loadGmsh(iter->second.meshFile);

    b1 = b2;
    b2 = parseFreqFile(iter->second.freqFile);

    v1 = v2;
    v2 = loadSurfaceDatFile(b2,
                            iter->second.datFile,
                            m2);

    kd1 = kd2;
    kd2 = createKDTree(m2);
}

void
projectToPlane()
{
}

// Run this once per timestep
void
updateTime(REAL time)
{
    // Update all oscillators to the correct time (one timestep after this time)
    updateOscillators(time);

    // Load new surface velocity data if necessary
    updateSurfaceVelocityData(time);

    // Compute velocity before and after current time step
    computeVelocities(time);

    // Compute acceleration
    Eigen::VectorXd accel = (velT2 - velT1) / (2 * dt);

    // Project to plane
    // This step is only necessary until the wavesolver handles deforming geometry
    projectToPlane();
}

void
getAcceleration(vec3 position, vec3 normal, REAL time)
{
    // Update time step if necessary
    updateTime(time);

    // Get velocity before and after current time using MLS interpolators

    // Finite difference to get acceleration
}

int
main (int argc,
      char **argv)
{
    if (argc != 4)
    {
        std::cout << "Usage: " << argv[0] << " <freq type c=capaticance m=minnaert> timeLength <bubble info file>" << std::endl;
        return 1;
    }

    std::string fType(argv[1]);

    FreqType ft;

    if (fType == "c")
    {
        ft = CAPACITANCE;
    }
    else if (fType == "m")
    {
        ft = MINNAERT;
    }
    else
    {
        std::cerr << "Invalid frequency type: " << fType << std::endl;
        return 1;
    }


    sampleFreq = 192000;
    dt = 1. / sampleFreq;

    std::map<int, Bubble> bubbles = parseConfigFile(argv[3], ft);
    std::vector<Oscillator> oscillators = makeOscillators(bubbles);

    // Load meshes and solution data
    // List of time steps, bubble solution data at each time step

    BubbleSurfaceDataIndex bubbleSurfaceData = loadSurfaceData();

    // Total time to step for
    timeLength = std::atof(argv[2]);
    buffer.resize(static_cast<int>(timeLength*sampleFreq + 0.5));

    RowVectorXd totalBuffer = RowVectorXd::Zero(static_cast<int>(timeLength*sampleFreq + 0.5));

    double maxVal = 0;
    int maxInd = 0;

    for ( int i = 0; i < static_cast<int>(oscillators.size()); ++i )
    {
        std::cout << "\rOsc " << i << " of " << oscillators.size() << "    ";
        std::cout.flush();

        //if ( i != 80793 ) continue;

        buffer.setZero();

        runOscillator(i,
                      oscillators[i]);

        // Now save the buffer if desired
        std::ostringstream os;
        os << "bubSound-" << std::setw(5) << std::setfill('0') << i << ".wav";

        double n = buffer.lpNorm<Infinity>();

        if (n > 0)
        {
            RowVectorXd newBuffer  = buffer / n * 1.01;
            if (0)
            {
                writeSoundFile(os.str(),
                               newBuffer,
                               sampleFreq);
            }

            if (n > maxVal)
            {
                maxVal = n;
                maxInd = i;
            }

            // And add to the total buffer
            totalBuffer += buffer;
        }
    }
    std::cout << std::endl;

    std::cout << "maxVal: " << maxVal << " maxInd: " << maxInd << std::endl;

    // Write the total buffer
    double n = totalBuffer.lpNorm<Infinity>();
    //n = 1.47714e-08;
    if (n > 0)
    {
        totalBuffer /= n * 1.01;
        writeSoundFile("totalSound.wav",
                       totalBuffer,
                       sampleFreq);
    }


    return 0;
}

