#ifndef _OSCILLATOR_HPP
#define _OSCILLATOR_HPP

#include <vector>
#include <map>
#include <Eigen/Dense>
#include "timeSeries.hpp"
#include "Bubble.hpp"
#include "bubbleForcing.hpp"

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

#endif // _OSCILLATOR_HPP

