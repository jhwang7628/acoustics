#ifndef _OSCILLATOR_HPP
#define _OSCILLATOR_HPP

#include <vector>
#include <map>
#include <Eigen/Dense>
#include "timeSeries.hpp"

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

#endif // _OSCILLATOR_HPP

