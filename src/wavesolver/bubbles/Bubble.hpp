#ifndef _BUBBLE_HPP
#define _BUBBLE_HPP

#include <vector>
#include <map>

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

    static EventType
    parseEventType(char type)
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

#endif // _BUBBLE_HPP

