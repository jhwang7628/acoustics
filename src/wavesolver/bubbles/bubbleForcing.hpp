#ifndef _BUBBLE_FORCING_HPP
#define _BUBBLE_FORCING_HPP

#include <random>
#include <memory>
#include <map>
#include <cmath>
#include "timeSeries.hpp"
#include "constants.h"
#include "Bubble.hpp"

static const double fTimeCutoff = 0.0006;

static std::default_random_engine s_forcingRnd;
//static std::normal_distribution<double> s_eta(0.84, 0.2);
//static std::normal_distribution<double> s_frac(0.64, 0.1);
static std::uniform_real_distribution<double> s_eta(0.4, 1.5);
static std::uniform_real_distribution<double> s_frac(0.4, 0.8);

class ForcingFunction
{
public:
    virtual
    ~ForcingFunction ()
    {
    }

    virtual double
    value (double t) = 0;
};

class CzerskiJetForcing : public ForcingFunction
{
public:
    CzerskiJetForcing (double r,
                       double cutoff = std::numeric_limits<double>::infinity(),
                       double eta = ETA,
                       bool useLaplace = false,
                       bool useModulation = false,
                       double multiplier = 1.0)
        : m_r(r),
          //m_cutoff (std::min(cutoff, 0.5 / (3.0 / r))), // 1/2 minnaert period
          m_cutoff (std::min(fTimeCutoff, std::min(cutoff, 0.5 / (3.0 / r)))), // 1/2 minnaert period
          //m_cutoff (std::min(cutoff, 0.0004)), // 1/2 minnaert period
          //m_eta(eta),
          m_eta(s_eta(s_forcingRnd)),
          m_multiplier(multiplier),
          m_useLaplace(useLaplace),
          m_laplaceDone(false),
          m_useModulation(useModulation)
    {
    }

    double
    modulation(double t)
    {
        return 0.5 - 1.0/M_PI * atan(5000.0 * (t - m_r));
        return 0.5 * erfc(24. * t/m_r - 5.);
    }

    double
    value (double t)
    {
        //return fittedValue(t);
        if (m_useLaplace)
        {
            return m_multiplier * (jetValue(t) - laplaceVal(t));
        }
        else
        {
            return m_multiplier * jetValue(t);
        }
        //return gaussianValue(t);
    }

    double
    laplaceVal(double t)
    {
        if (m_useLaplace && !m_laplaceDone)
        {
            m_laplaceDone = true;
            return 2 * SIGMA / m_r;
        }

        return 0;
    }

    double
    jetValue (double t)
    {
        if (t > m_cutoff)
        {
            return 0;
        }

        double jval = -9 * GAMMA * SIGMA * m_eta * (ATM + 2 * SIGMA/m_r) * sqrt(1 + m_eta*m_eta) /
                     (4 * RHO_WATER * RHO_WATER * m_r * m_r * m_r * m_r * m_r) * t*t;

        // Convert to radius (instead of fractional radius)
        jval *= m_r;

        // Convert to pressure
        double mrp = RHO_WATER * m_r;
        jval *= mrp;

        // Convert to force
        //double factor = 4 * M_PI * m_r * m_r;
        //jval *= factor;

        if (m_useModulation)
        {
            return jval * modulation(t);
        }
        else
        {
            return jval;
        }
    }

    double
    gaussianValue (double t)
    {
        double center = m_cutoff;
        double mag = jetValue(center);
        double stdDev = center / 3.0;

        return mag * exp(- (t - center) * (t - center) / (2 * stdDev * stdDev));
    }

    double
    fittedValue (double t)
    {
        if (t > 0.05)
        {
            return 0;
        }

        double fittedR = m_r; //0.002;

        double theta = 40 * M_PI / 180.0;
        double pEq = ATM + 2 * SIGMA/fittedR;
        double top = -9 * GAMMA * SIGMA * pEq * tan(theta) * tan(theta);
        double bottom = 4 * RHO_WATER * RHO_WATER * pow(fittedR, 5) * sin(theta);

        double jval = top / bottom * pow(t, 469.5*t + 1.955);

        // Convert to radius (instead of fractional radius)
        jval *= fittedR;

        // Convert to pressure
        double mrp = RHO_WATER * fittedR;
        jval *= mrp;

        return jval * modulation(t);
    }

private:
    double m_r;
    double m_cutoff;
    double m_eta;

    double m_multiplier;

    bool m_useLaplace;
    bool m_laplaceDone;
    bool m_useModulation;
};

class MergeForcing : public ForcingFunction
{
public:
    // r: new bubble radius
    // r1: parent bubble 1 radius
    // r2: parent bubble 2 radius
    MergeForcing (double r,
                  double r1,
                  double r2,
                  double cutoff = std::numeric_limits<double>::infinity())
        : m_r(r)
    {
        // tlim is the time taken for the expanding radius
        // to reach a fixed fraction of the smaller parent
        // radius. Values in the paper ranged from 0.64 - 0.75

        // Equation 5 from Czerski 2011
        double frac = 0.64;
        frac = s_frac(s_forcingRnd);
        double factor = std::pow(2. * SIGMA * r1 * r2 / (RHO_WATER * (r1 + r2)),
                                 0.25);

        cutoff = std::min(cutoff, 0.5 / (3.0 / r)); // 1/2 minnaert period
        cutoff = std::min(cutoff, fTimeCutoff); // 1/2 minnaert period
        m_tlim = std::min(cutoff,
                          std::pow(frac * std::min(r1, r2) / 2. / factor,
                                   2));
    }

    double
    modulation (double t)
    {
        // From Czerski paper, doesn't go to zero fast enough
        return 0.5 - 1/M_PI * std::atan(3 * (t - m_tlim)/m_tlim);

        //return 0.5 * std::erfc(4000. * (t - m_tlim));
    }

    double
    value (double t)
    {
        if (t > m_tlim)
        {
            return 0;
        }

        double jval = 6 * SIGMA * GAMMA * (ATM + 2 * SIGMA/m_r) * t*t /
                      (RHO_WATER * RHO_WATER * m_r * m_r * m_r * m_r * m_r);

        // Convert to radius (instead of fractional radius)
        jval *= m_r;

        // Convert to pressure
        double mrp = RHO_WATER * m_r;
        jval *= mrp;

        //return jval;
        return jval * modulation(t);
    }

private:
    double m_r;
    double m_tlim;
};

class ZeroForcing : public ForcingFunction
{
public:
    double
    value (double t)
    {
        return 0;
    }
};

class PressureForcing : public ForcingFunction
{
public:
    PressureForcing (const std::vector<double>& times,
                     const std::vector<double>& pVals,
                     double sampleLength)
        : m_length(sampleLength)
    {
        std::vector<double> nTimes;
        for (size_t i = 0; i < times.size(); ++i)
        {
            nTimes.push_back(times[i] - times[0]);
        }

        m_series.setData(Eigen::Map<const Eigen::VectorXd>(nTimes.data(), nTimes.size()),
                         Eigen::Map<const Eigen::VectorXd>(pVals.data(), pVals.size()));
    }

    double
    value (double t)
    {
        //if (t > 0.0015) return 0;

        return -(m_series.interp(t) - m_series.interp(t - m_length));
    }

private:
    double m_length;
    TimeSeries<double> m_series;
};

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

#endif // _BUBBLE_FORCING_HPP

