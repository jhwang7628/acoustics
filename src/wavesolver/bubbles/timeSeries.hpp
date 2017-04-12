#ifndef _TIME_SERIES_HPP
#define _TIME_SERIES_HPP

#include <Eigen/Dense>
#include "cubicInterp.hpp"
//#include <mitchellNetravali.hpp>
#include <math/MLSModeInterpolator.hpp>

/**
 * A general time series class, which allows
 * interpolation of data. Assumes that data
 * will be accessed linearly.
 */
template<typename Scalar>
class TimeSeries
{
public:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 2> DataVector;

    TimeSeries()
         : m_mlsRadius(0.006)
         //m_interp(1./3., 1./3.),
    {
        m_mls.setLookupOrder(1);
    }

    template<int Size>
    TimeSeries(const Eigen::Matrix<Scalar, Size, 1> &times,
               const Eigen::Matrix<Scalar, Size, 1> &data)
         : m_mlsRadius(0.006)
         //m_interp(1./3., 1./3.),
    {
        m_mls.setLookupOrder(1);
        m_data.resize(times.rows(), 2);
        m_data.col(0) = times;
        m_data.col(1) = data;

        for (size_t i = 0; i < times.rows(); ++i)
        {
            Point p;
            p << times(i);
            m_times.push_back(p);

            p << data(i);
            m_dataVals.push_back(p);
        }
    }

    template<int Size>
    TimeSeries(const DataVector &data)
        : m_data(data),
          m_mlsRadius(0.006)
          //m_interp(1./3., 1./3.),
    {
        m_mls.setLookupOrder(1);
        for (size_t i = 0; i < data.rows(); ++i)
        {
            Point p;
            p << data(i, 0);
            m_times.push_back(p);

            p << data(i, 1);
            m_dataVals.push_back(p);
        }
    }

    template<typename Expression>
    void
    setData(const Expression &times,
            const Expression &data)
    {
        //int size = std::max(times.rows(),
                            //times.cols());

        int size = times.rows();
        m_data.resize(size, 2);
        m_data.col(0) = times;
        m_data.col(1) = data;

        for (int i = 0; i < m_data.rows(); ++i)
        {
            Point p;
            p << m_data(i, 0);
            m_times.push_back(p);

            p << m_data(i, 1);
            m_dataVals.push_back(p);
        }
    }

    const DataVector&
    getData() const
    {
        return m_data;
    }

    /**
     * Interpolation of data at time t.
     */
    Scalar
    interp(Scalar t, bool print = false) const
    {
        if (t < m_data(0,0))
        {
            return m_data(0,1);
        }

        if (t > m_data(m_data.rows()-1, 0))
        {
            return m_data(m_data.rows()-1, 1);
        }

        //return linearInterp(t, print);
        //return movingAverage(t, print);

        return m_interp(t,
                        m_data);

        typedef typename MLSModeInterpolator<Scalar, 1, 1>::MLSPoint Point;

        Point pos;
        pos << t;
        Point val = m_mls.lookup(pos,
                                 m_times,
                                 m_dataVals,
                                 m_mlsRadius);

        if (print)
        {
            std::cout << pos << " " << val << std::endl;
            for (size_t i = 0; i < m_times.size(); ++i)
            {
                std::cout << m_times[i] << " ";
            }

            std::cout << std::endl;

            for (size_t i = 0; i < m_dataVals.size(); ++i)
            {
                std::cout << m_dataVals[i] << " ";
            }

            std::cout << std::endl;
        }

        return val(0);

        //int end = m_data.cols()-1;

        //if (t > m_data(0,end))
        //{
            //return m_data(1,end);
        //}

        //return spline( (t - m_data(0,0)) / (m_data(0,end) - m_data(0,0)))(1);
    }

    /**
     * Linear interpolation of data at time t.
     */
    Scalar
    linearInterp(Scalar t,
                 bool print = false) const
    {
        int winStart = 0, winEnd = 0;
        setWindow(t, winStart, winEnd);

        double dist = (t - m_data(winStart, 0)) / (m_data(winEnd, 0) - m_data(winStart, 0));//m_times[m_winStart];

        double val = (1.0 - dist) * m_data(winStart, 1) + dist * m_data(winEnd, 1);

        if (print)
        {
            std::cout << m_data.transpose() << std::endl;
            //std::cout << winStart << " " << winEnd << " " << dist << std::endl;
            //std::cout << m_data(winStart, 0) << " " << m_data(winEnd, 0) << " " << val << std::endl;
        }

        return val;
    }

    Scalar
    movingAverage(Scalar t,
                  bool print = false) const
    {
        int winStart = 0, winEnd = 0;
        setWindow(t, winStart, winEnd);

        int winSz = 4;

        if (print)
        {
            std::cout << "Start: " << winStart << " " << winEnd << std::endl;
        }

        //winStart -= std::min<int>(winSz, winStart);
        //winEnd += std::min<int>(winSz, m_data.rows() - winEnd - 1);

        winStart -= winSz;
        winEnd += winSz;

        if (print)
        {
            std::cout << "Adj: " << winStart << " " << winEnd << std::endl;
        }

        int n = winEnd - winStart;
        double sum = 0;
        for (int i = winStart; i < winEnd; ++i)
        {
            if (i < 0)
                sum += m_data(0, 1);
            else if (i >= m_data.rows())
                sum += m_data(m_data.rows()-1, 1);
            else
                sum += m_data(i, 1);
        }

        return sum / n;
    }

    // TODO: different window lengths for different interpolation methods

private:
    DataVector m_data;
    CubicInterp<Scalar, CatmullRom<Scalar, ZeroEndPt<Scalar> > > m_interp;
    //MitchellNetravali<Scalar> m_interp;

    typedef MLSModeInterpolator<Scalar, 1, 1> MLSInterp;
    typedef typename MLSInterp::MLSPoint Point;
    MLSInterp m_mls;

    std::vector<Point> m_times;
    std::vector<Point> m_dataVals;

    const double m_mlsRadius;

    void
    setWindow(Scalar t,
              int &winStart,
              int &winEnd) const
    {
        while (winStart < m_data.rows() - 1 && m_data(winStart + 1, 0) < t)
        {
            ++winStart;
        }

        winEnd = winStart;

        while (winEnd < m_data.rows() - 1 && m_data(winEnd, 0) <= t)
        {
            ++winEnd;
        }
    }
};

#endif // _TIME_SERIES_HPP

