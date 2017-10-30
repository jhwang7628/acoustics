#ifndef _SPLINE_DERIVATIVES_HPP
#define _SPLINE_DERIVATIVES_HPP

#include <Eigen/Dense>

/**
 * Calculates one sided forward differences at the end
 * points.
 */
template<typename Scalar>
class OneSideFD
{
public:
    typedef Eigen::Matrix<Scalar, 1, 2> DiffMat;

    enum
    {
        UseEndPoints = 1
    };

    Scalar
    operator() (const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data,
                int index) const
    {
        assert (index == 0 || index == data.rows() - 1);

        if (index == 0)
        {
            DiffMat diff = data.row(index + 1) - data.row(index);
            return diff(1) / diff(0);
        }
        else
        {
            DiffMat diff = data.row(index) - data.row(index - 1);

            return diff(1) / diff(0);
        }
    }
};

/**
 * Sets end point derivatives to 0
 */
template<typename Scalar>
class ZeroEndPt
{
public:
    enum
    {
        UseEndPoints = 1
    };

    Scalar
    operator() (const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data,
                int index) const
    {
        assert (index == 0 || index == data.rows() - 1);
        return 0;
    }
};

/**
 * Does not evaluate derivatives at ends of curve
 * (default for Catmull-Rom)
 */
template<typename Scalar>
class
NoEndPt
{
public:
    enum
    {
        UseEndPoints = 0
    };

    Scalar
    operator() (const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data,
                int index) const
    {
        throw std::runtime_error("Cannot evaluate end segments with this spline");
    }
};

/**
 * Uses finite differences to set derivatives
 * of the splines
 */
template<typename Scalar, typename EndPt = OneSideFD<Scalar> >
class SplineFD
{
public:
    enum
    {
        UseEndPoints = EndPt::UseEndPoints
    };

    typedef Eigen::Matrix<Scalar, 1, 2> DiffMat;

    /**
     * Locations in first column, data in second
     */
    Scalar
    operator() (const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data,
                int index) const
    {
        assert (index >= 0 && index < data.rows());

        if (index > 0 && index < data.rows() - 1)
        {
            // Internal point
            DiffMat d1 = data.row(index + 1) - data.row(index);
            DiffMat d2 = data.row(index) - data.row(index - 1);

            return 0.5 * (d1(1) / d1(0) + d2(1) / d2(0));
        }
        else
        {
            // End point
            return endPt(data,
                         index);
        }
    }

private:
    EndPt endPt;
};

/**
 * Creates a Catmull-Rom spline
 */
template<typename Scalar, typename EndPt = NoEndPt<Scalar> >
class CatmullRom
{
public:
    enum
    {
        UseEndPoints = EndPt::UseEndPoints
    };

    typedef Eigen::Matrix<Scalar, 1, 2> DiffMat;

    /**
     * Locations in first column, data in second
     */
    Scalar
    operator() (const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data,
                int index) const
    {
        assert (index >= 0 && index < data.rows());

        if (index > 0 && index < data.rows() - 1)
        {
            // Internal point
            DiffMat d1 = data.row(index + 1) - data.row(index - 1);

            return d1(1) / d1(0);
        }
        else
        {
            // End point
            return endPt(data,
                         index);
        }
    }

private:
    EndPt endPt;
};
#endif // _SPLINE_DERIVATIVES_HPP

