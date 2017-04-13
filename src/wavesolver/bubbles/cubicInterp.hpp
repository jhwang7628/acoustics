#ifndef _CUBIC_INTERP_HPP
#define _CUBIC_INTERP_HPP

#include <stdexcept>
#include <Eigen/Dense>
#include "splineDerivatives.hpp"

/**
 * A cubic hermite spline interpolator
 *
 * Diff is the function object to calculate derivatives
 */
template<typename Scalar, typename Diff = SplineFD<Scalar> >
class CubicInterp
{
public:
    /**
     * Points are locations, i.e., times
     * and data is the data
     */
    CubicInterp ()
    {
    }

    /**
     * Interpolates the data at point p
     * Points (location) are in first column of data
     * The data values are in the second column
     */
    Scalar
    operator() (Scalar p,
                const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data) const
    {
        int index = findIndex(p,
                              data);

        if (index == data.rows() - 1)
        {
            return data(index, 1);
        }

        if (!diff.UseEndPoints && (index == 0 || index == data.rows() - 2))
        {
            throw std::runtime_error("Cannot use end intervals with this type of interpolation");
        }

        Scalar m0 = diff(data,
                         index);

        Scalar m1 = diff(data,
                         index + 1);

        double diff = data(index + 1, 0) - data(index, 0);
        double t = (p - data(index, 0)) / diff;

        double t2 = t * t;
        double t3 = t * t2;

        return (2 * t3 - 3 * t2 + 1) * data(index, 1) +
               (t3 - 2 * t2 + t)     * diff * m0 +
               (-2 * t3 + 3 * t2)    * data(index + 1, 1) +
               (t3 - t2)             * diff * m1;
    }

private:
    Diff diff;

    /**
     * Return the index of the first point of the index
     * that p is in
     */
    int
    findIndex(Scalar p,
              const Eigen::Matrix<Scalar, Eigen::Dynamic, 2> &data) const
    {
        int i = 0;

        if (data(i, 0) > p)
        {
            throw std::runtime_error("interpolation point before range");
        }

        while (i < data.rows() - 1 && data(i + 1, 0) < p)
        {
            ++i;
        }

        //if (i == data.rows() - 1)
        //{
            //throw std::runtime_error("interpolation point after range");
        //}

        return i;
    }
};

#endif // _CUBIC_INTERP_HPP

