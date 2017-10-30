#ifndef _MITCHELL_NETRAVALI_HPP
#define _MITCHELL_NETRAVALI_HPP

#include <Eigen/Dense>

/**
 * A 1d Mitchell-Netravali cubic filter
 */
template<typename Scalar>
class MitchellNetravali
{
public:
    MitchellNetravali(Scalar B = 1./3.,
                      Scalar C = 1./3.)
        : a1(12 - 9*B - 6*C),
          a2(-18 + 12*B + 6*C),
          a3(6 - 2*B),
          a4(-B - 6*C),
          a5(6*B + 30*C),
          a6(-12*B - 48*C),
          a7(8*B + 24*C)
    {
    }

    /**
     * Interpolates the data at point p
     * Points (location) are in first column of data
     * The data values are in the second column
     */
    Scalar
    operator() (Scalar p,
                const Matrix<Scalar, Dynamic, 2> &data) const
    {
        assert(data.rows() > 4);

        int index = findIndex(p,
                              data);

        if (index == data.rows() - 1)
        {
            return data(data.rows() - 1, 1);
        }

        int start = std::max<int>(index - 1, 0);
        int end = std::min<int>(index + 2, data.rows() - 1);

        Scalar sum = 0;
        Scalar ts = data(1,0) - data(0,0);

        for ( ; start <= end; ++start)
        {
            Scalar d = (p - data(start, 0)) / ts;
            Scalar val = data(start, 1);
            Scalar fVal = filterVal(d);
            sum += val * fVal;
        }

        return sum / 6.;
    }

    Scalar
    filterVal(Scalar x) const
    {
        Scalar fx = fabs(x);
        Scalar fx2 = fx * fx;
        Scalar fx3 = fx2 * fx;

        if (fx > 2)
        {
            return 0;
        }

        if (fx < 1)
        {
            return a1 * fx3 + a2 * fx2 + a3;
        }
        else
        {
            return a4 * fx3 + a5 * fx2 + a6 * fx + a7;
        }
    }

private:
    const Scalar a1, a2, a3, a4, a5, a6, a7;

    int
    findIndex(Scalar p,
              const Matrix<Scalar, Dynamic, 2> &data) const
    {
        int index = 0;

        if (data(0, 0) > p)
        {
            throw std::runtime_error("interpolation point before range");
        }

        while (index < data.rows() - 1 && data(index + 1, 0) < p)
        {
            ++index;
        }

        return index;
    }
};

#endif // _MITCHELL_NETRAVALI_HPP

