#ifndef INTERP_CUBIC_SPLINE_HPP
#   define INTERP_CUBIC_SPLINE_HPP

#include <vector>
#include "BasicInterp.hpp"

#include <cassert>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

/*
 * cubic spline interpolation (based on Numeric Recipes 3rd Ed)
 */
template <typename _T, bool _CopyData,
          class _Acc = BisecAcc<_T> >
class CSpline : public BasicInterp<_T, _CopyData>
{
    public:
        CSpline() { }

        CSpline(size_t n, const _T* vx, const _T* vy, _T yp1 = 1e+99, _T ypn = 1e+99)
        {
            assert(n > 1);
            init(n, vx, vy, yp1, ypn);
        }

        CSpline(const std::vector<_T>& vx, const std::vector<_T>& vy, 
                _T yp1 = 1e+99, _T ypn = 1e+99)
        {
            assert(vx.size() == vy.size());
            init(vx.size(), &vx[0], &vy[0], yp1, ypn);
        }

        /*
         * evaluated the interpolated value
         */
        _T eval(_T x) const;
        /*
         * compute the second derivative at each tabulated value x_i
         * if yp1 or ypn is large than 0.99e+99, then the second derivative at
         * the boundary is zero, leading to the natural cubic spline
         *
         * yp1, ypn are the first derivative values at boundary points
         */
        void init(size_t n, const _T* xv, const _T* yv, _T yp1 = 1e+99, _T ypn = 1e+99);

    private:
        using BasicInterp<_T, _CopyData>::m_xx;
        using BasicInterp<_T, _CopyData>::m_yy;
        using BasicInterp<_T, _CopyData>::m_n;

        std::vector<_T>  m_yp2;  // the second derivative at each tabulated value x_i
        _Acc             _acc;
};

///////////////////////////////////////////////////////////////////////////////

template <typename _T, bool _CopyData, class _Acc>
_T CSpline<_T, _CopyData, _Acc>::eval(_T x) const
{
    int klo = _acc(m_n, m_xx, x);
    int khi = klo + 1;
    _T  h   = m_xx[khi] - m_xx[klo];

    if ( fabs(h) < 1E-10 ) 
    {
        fprintf(stderr, "x value should be monotonically increased");
        exit(1);
    }

    _T invh = (_T)1 / h;
    _T a = (m_xx[khi] - x) * invh;
    _T b = (x - m_xx[klo]) * invh;

    return a*m_yy[klo] + b*m_yy[khi] + ((a*a*a-a)*m_yp2[klo]
            + (b*b*b-b)*m_yp2[khi])*(h*h)/6.0;
}

template <typename _T, bool _CopyData, class _Acc>
void CSpline<_T, _CopyData, _Acc>::init(
        size_t n, const _T* xv, const _T* yv, _T yp1, _T ypn)
{
    //// init data in BasicInterp
    this->init_data(n, xv, yv);

    std::vector<_T> u(n);    // aux vector
    m_yp2.resize(n);

    //// lower boundary condition
    if ( yp1 > 0.99e+99 )
        m_yp2[0] = u[0] = 0;
    else
    {
        m_yp2[0] = -0.5;
        u[0] = ((_T)3/(xv[1]-xv[0])) * ((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
    }

    _T p, qn, sig, un;
    //// tridiagonal decomposition
    for(size_t i = 1;i < n-1;++ i) 
    {
        sig   = (xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
        p     = sig*m_yp2[i-1] + (_T)2;
        m_yp2[i] = (sig-(_T)1)/p;
        u[i]  = (yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
        u[i]  = ((_T)6*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1]) / p;
    }

    if ( ypn > 0.99e+99 )
        qn = un = 0;
    else
    {
        qn = 0.5;
        un = (3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
    }

    //// update the second derivative
    m_yp2[n-1]=(un-qn*u[n-2]) / (qn*m_yp2[n-2]+(_T)1);
    for(int k = (int)n-2;k >= 0;-- k)
        m_yp2[k] = m_yp2[k]*m_yp2[k+1] + u[k];
}

#ifdef USE_NAMESPACE
}
#endif

#endif
