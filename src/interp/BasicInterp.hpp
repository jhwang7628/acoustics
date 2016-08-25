#ifndef INTERP_BASIC_INTERP_HPP
#   define INTERP_BASIC_INTERP_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

//##############################################################################
// This accessor applies when all samples are with fixed sampling width. The
// data can be monotonically increasing or decreasing. O(1) search.
//##############################################################################
template <typename _T>
struct FixAcc
{
    int operator() (int n, const _T* xx, _T x) const 
    {
        if ( n < 2 )
            throw std::runtime_error("**ERROR** data size is too short:" + std::to_string(n)); 

        if ( x < fmin(xx[0], xx[n-1]) || x > fmax(xx[0], xx[n-1]) )
            throw std::runtime_error("**ERROR** given value is out of range");

        if (x == xx[n-1])
            return n-2; // consistent with what was returned by BisecAcc
        else 
            return (int)((x - xx[0]) / (xx[1] - xx[0]));
    }
}; 

//##############################################################################
// This accessor applies when samples are monotonically increasing or
// decreasing. The samples can have varying sampling width. This incurs O(log
// n) search cost. 
//##############################################################################
template <typename _T>
struct BisecAcc
{
    int operator() (int n, const _T* xx, _T x) const
    {
        if ( n < 2 )
            throw std::runtime_error("**ERROR** data size is too short:" + std::to_string(n)); 

        if ( x < fmin(xx[0], xx[n-1]) || x > fmax(xx[0], xx[n-1]) )
            throw std::runtime_error("**ERROR** given value is out of range");

        bool ascnd = xx[n-1] >= xx[0];
        int jl = 0, jm, ju = n-1;
        while ( ju - jl > 1 )
        {
            jm = (ju + jl) >> 1;
            if ( (x >= xx[jm]) == ascnd )
                jl = jm;
            else
                ju = jm;
        }
        return jl;
    }
};

template <typename _T, bool _CopyData>
class BasicInterp
{
    public:
        BasicInterp():m_xx(NULL), m_yy(NULL), m_n(0) { }

        ~BasicInterp()
        {
            if ( _CopyData )
            {
                delete []m_xx;
                delete []m_yy;
            }
        }

    protected:
        void init_data(size_t n, const _T* xx, const _T* yy)
        {
            if ( _CopyData )
            {
                delete []m_xx;
                delete []m_yy;

                m_n  = n;
                m_xx = new _T[n];
                m_yy = new _T[n];
                memcpy(const_cast<_T*>(m_xx), xx, sizeof(_T)*n);
                memcpy(const_cast<_T*>(m_yy), yy, sizeof(_T)*n);
            }
            else
            {
                m_xx = xx;
                m_yy = yy;
                m_n  = n;
            }
        }

    protected:
        const _T*   m_xx;
        const _T*   m_yy;
        size_t      m_n;
};

#ifdef USE_NAMESPACE
}
#endif
#endif
