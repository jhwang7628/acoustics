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

template <typename _T>
struct BisecAcc
{
    int operator() (int n, const _T* xx, _T x) const
    {
        if ( n < 2 )
        {
            fprintf(stderr, "ERROR: data size is too short %d\n", n);
            exit(1);
        }

        if ( x < fmin(xx[0], xx[n-1]) || x > fmax(xx[0], xx[n-1]) )
        {
            fprintf(stderr, "ERROR: given value is out of range\n");
            exit(1);
        }

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
