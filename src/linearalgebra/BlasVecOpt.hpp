#ifndef LINEARALGEBRA_BLAS_VEC_OPT_HPP
#   define LINEARALGEBRA_BLAS_VEC_OPT_HPP

#include <assert.h>
#ifdef USE_MKL
#   include <mkl.h>
#else
#   if defined(__APPLE__) || defined(MACOSX)
#       include <vecLib/cblas.h>
#   endif
#endif
#include <vector>
#include <math.h>
#include "Vector3.hpp"

struct BlasVecOpt
{
    /*
     * dst = src
     */
    static inline void copy_from(std::vector<double>& dst, const std::vector<double>& src)
    {
        dst.resize(src.size());
        cblas_dcopy((int)src.size(), &src[0], 1, &dst[0], 1);
    }

    static inline void copy_from(double* dst, const double* src, int n)
    {
        cblas_dcopy(n, src, 1, dst, 1);
    }

    static inline void copy_from(float* dst, const float* src, int n)
    {
        cblas_scopy(n, src, 1, dst, 1);
    }

    static inline void scale(std::vector<double>& vec, double s)
    {
        cblas_dscal(vec.size(), s, &vec[0], 1);
    }

    static inline void scale(std::vector< Vector3<double> >& vec, double s)
    {
        cblas_dscal(vec.size()*3, s, (double *)&vec[0], 1);
    }

    static inline void substract_from(std::vector< Vector3<double> >& out, 
                                      const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, -1, (double *)&vec[0], 1, (double *)&out[0], 1);
    }

    /*
     * out = out + vec
     */
    static inline void add_to(std::vector< Vector3<double> >& out, 
                              const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, 1, (double *)&vec[0], 1, (double *)&out[0], 1);
    }

    static inline double norm(const std::vector< Vector3<double> >& vec)
    {
        return cblas_dnrm2(vec.size()*3, (const double *)&vec[0], 1);
    }

    static inline void axpy(std::vector< Vector3<double> >& out, double alpha, 
                     const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, alpha, (const double *)&vec[0], 1, (double *)&out[0], 1);
    }
};

#endif
