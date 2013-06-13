#include "PardisoMatrix.hpp"
#include <stdio.h>

#include <map>
#ifdef USE_MKL
#   include <mkl.h>
#else
#   if defined(__APPLE__) || defined(MACOSX)
#       include <vecLib/cblas.h>
#   endif
#endif

#include "utils/math.hpp"

template <>
void PardisoMatrix<double>::add(const PardisoMatrix<double>& mat)
{
    if ( m_data.size() != mat.m_data.size() ||
         m_rowIdx.size() != mat.m_rowIdx.size() )
    {
        fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
        exit(1);
    }

    cblas_daxpy(m_data.size(), 1, &mat.m_data[0], 1, &m_data[0], 1);
}

template<>
void PardisoMatrix<double>::axpy(double beta, const PardisoMatrix<double>& mat)
{
    if ( m_data.size() != mat.m_data.size() ||
         m_rowIdx.size() != mat.m_rowIdx.size() )
    {
        fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
        exit(1);
    }

    cblas_daxpy(m_data.size(), beta, &mat.m_data[0], 1, &m_data[0], 1);
}

template<>
bool PardisoMatrix<double>::operator == (const PardisoMatrix<double>& m2) const
{
    if ( m_cols.size() != m2.m_cols.size() ||
         m_rowIdx.size() != m2.m_rowIdx.size() ||
         m_data.size() != m2.m_data.size() )
        return false;

    for(int i = (int)m_cols.size()-1;i >= 0;-- i)
        if ( m_cols[i] != m2.m_cols[i] ) return false;
    for(int i = (int)m_rowIdx.size()-1;i >= 0;-- i)
        if ( m_rowIdx[i] != m2.m_rowIdx[i] ) return false;
    for(int i = (int)m_data.size()-1;i >= 0;-- i)
        if ( M_ABS(m_data[i]-m2.m_data[i]) > 1E-12 ) return false;
    return true;
}

template<>
bool PardisoMatrix<float>::operator == (const PardisoMatrix<float>& m2) const
{
    if ( m_cols.size() != m2.m_cols.size() ||
         m_rowIdx.size() != m2.m_rowIdx.size() ||
         m_data.size() != m2.m_data.size() )
        return false;

    for(int i = (int)m_cols.size()-1;i >= 0;-- i)
        if ( m_cols[i] != m2.m_cols[i] ) return false;
    for(int i = (int)m_rowIdx.size()-1;i >= 0;-- i)
        if ( m_rowIdx[i] != m2.m_rowIdx[i] ) return false;
    for(int i = (int)m_data.size()-1;i >= 0;-- i)
        if ( M_ABS(m_data[i]-m2.m_data[i]) > (float)1E-8 ) return false;
    return true;
}

