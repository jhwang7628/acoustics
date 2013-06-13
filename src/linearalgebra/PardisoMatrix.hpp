/******************************************************************************
 *  File: PardisoMatrix.hpp
 *  A sparse matrix used for Intel MKL PARDISO and DSS solver
 *  Copyright (c) 2007 by Changxi Zheng
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#ifndef LINEARALGEBRA_PARDISO_MATRIX_HPP
#   define LINEARALGEBRA_PARDISO_MATRIX_HPP

#include <string.h>
#include <map>
#include <assert.h>
#include "DiagonalMatrix.hpp"
#include "Matrix.hpp"
#include "Vector3.hpp"

/*!
 * Sparse Matrix
 * The data structure of this sparse matrix is particularly optimized
 * for intel MKL PARDISO or DSS solver call. 
 *
 * NOTE: it assumes the diagonal elements is always non-zero
 */

struct PardisoMatrixIO;

template <typename T>
class PardisoMatrix
{
    friend struct PardisoMatrixIO;

    public:
        PardisoMatrix(int m, int n):
                m_numCols(n), m_numRows(m), m_isSym(false),
                m_isPosDef(false)
        {  m_rows.resize(m); }

        PardisoMatrix(int s, bool isSym = false):
                m_numCols(s), m_numRows(s), m_isSym(isSym),
                m_isPosDef(false)
        {  m_rows.resize(s); }

        PardisoMatrix(int s, bool isSym, bool isPosDef):
                m_numCols(s), m_numRows(s), m_isSym(isSym),
                m_isPosDef(isPosDef)
        {  m_rows.resize(s); }

        PardisoMatrix():m_numCols(0), m_numRows(0), m_isSym(false),
                m_isPosDef(false)
        {}

        void clear()
        {
            for(size_t i = 0;i < m_rows.size();++ i) m_rows[i].clear();
        }

        void resize(int s, bool isSym = false)
        {
            m_numCols = m_numRows = s;
            m_isSym = isSym;
            m_rows.resize(s);
            clear();
        }

        void resize(int s, bool isSym, bool isPosDef)
        {
            m_numCols = m_numRows = s;
            m_isSym = isSym;
            m_isPosDef = isPosDef;
            m_rows.resize(s);
            clear();
        }

        void add(const PardisoMatrix<T>& mat)
        {
            if ( m_data.size() != mat.m_data.size() ||
                 m_rowIdx.size() != mat.m_rowIdx.size() )
            {
                fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
                exit(1);
            }

            for(size_t i = 0;i < m_data.size();++ i)
                m_data[i] += mat.m_data[i];
        }

        /*!
         * Axpy operator: compute this = this + beta * mat
         * It assumes the this matrix and the mat have the same sparsity pattern.
         */
        void axpy(T beta, const PardisoMatrix<T>& mat)
        {
            if ( m_data.size() != mat.m_data.size() ||
                 m_rowIdx.size() != mat.m_rowIdx.size() )
            {
                fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
                exit(1);
            }

            for(size_t i = 0;i < m_data.size();++ i)
                m_data[i] += beta * mat.m_data[i];
        }

        void axpy(T beta, const DiagonalMatrix<T>& mat)
        {
            if ( mat.size() != m_numCols || m_numCols != m_numRows )
            {
                fprintf(stderr, "ERROR: the size of diagonal matrix is inconsistent with this PardisoMatrix\n");
                exit(1);
            }

#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 8000) shared(beta, mat)
#endif
            for(int i = 0;i < m_numCols;++ i)
                *(m_rows[i][i]) += beta * mat[i];
        }

        void axpy(T beta, const DiagonalMatrix< Vector3<T> >& mat)
        {
            if ( mat.size()*3 != m_numCols || m_numCols != m_numRows )
            {
                fprintf(stderr, "ERROR: the size of diagonal matrix is inconsistent with this PardisoMatrix\n");
                exit(1);
            }

            const T* ptr = (const T*)mat.data();
#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 8000) shared(beta, ptr)
#endif
            for(int i = 0;i < m_numCols;++ i)
                *(m_rows[i][i]) += beta * ptr[i];
        }

        /*! 
         * indicate the mat[nr][nc] is nonzero
         * the benifit is, for example,
         * as long as we don't remesh, and the mesh isn't split during
         * the simulation, the sparity pattern of stiffness matrix 
         * doesn't change.
         */
        void set_nonzero(int nr, int nc);
        void generate_pattern();        // $$TESTED

        void add(int nr, int nc, T v)
        {
            if ( m_isSym && nr > nc )
            {
                /*
                 * Here we assume mat[nr][nc] += v only applies to the
                 * upper trangle part of the matrix when this is a 
                 * symmetric matrix
                 */
                return;
            }
            assert(m_rows[nr].count(nc));
            *(m_rows[nr][nc]) += v;
        }

        /*! make all the elements be zero */
        void zeros()
        {
            memset(&m_data[0], 0, sizeof(T)*m_data.size());
        }

        void multiply(const std::vector< Vector3<T> >& in, 
                      std::vector< Vector3<T> >& out) const;

        void multiply(const T* in, T* out) const;

        bool check_symmetry();

        void dump() const;
        void to_matrix(Matrix<T>& mat) const;

        /* ================ Retrival Methods =============== */
        int num_rows() const
        { return m_numRows; }
        int num_cols() const
        { return m_numCols; }
        int num_nonzeros() const
        { return m_data.size(); }

        const T* data() const
        { return &m_data[0]; }

        const int* col_indices() const
        { return &m_cols[0]; }

        const int* row_indices() const
        { return &m_rowIdx[0]; }

        T* data() 
        { return &m_data[0]; }

        int* col_indices() 
        { return &m_cols[0]; }

        int* row_indices()
        { return &m_rowIdx[0]; }

        bool is_symmetric() const 
        { return m_isSym; }

        bool is_positive_definite() const
        { return m_isPosDef; }

        bool& is_positive_definite() 
        { return m_isPosDef; }

        bool operator == (const PardisoMatrix<T>& m2) const;

    private:
        typedef std::map<int, T*> RowMap;

        int                 m_numCols;
        int                 m_numRows;
        int                 m_numNonZero;
        std::vector<RowMap> m_rows;
        /*! indicate if it is symmetric matrix */
        bool                m_isSym;
        bool                m_isPosDef;

        /* the data structure for Pardiso storage format */
        std::vector<T>      m_data;
        std::vector<int>    m_cols;
        std::vector<int>    m_rowIdx;
};

///////////////////////////////////////////////////////////////////////////////
template <typename T>
bool PardisoMatrix<T>::check_symmetry()
{
    if ( m_numCols != m_numRows ) return false;
    if ( m_isSym ) return true;

    for(int i = 0;i < m_numRows;++ i)
    {
        typename RowMap::iterator end = m_rows[i].end();
        for(typename RowMap::iterator it = m_rows[i].begin();
                it != end;++ it)
        {
            if ( it->first < i ) continue;
            if ( !m_rows[it->first].count(i) || 
                 M_ABS(*m_rows[it->first][i] - *it->second) > 
                 1E-10*M_MAX(M_ABS(*m_rows[it->first][i]), M_ABS(*it->second)) )
                return false;
        }
    }
    return true;
}

/*
 * If the matrix is symmetric, it stores only the upper part of it
 */
template <typename T>
void PardisoMatrix<T>::set_nonzero(int nr, int nc)
{
    if ( m_isSym && nr > nc ) 
        m_rows[nc][nr] = NULL;
    else
        m_rows[nr][nc] = NULL;
}

/*
 * NOTE that in order to ease to work with Fortran routines, the indices
 *      in both m_cols and m_rowIdx are all 1-based.
 */
template <typename T>
void PardisoMatrix<T>::generate_pattern()
{
    m_data.clear();
    m_cols.clear();
    m_rowIdx.clear();

    // # of nonzero elements
    m_numNonZero = 0;
    for(int i = 0;i < m_numRows;++ i)
        m_numNonZero += m_rows[i].size();
    m_data.resize(m_numNonZero);

    for(int i = 0, ptr = 0;i < m_numRows;++ i)
    {
        if ( m_rows.empty() )
        {
            fprintf(stderr, "Warning: With zeros in the %d-th row\n, cannot generate sparsity pattern for PARDISO call", i);
            return;
        }

        typename RowMap::iterator it  = m_rows[i].begin();
        typename RowMap::iterator end = m_rows[i].end();
        m_rowIdx.push_back(ptr + 1);
        for(;it != end;++ it)
        {
            m_cols.push_back(it->first + 1);    // 1-based column index
            it->second = &(m_data[ptr]);
            ++ ptr;
        }
    }

    m_rowIdx.push_back(m_data.size() + 1);      // 1-based row index
}

template <typename T>
void PardisoMatrix<T>::multiply(const std::vector< Vector3<T> >& in, 
                                std::vector< Vector3<T> >& out) const
{
    assert(in.size()*3 == m_numCols);
    if ( out.size() != m_numRows ) out.resize(m_numRows);

    const T* inPtr = (const T*)&in[0];
    T* outPtr = (T*)out[0];
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic) shared(inPtr, outPtr)
#endif
    for(int i = 0;i < m_numRows;++ i)
    {
        outPtr[i] = 0;
        typename RowMap::const_iterator end = m_rows[i].end();
        for(typename RowMap::const_iterator it = m_rows[i].begin();
                it != end;++ it)
            outPtr[i] += inPtr[it->first]*(*(it->second));
    }
}

template <typename T>
void PardisoMatrix<T>::multiply(const T* in, T* out) const
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic) shared(inPtr, outPtr)
#endif
    for(int i = 0;i < m_numRows;++ i)
    {
        out[i] = 0;
        typename RowMap::const_iterator end = m_rows[i].end();
        for(typename RowMap::const_iterator it = m_rows[i].begin();
                it != end;++ it)
            out[i] += in[it->first]*(*(it->second));
    }
}

template <typename T>
void PardisoMatrix<T>::to_matrix(Matrix<T>& mat) const
{
    mat.zeros();

    int rowId = 0;
    for(size_t i = 0;i < m_data.size();++ i)
    {
        int colId = m_cols[i] - 1;
        if ( m_rowIdx[rowId+1] == i+1 ) ++ rowId;
        mat.set(colId, rowId, m_data[i]);
    }
}

template <typename T>
void PardisoMatrix<T>::dump() const
{
    using namespace std;

    cout << "data: { ";
    for(int i = 0;i < m_data.size();++ i)
        cout << m_data[i] << ' ';
    cout << '}' << endl;
    cout << "col: { ";
    for(int i = 0;i < m_cols.size();++ i)
        cout << m_cols[i] << ' ';
    cout << '}' << endl;
    cout << "row: { ";
    for(int i = 0;i < m_rowIdx.size();++ i)
        cout << m_rowIdx[i] << ' ';
    cout << '}' << endl;
}

template <>
void PardisoMatrix<double>::add(const PardisoMatrix<double>& mat);

template<>
void PardisoMatrix<double>::axpy(double beta, const PardisoMatrix<double>& mat);

template<>
bool PardisoMatrix<double>::operator == (const PardisoMatrix<double>& m2) const;

template<>
bool PardisoMatrix<float>::operator == (const PardisoMatrix<float>& m2) const;

#endif
