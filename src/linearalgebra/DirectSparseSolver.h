/******************************************************************************
 *  File: DirectSparseSolver.hpp
 *  A linear solver using Intel MKL DSS routines
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
#ifndef LINEARALGEBRA_DIRECT_SPARSE_SOLVER_H
#   define LINEARALGEBRA_DIRECT_SPARSE_SOLVER_H

#include <mkl_dss.h>
#include "PardisoMatrix.hpp"

/*!
 * This is basically a wrap of MKL's Direct Sparse Solver
 *
 * It is used for the case where A is the same, and there are multiple b for
 * the linear system Ax = b. Therefore we can precompute the factorization of
 * A.
 */
class DirectSparseSolver
{
    public:
        DirectSparseSolver():m_opt(MKL_DSS_DEFAULTS), m_allocated(false)
        { }

        ~DirectSparseSolver()
        {
            if ( m_allocated ) dss_delete(m_handle, m_opt);
        }

        /*!
         * This method first free any previous allocated memory
         * - call dss_create to create a new handle for solver
         * - call dss_define_structure to tell the matrix structure
         *   NOTE: only MKL_DSS_SYMMETRIC / MKL_DSS_NON_SYMMETRIC would be passed
         *         to define the matrix structure
         *
         * - call dss_reorder to compute permutation matrix
         * - call dss_factor_real to factorize the given matrix
         *   NOTE: only MKL_DSS_POSITIVE_DEFINITE / MKL_DDS_INDEFINITE would be
         *         passed for factoring
         */
        void load_matrix(const PardisoMatrix<double>* A);

        /* 
         * solve the system with pre-loaded matrix A
         * it solves Ax = b, and the results is in x
         */
        void solve(const double* b, double *x);

        int& options()
        {  return m_opt; }
        int options() const
        {  return m_opt; }

    private:
        _MKL_DSS_HANDLE_t   m_handle;
        int                 m_opt;
        bool                m_allocated;
};
#endif
