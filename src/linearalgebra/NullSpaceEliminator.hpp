/******************************************************************************
 *  File: NullSpaceEliminator.hpp
 *  Classes to project out the null space from a given vector
 *  Copyright (c) 2009 by Changxi Zheng
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
#ifndef LINEARALGEBRA_NS_ELIMINATOR_HPP
#   define LINEARALGEBRA_NS_ELIMINATOR_HPP

#include <vector>
#include <mkl.h>
#include <mkl_lapack.h>
#include <stdlib.h>
#include "utils/macros.h"

/*
 * Given a matrix which spans the null space of a system, this
 * class is to eliminate the projection in that null space for a
 * given vector.
 *
 * The null space is given by a m x n matrix N (where m > n),
 * we compute the QR factorization first:
 *                           N = Q*R
 * where Q is the orthogonal matrix represent the basis of N
 *
 * then given the vector f, we can eliminate the null space by
 *                    f - Q*Q'*f
 *
 * When use it:
 *  1. Create an instance with given M and N
 *  2. call column_ptr() get the data pointer, and initialize the 
 *     matrix data
 *  3. call init() to compute the QR factorization
 *  4. call eliminate() to eliminate the null space for a given
 *     vector
 */
template <typename T>
class QrNSEliminator
{
    public:
        QrNSEliminator(int m, int n):side_('L'), m_(m), n_(n), work_(8)
        { 
            MSG_ASSERT(m > n, "M should be larger than N of the given matrix");
            // allocate the data space
            data_.resize(m_*n_);
        }

        /*
         * initialize it by computing the QR factorization
         */
        void init();

        /*
         * NOTE: make sure both invec and outvec have the length m
         */
        void eliminate(const T* invec, T* outvec);

        T* column_ptr(int icol)
        {   return &data_[icol*m_]; }

    private:
        char                side_;
        int                 m_, n_;
        int                 lwork_;

        std::vector<T>      data_;
        std::vector<T>      tau_;
        std::vector<T>      work_;
};

///////////////////////////////////////////////////////////////////////////////

template <>
void QrNSEliminator<double>::init()
{
    tau_.resize(n_);

    // first query for the size of work space
    int info;
    lwork_ = -1;
    dgeqrf(&m_, &n_, &data_[0], &m_, &tau_[0], &work_[0], &lwork_, &info);
    if ( info )
    {
        fprintf(stderr, "ERROR: dgeqrf failed (query mode): %d parameter is illegal\n",
                -info);
        exit(1);
    }

    lwork_ = (int)work_[0];
    printf("INFO: workspace with length [%d] for dgeqrf\n", lwork_);
    work_.resize(lwork_);

    // compute the QR factorization
    printf("INFO: computing the QR factorization\n");
    dgeqrf(&m_, &n_, &data_[0], &m_, &tau_[0], &work_[0], &lwork_, &info);
    if ( info )
    {
        fprintf(stderr, "ERROR: dgeqrf failed: %d parameter is illegal\n",
                -info);
        exit(1);
    }

    lwork_ = -1;
    // query for the size of work 
    dorgqr(&m_, &n_, &n_, &data_[0], &m_, &tau_[0], &work_[0], &lwork_, &info);
    if ( info )
    {
        fprintf(stderr, "ERROR: dorgqr failed (query mode): %d parameter is illegal\n",
                -info);
        exit(1);
    }

    lwork_ = std::max(m_, (int)work_[0]);
    work_.resize(lwork_);
    printf("INFO: workspace with length [%d]\n", lwork_);

    printf("INFO: compute %d orthogonal basis\n", n_);
    // compute the n_ leading column of Q
    dorgqr(&m_, &n_, &n_, &data_[0], &m_, &tau_[0], &work_[0], &lwork_, &info);
    if ( info )
    {
        fprintf(stderr, "ERROR: dorgqr failed: %d parameter is illegal\n",
                -info);
        exit(1);
    }
}

template <>
void QrNSEliminator<double>::eliminate(const double* invec, double* outvec)
{
    // compute  work_ = Q'*invec
    cblas_dgemv(CblasColMajor, CblasTrans, m_, n_, 1, &data_[0], m_,
                invec, 1, 0, &work_[0], 1);

    // compute outvec = Q*work_
    cblas_dgemv(CblasColMajor, CblasNoTrans, m_, n_, 1, &data_[0], m_,
                &work_[0], 1, 0, outvec, 1); 

    // compute  outvec <- invec - outvec
    cblas_dscal(m_, -1, outvec, 1);    // outvec <- -outvec
    cblas_daxpy(m_, 1, invec, 1, outvec, 1);
}

#endif
