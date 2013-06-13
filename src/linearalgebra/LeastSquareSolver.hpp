#ifndef LEAST_SQUARE_SOLVER
#   define LEAST_SQUARE_SOLVER

#include <vector>
#include <mkl.h>
#include <mkl_lapack.h>
#include <gsl/gsl_linalg.h>

#include "generic/trivial_type.hpp"

/*!
 * This header file contains the classes for solving linear least-square 
 * problem.
 */
template<typename T>
class NormalEquSolver
{
    public:
        typedef T       value_type;

        void solve(const T* A, const T* b, T* x);
        void solve(const T* U, const T* A, const T* b, T* x);

        int m() { return m_m; }
        int n() { return m_n; }
        void set_size(int m, int n)
        {
            m_m = m;
            m_n = n;
        }

    private:
        int m_m;
        int m_n;
};

/*!
 * Least-square using normal equation with double precison complex values.
 *
 * NOTE: Both U and A should be stored in row major
 */
template<>
void NormalEquSolver< std::complex<double> >::
solve(const std::complex<double>* U, 
      const std::complex<double>* A, 
      const std::complex<double>* b, 
      std::complex<double>* x)
{
    const std::complex<double> alpha(1, 0);
    const std::complex<double> beta(0, 0);
    
    std::complex<double>  r[m_n]; 
    std::complex<double> NM[m_n][m_n];

    // construct the normal equation
    // r <-- U' * b    U: nrow x width 
    cblas_zgemv(CblasRowMajor, CblasConjTrans, m_m, m_n, &alpha,
                U, m_n, b, 1, &beta, r, 1);
    // NM <-- U' * A 
    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, m_n, m_n,
                m_m, &alpha, U, m_n, A, m_n, &beta, NM, m_n);

    // solve it using gsl routine
    int s;
    gsl_matrix_complex_view m   = gsl_matrix_complex_view_array((double*)NM, m_n, m_n);
    gsl_vector_complex_view rhs = gsl_vector_complex_view_array((double*)r, m_n);
    gsl_vector_complex_view sol = gsl_vector_complex_view_array((double*)x, m_n);
    gsl_permutation*        p   = gsl_permutation_alloc(m_n);

    gsl_linalg_complex_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_complex_LU_solve(&m.matrix, p, &rhs.vector, &sol.vector);
    gsl_permutation_free(p);
}

///////////////////////////////////////////////////////////////////////////////

/*!
 * Ridge regression LS solver
 */
template<typename T>
class RidgeRegrSolver
{
    public:
        typedef T       value_type;
        typedef typename carbine::TrivialType<T>::type trivial_type;

        RidgeRegrSolver():m_eps((trivial_type)1E-8)
        {}

        void solve(const T* A, const T* b, T* x);
        void solve(const T* U, const T* A, const T* b, T* x);

        int m() { return m_m; }
        int n() { return m_n; }
        void set_size(int m, int n)
        {
            m_m = m;
            m_n = n;
        }
        trivial_type& epsilon() { return m_eps; }

    private:
        int         m_m;
        int         m_n;
        trivial_type m_eps;
};

template<>
void RidgeRegrSolver< std::complex<double> >::
solve(const std::complex<double>* U, 
      const std::complex<double>* A, 
      const std::complex<double>* b, 
      std::complex<double>* x)
{
    const std::complex<double> alpha(1, 0);
    const std::complex<double> beta(0, 0);
    
    std::complex<double>  r[m_n]; 
    std::complex<double> NM[m_n][m_n];

    // construct the normal equation
    // r <-- U' * b    U: nrow x width 
    cblas_zgemv(CblasRowMajor, CblasConjTrans, m_m, m_n, &alpha,
                U, m_n, b, 1, &beta, r, 1);
    // NM <-- U' * A 
    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, m_n, m_n,
                m_m, &alpha, U, m_n, A, m_n, &beta, NM, m_n);
    // compute the Frobenius norm of NM
    double nrm = cblas_dznrm2(m_n*m_n, NM, 1);
    // add epsilon on the diagonal elements
    for(int i = 0;i < m_n;++ i)
        NM[i][i].real() += m_eps * nrm;

    // solve it using gsl routine
    int s;
    gsl_matrix_complex_view m   = gsl_matrix_complex_view_array((double*)NM, m_n, m_n);
    gsl_vector_complex_view rhs = gsl_vector_complex_view_array((double*)r, m_n);
    gsl_vector_complex_view sol = gsl_vector_complex_view_array((double*)x, m_n);
    gsl_permutation*        p   = gsl_permutation_alloc(m_n);

    gsl_linalg_complex_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_complex_LU_solve(&m.matrix, p, &rhs.vector, &sol.vector);
    gsl_permutation_free(p);
}

///////////////////////////////////////////////////////////////////////////////
/*!
 * Solve least-square using QR factorization.
 * [  A ][x] = [b]
 * we assume the final matrix is full rank
 *
 * NOTE: the size of A should be at least m*n, and b should be at least m
 */
template <typename T>
class LsQrSolver
{
    public:
        typedef T        value_type;
        typedef typename carbine::TrivialType<T>::type trivial_type;

        LsQrSolver():m_m(0), m_n(0), m_lwork(0)
        {
            m_workspace.resize(2);
        }
        void set_size(int m, int n, int nrhs = 1);

        void solve(T* A, T* b);

    private:
        int         m_m;
        int         m_n;
        int         m_nrow;
        int         m_nrhs;
        int         m_lwork;

        std::vector<value_type>     m_workspace;
};

template<>
void LsQrSolver< std::complex<double> >::
set_size(int m, int n, int nrhs)
{
    m_m = m;
    m_n = n;
    m_nrhs = nrhs;
    m_nrow = m;

    // call the routine in query mode to get the optimal size of workspace
    std::complex<double> *a = NULL;
    std::complex<double> *b = NULL;

    int MINUS_ONE = -1;
    int ret;
    char TRANS = 'N';

    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)a, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &MINUS_ONE,
          &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgels (query mode) %d\n", ret);
        exit(1);
    }

    // allocate workspace
    m_workspace.resize((int)m_workspace[0].real());
    m_lwork = m_workspace.size();
}

template<>
void LsQrSolver< std::complex<double> >::
solve(std::complex<double>* A,
      std::complex<double>* b)
{
    int ret;
    char TRANS = 'N';
    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)A, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &m_lwork,
          &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgels (solve mode) %d\n", ret);
        exit(1);
    }

    //MKL_FreeBuffers();  /* free memory */
}
//-------------------------------------------------------------------------------------------------

template<>
void LsQrSolver<double>::set_size(int m, int n, int nrhs)
{
    m_m = m;
    m_n = n;
    m_nrhs = nrhs;
    m_nrow = m;

    // call the routine in query mode to get the optimal size of workspace
    double *a = NULL;
    double *b = NULL;

    int MINUS_ONE = -1;
    int ret;
    char TRANS = 'N';

    dgels(&TRANS, &m_m, &m_n, &m_nrhs, 
          a, &m_m, b, &m_m, &m_workspace[0], 
          &MINUS_ONE, &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on dgels (query mode) %d\n", ret);
        exit(1);
    }

    // allocate workspace
    m_workspace.resize((int)m_workspace[0]);
    m_lwork = m_workspace.size();
}

template<>
void LsQrSolver<double>::solve(double* A, double* b)
{
    int ret;
    char TRANS = 'N';
    dgels(&TRANS, &m_m, &m_n, &m_nrhs,
          A, &m_m, b, &m_m, &m_workspace[0], 
          &m_lwork, &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on dgels (query mode) %d\n", ret);
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
/*!
 * Solve least-square using Tikhonov regularization and QR.
 * [  A ][x] = [b]
 * [ aI ]      [0]
 * Because of the regularization, we assure the final matrix is full rank
 *
 * NOTE: the size of A should be at least (m+n)*n, and b should be at
 *       least (m + n)
 */
template <typename T>
class TikhonovQrSolver
{
    public:
        typedef T        value_type;
        typedef typename carbine::TrivialType<T>::type trivial_type;

        TikhonovQrSolver():m_m(0), m_n(0), m_lwork(0),
                m_eps((trivial_type)1E-8)
        {
            m_workspace.resize(2);
        }
        void set_size(int m, int n, int nrhs = 1);
        trivial_type& epsilon()    { return m_eps; }

        void solve(T* A, T* b);
        void solve(T* A, T* b, FILE* logfd, int tid);

    private:
        int         m_m;
        int         m_n;
        int         m_nrow;
        int         m_nrhs;
        int         m_lwork;
        trivial_type m_eps;

        std::vector<value_type>     m_workspace;
};

template<>
void TikhonovQrSolver< std::complex<double> >::
set_size(int m, int n, int nrhs)
{
    m_m = m + n;
    m_n = n;
    m_nrhs = nrhs;
    m_nrow = m;

    // call the routine in query mode to get the optimal size of workspace
    std::complex<double> *a = NULL;
    std::complex<double> *b = NULL;

    int MINUS_ONE = -1;
    int ret;
    char TRANS = 'N';

    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)a, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &MINUS_ONE,
          &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgels (query mode) %d\n", ret);
        exit(1);
    }

    // allocate workspace
    m_workspace.resize((int)m_workspace[0].real());
    m_lwork = m_workspace.size();
}

void TikhonovQrSolver< std::complex<double> >::
solve(std::complex<double>* A,
      std::complex<double>* b)
{
    using namespace std; 

    //// estimate the maximum singular value by forbenius norm
    double nrm = 0;
    complex<double> tv;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += m_m)
    {
        cblas_zdotc_sub(m_nrow, &A[i*m_m], 1, &A[i*m_m], 1, &tv);
        nrm += tv.real();
        memset(&A[rr], 0, sizeof(complex<double>)*m_n);
    }
    //// construct the new left-hand matrix A
    nrm = sqrt(nrm) * m_eps;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += (m_m+1))
        A[rr].real() = nrm;
    memset(&b[m_nrow], 0, sizeof(complex<double>)*m_n);

    int ret;
    char TRANS = 'N';
    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)A, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &m_lwork,
          &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgels (query mode) %d\n", ret);
        exit(1);
    }

    //MKL_FreeBuffers();  /* free memory */
}

void TikhonovQrSolver< std::complex<double> >::
solve(std::complex<double>* A,
      std::complex<double>* b,
      FILE* logfd, int tid)
{
    using namespace std; 

    fprintf(logfd, ">>>> MSG: solve 0 ... thread #%d\n", tid);
    fflush(logfd);

    //// estimate the maximum singular value by forbenius norm
    double nrm = 0;
    complex<double> tv;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += m_m)
    {
        cblas_zdotc_sub(m_nrow, &A[i*m_m], 1, &A[i*m_m], 1, &tv);
        nrm += tv.real();
        memset(&A[rr], 0, sizeof(complex<double>)*m_n);
    }
    //// construct the new left-hand matrix A
    nrm = sqrt(nrm) * m_eps;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += (m_m+1))
        A[rr].real() = nrm;
    memset(&b[m_nrow], 0, sizeof(complex<double>)*m_n);

    fprintf(logfd, ">>>> MSG: solve 1 ... thread #%d\n", tid);
    fflush(logfd);

    int ret;
    char TRANS = 'N';
    zgels(&TRANS, &m_m, &m_n, &m_nrhs,
          (MKL_Complex16*)A, &m_m,
          (MKL_Complex16*)b, &m_m,
          (MKL_Complex16*)&m_workspace[0], &m_lwork,
          &ret);

    fprintf(logfd, ">>>> MSG: solve 2 ... thread #%d\n", tid);
    fflush(logfd);

    if ( ret )
    {
        fprintf(logfd, "ERROR: on zgels (query mode) %d\n", ret);
        exit(1);
    }

    MKL_FreeBuffers();  /* free memory */

//int PTR;
//int vv = MKL_MemStat(&PTR);
//fprintf(logfd, ">>>> MSG: %ld bytes are allocated in %d buffers by MKL thread #%d\n", (long)vv, PTR, tid);
//fflush(logfd);
}
///////////////////////////////////////////////////////////////////////////////
/*!
 * Solve least-square using Tikhonov regularization and SVD. It is pretty similar
 * to the TSVD method, except that it constructs the least-square problem as
 * [  A ][x] = [b]
 * [ aI ]      [0]
 *
 * NOTE: the size of A should be at least (m+n)*n, and b should be at
 *       least (m + n)
 * Matrice are in column-major following the fortran style.
 */
template <typename T>
class TikhonovSvdSolver
{
    public:
        typedef T        value_type;
        typedef typename carbine::TrivialType<T>::type trivial_type;

        TikhonovSvdSolver():m_eps((trivial_type)1E-8),
                m_m(0), m_n(0), m_truncation(1E-8), 
                m_lwork(0), m_iwork(0), m_rwork(0)
        {
            m_lworkSpace.resize(2);
            m_rworkSpace.resize(2);
            m_iworkSpace.resize(2);
        }

        void solve(T* A, T* b);
        int  m() const { return m_m - m_n; }
        int  n() const { return m_n; }
        void set_size(int m, int n, int nrhs = 1)
        {
                m_m = m + n;
                m_n = n;
                m_nrow = m;
                m_nrhs = nrhs;
        }

        trivial_type& truncation() { return m_truncation; }
        trivial_type& epsilon()    { return m_eps; }

        const std::vector<trivial_type>& singular_values() const
        {   return m_singulars; }

        trivial_type min_singular_value() const;

        int rank() const { return m_rank; }

    private:
        int         m_m;
        int         m_n;
        int         m_nrow;
        int         m_nrhs;
        trivial_type m_truncation;
        trivial_type m_eps;
        int         m_rank;
        int         m_lwork, m_iwork, m_rwork;      // size of workspace for ?gelsd routine

        std::vector<value_type>     m_lworkSpace;   // workspace given to ?gelsd routine
        std::vector<trivial_type>    m_rworkSpace;
        std::vector<int>            m_iworkSpace;
        std::vector<trivial_type>    m_singulars;
};

template<>
double TikhonovSvdSolver< std::complex<double> >::
min_singular_value() const
{
    return m_singulars[m_rank-1];
}

/*!
 * \param m number of rows for the original problem
 * \param n number of columns for the original problem
 */
template<>
void TikhonovSvdSolver< std::complex<double> >::
set_size(int m, int n, int nrhs)
{
    m_m = m + n;
    m_n = n;
    m_nrhs = nrhs;
    m_nrow = m;

    // call the routine in query mode to get the optimal size of workspace
    std::complex<double> *a = NULL;
    std::complex<double> *b = NULL;
    int MINUS_ONE = -1;
    int ret;

    zgelsd(&m_m, &m_n, &m_nrhs, 
           (MKL_Complex16*)a, &m_m, 
           (MKL_Complex16*)b, &m_m, 
           &m_singulars[0], &m_truncation, &m_rank, 
           (MKL_Complex16*)&m_lworkSpace[0], &MINUS_ONE, 
           &m_rworkSpace[0], &m_iworkSpace[0], &ret);

    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgelsd (query mode) %d\n", ret);
        exit(1);
    }

    // allocate memory for workspace
    m_lworkSpace.resize((int)m_lworkSpace[0].real());
    m_lwork = m_lworkSpace.size();

    const int S = std::min(m, n);
    int dummy = -1;
    int ispec = 9;
    int smlsiz = ilaenv(&ispec, "zgelsd", "N", &m_m, &m_n, &m_nrhs, &dummy);
    int nlvl   = S / (smlsiz + 1) + 1;  // approximate nlvl without log_2
    m_iwork = 3*S*nlvl + 11*S;
    m_iworkSpace.resize(m_iwork);

    m_rwork = 10*S + 2*S*smlsiz + 8*n*nlvl + 3*smlsiz*m_nrhs + (smlsiz+1)*(smlsiz+1);
    m_rworkSpace.resize(m_rwork);

    m_singulars.resize(S);
}

/*!
 * NOTE: input data should be fortran-stype, i.e. data stored colume-major
 * rather than row-major
 */
template<>
void TikhonovSvdSolver< std::complex<double> >::
solve(std::complex<double>* A, 
      std::complex<double>* b)
{
    using namespace std;

    //// estimate the maximum singular value by forbenius norm
    double nrm = 0;
    complex<double> tv;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += m_m)
    {
        cblas_zdotc_sub(m_nrow, &A[i*m_m], 1, &A[i*m_m], 1, &tv);
        nrm += tv.real();
        memset(&A[rr], 0, sizeof(complex<double>)*m_n);
    }
    //// construct the new left-hand matrix A
    nrm = sqrt(nrm) * m_eps;
    for(int i = 0, rr = m_nrow;i < m_n;++ i, rr += (m_m+1))
        A[rr].real() = nrm;
    memset(&b[m_nrow], 0, sizeof(complex<double>)*m_n);


    int ret;
    zgelsd(&m_m, &m_n, &m_nrhs,
           (MKL_Complex16*)A, &m_m,
           (MKL_Complex16*)b, &m_m,
           &m_singulars[0], &m_truncation, &m_rank,
           (MKL_Complex16*)&m_lworkSpace[0], &m_lwork,
           &m_rworkSpace[0], &m_iworkSpace[0], &ret);

    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgelsd (query mode) %d\n", ret);
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
/*!
 * Truncated SVD solver is used for rank deficiency case. 
 * It calls LAPACK routine, ?gelsd, using SVD and divide-and-conquer method
 * to solve the problem.
 *
 * It has an extra method, truncation(), to get/set truncation, which is 
 * used as the parameter, rcond, in ?gelsd LAPACK routine.
 *
 * NOTE: this solver is not thread-safe, you have to create an instance of
 * this solver for each thread.
 * Matrice are in column-major following the fortran style.
 */
template<typename T>
class TruncatedSvdSolver
{
    public:
        typedef T        value_type;
        typedef typename carbine::TrivialType<T>::type trivial_type;

        TruncatedSvdSolver():
                m_m(0), m_n(0), m_truncation(1E-8), 
                m_lwork(0), m_iwork(0), m_rwork(0)
        { 
            m_lworkSpace.resize(2);
            m_rworkSpace.resize(2);
            m_iworkSpace.resize(2);
        }

        /*! both A and b are overwritten at output */
        void solve(T* A, T* b);

        int m() const { return m_m; }
        int n() const { return m_n; }
        void set_size(int m, int n, int nrhs = 1)
        {
            m_m = m;
            m_n = n;
            m_nrhs = nrhs;
        }

        trivial_type& truncation()
        {   return m_truncation; }

        const std::vector<trivial_type>& singular_values()
        {   return m_singulars; }

        int rank() const { return m_rank; }

    private:
        int         m_m;
        int         m_n;
        int         m_nrhs;
        trivial_type m_truncation;
        int         m_rank;
        int         m_lwork, m_iwork, m_rwork;      // size of workspace for ?gelsd routine

        std::vector<value_type>     m_lworkSpace;   // workspace given to ?gelsd routine
        std::vector<trivial_type>    m_rworkSpace;
        std::vector<int>            m_iworkSpace;
        std::vector<trivial_type>    m_singulars;
};

template<>
void TruncatedSvdSolver< std::complex<double> >::
set_size(int m, int n, int nrhs)
{
    m_m = m;
    m_n = n;
    m_nrhs = nrhs;

    // call the routine in query mode to get the optimal size of workspace
    std::complex<double> *a = NULL;
    std::complex<double> *b = NULL;
    int MINUS_ONE = -1;
    int ret;

    zgelsd(&m_m, &m_n, &m_nrhs, 
           (MKL_Complex16*)a, &m_m, 
           (MKL_Complex16*)b, &m_m, 
           &m_singulars[0], &m_truncation, &m_rank, 
           (MKL_Complex16*)&m_lworkSpace[0], &MINUS_ONE, 
           &m_rworkSpace[0], &m_iworkSpace[0], &ret);

    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgelsd (query mode) %d\n", ret);
        exit(1);
    }

    // allocate memory for workspace
    m_lworkSpace.resize((int)m_lworkSpace[0].real());
    m_lwork = m_lworkSpace.size();

    const int S = std::min(m, n);
    int dummy = -1;
    int ispec = 9;
    int smlsiz = ilaenv(&ispec, "zgelsd", "N", &m_m, &m_n, &m_nrhs, &dummy);
    int nlvl   = S / (smlsiz + 1) + 1;  // approximate nlvl without log_2
    m_iwork = 3*S*nlvl + 11*S;
    m_iworkSpace.resize(m_iwork);

    m_rwork = 10*S + 2*S*smlsiz + 8*n*nlvl + 3*smlsiz*m_nrhs + (smlsiz+1)*(smlsiz+1);
    m_rworkSpace.resize(m_rwork);

    m_singulars.resize(S);
}

/*!
 * NOTE: input data should be fortran-stype, i.e. data stored colume-major
 * rather than row-major
 */
template<>
void TruncatedSvdSolver< std::complex<double> >::
solve(std::complex<double>* A, 
      std::complex<double>* b)
{
    int ret;
    zgelsd(&m_m, &m_n, &m_nrhs,
           (MKL_Complex16*)A, &m_m,
           (MKL_Complex16*)b, &m_m,
           &m_singulars[0], &m_truncation, &m_rank,
           (MKL_Complex16*)&m_lworkSpace[0], &m_lwork,
           &m_rworkSpace[0], &m_iworkSpace[0], &ret);

    if ( ret )
    {
        fprintf(stderr, "ERROR: on zgelsd (query mode) %d\n", ret);
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
template<typename T>
typename carbine::TrivialType<T>::type compute_fitting_error(
        const int order, 
        const T *U, const T *b, const T *x, 
        const int nrow, const int width);

template<>
double compute_fitting_error< std::complex<double> >(
        const int order,
        const std::complex<double> *A, 
        const std::complex<double> *b,
        const std::complex<double> *x, 
        const int nrow, const int width)
{
    using namespace std;

    const complex<double> beta(0, 0);
    complex<double>       alpha(1, 0);
    complex<double>       residual[nrow];

    switch (order)
    {
        case CblasRowMajor:
            cblas_zgemv(CblasRowMajor, CblasNoTrans, nrow, width, &alpha, 
                        A, width, x, 1, &beta, residual, 1);
            break;
        case CblasColMajor:
            cblas_zgemv(CblasColMajor, CblasNoTrans, nrow, width, &alpha, 
                        A, nrow, x, 1, &beta, residual, 1);
            break;
        default:
            fprintf(stderr, "Unknown matrix layout order!\n");
            exit(1);
    }

    alpha.real() = -1;
    cblas_zaxpy(nrow, &alpha, b, 1, residual, 1);

    return cblas_dznrm2(nrow, residual, 1) / cblas_dznrm2(nrow, b, 1);
}

template<typename T>
typename carbine::TrivialType<T>::type compute_residual(
        const int order, 
        const T *U, T *b, const T *x, 
        const int nrow, const int width);

/*
 * NOTE: after returning from this method, the vector b is overwritten by 
 *       the residual. i.e. b = Ax - b
 */
template<>
double compute_residual< std::complex<double> >(
        const int order,
        const std::complex<double> *A, 
              std::complex<double> *b,
        const std::complex<double> *x, 
        const int nrow, const int width)
{
    using namespace std;

    const complex<double> beta(-1, 0);
    complex<double>       alpha(1, 0);

    double nrm_b = cblas_dznrm2(nrow, b, 1);

    switch (order)
    {
        case CblasRowMajor:
            cblas_zgemv(CblasRowMajor, CblasNoTrans, nrow, width, &alpha, 
                        A, width, x, 1, &beta, b, 1);
            break;
        case CblasColMajor:
            cblas_zgemv(CblasColMajor, CblasNoTrans, nrow, width, &alpha, 
                        A, nrow, x, 1, &beta, b, 1);
            break;
        default:
            fprintf(stderr, "Unknown matrix layout order!\n");
            exit(1);
    }

    return cblas_dznrm2(nrow, b, 1) / nrm_b;
}

#endif
