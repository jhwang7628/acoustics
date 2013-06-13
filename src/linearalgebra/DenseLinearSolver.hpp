#ifndef DENSE_LINEAR_SOLVER
#   define DENSE_LINEAR_SOLVER

#include <vector>
#include <mkl.h>
#include <mkl_lapack.h>

template <typename T>
class SymmLinearSolver
{
    public:
        typedef T       value_type;

        SymmLinearSolver():n_(0), lwork_(0)
        {   workspace_.resize(2); }

        void set_size(int n, char uplow='L', int nrhs = 1);
        void solve(T* A, T* b);

    private:
        int     n_;
        int     lwork_;
        int     nrhs_;
        char    uplo_;
        std::vector<int>            ipiv_;
        std::vector<value_type>     workspace_;
};

template<>
void SymmLinearSolver<double>::set_size(int n, char uplo, int nrhs)
{
    n_ = n;
    nrhs_ = nrhs;
    ipiv_.resize(n);
    uplo_ = uplo;

    double *a = NULL, *b = NULL;
    int MINUS_ONE = -1;
    int ret;

    dsysv(&uplo, &n_, &nrhs_, a, &n_, &ipiv_[0], b, &n_, 
          &workspace_[0], &MINUS_ONE, &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on dsysv(query mode) %d\n", ret);
        exit(1);
    }

    lwork_ = (int)workspace_[0];
    workspace_.resize(lwork_);
}

template<>
void SymmLinearSolver<double>::solve(double* A, double* b)
{
    int ret;
    dsysv(&uplo_, &n_, &nrhs_, A, &n_, &ipiv_[0], b, &n_,
          &workspace_[0], &lwork_, &ret);
    if ( ret )
    {
        fprintf(stderr, "ERROR: on dsysv(solve mode) %d\n", ret);
        exit(1);
    }
}

#endif
