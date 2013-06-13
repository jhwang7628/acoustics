#include <gtest/gtest.h>
#include "linearalgebra/LeastSquareSolver.hpp"
#include "linearalgebra/DenseLinearSolver.hpp"

using namespace std;

TEST(TestLsQrSolver, test_double)
{
    double A[] = {9, 8, 12, 14, 15, 5, 13, 13, 3, 2, 9, 19, 6, 11, 7};
    double b[] = {2, 7, 7, 3, 5};
    LsQrSolver<double> sol;
    sol.set_size(5,3);
    sol.solve(A, b);
    for(int i = 0;i < 3;++ i)
        cout << b[i] << endl;
}

TEST(TestLsQRSolver, test_linearsolver)
{
    double A[] = {
           0.814723686393179,  0.905791937075619,  0.126986816293506,  0.913375856139019,
           0.632359246225410,   0.097540404999410,   0.278498218867048,   0.546881519204984,
           0.957506835434298,   0.964888535199277,   0.157613081677548,   0.970592781760616,
           0.957166948242946,   0.485375648722841,   0.800280468888800,   0.141886338627215 };
    double b[] = {   0.421761282626275, 0.915735525189067, 0.792207329559554, 0.959492426392903 };
    int ipiv[4];
    int N = 4, nrhs = 1;
    int info;
    dgesv(&N, &nrhs, A, &N, ipiv, b, &N, &info);
    if ( info != 0 )
        cerr << "ERROR: " << info << endl;
    else
        for(int i = 0;i < 4;++ i)
            cout << i << " :  " << b[i] << endl;
}

TEST(TestLsQRSolver, test_symm_solver)
{
    double A[] = { 16, 15, 10, 18, 0, 1, 11, 9, 0, 0, 2, 17, 0, 0, 0, 2 };
    double b[] = {1, 2, 3, 4};
    SymmLinearSolver<double> sol;
    sol.set_size(4);
    sol.solve(A, b);
    for(int i = 0;i < 4;++ i)
        cout << b[i] << endl;
}

TEST(TestLsQRSolver, test_symm_solver_packed)
{
    double A[] = { 16, 15, 10, 18, 1, 11, 9, 2, 17, 2 };
    double b[] = {1, 2, 3, 4};

    char uplo = 'L';
    int N = 4, nrhs = 1;
    int ipiv[4];
    int info;
    dspsv(&uplo, &N, &nrhs, A, ipiv, b, &N, &info);
    if ( info != 0 )
        cerr << "ERROR: " << info << endl;
    else
        for(int i = 0;i < 4;++ i)
            cout << i << " :  " << b[i] << endl;
}
