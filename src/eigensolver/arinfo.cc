#include <stdio.h>

// Auto-generated from ARPACK sources

static void dnaupd_info_string(int info)
{
    if (info == 0) {
        printf("Normal exit.\n");
    }
    if (info == 1) {
        printf("Maximum number of iterations taken.\n");
        printf("All possible eigenvalues of OP has been found. IPARAM(5)\n");
        printf("returns the number of wanted converged Ritz values.\n");
    }
    if (info == 2) {
        printf("No longer an informational error. Deprecated starting\n");
        printf("with release 2 of ARPACK.\n");
    }
    if (info == 3) {
        printf("No shifts could be applied during a cycle of the\n");
        printf("Implicitly restarted Arnoldi iteration. One possibility\n");
        printf("is to increase the size of NCV relative to NEV.\n");
        printf("See remark 4 below.\n");
    }
    if (info == -1) {
        printf("N must be positive.\n");
    }
    if (info == -2) {
        printf("NEV must be positive.\n");
    }
    if (info == -3) {
        printf("NCV-NEV >= 2 and less than or equal to N.\n");
    }
    if (info == -4) {
        printf("The maximum number of Arnoldi update iteration\n");
        printf("must be greater than zero.\n");
    }
    if (info == -5) {
        printf("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
    }
    if (info == -6) {
        printf("BMAT must be one of 'I' or 'G'.\n");
    }
    if (info == -7) {
        printf("Length of private work array is not sufficient.\n");
    }
    if (info == -8) {
        printf("Error return from LAPACK eigenvalue calculation;\n");
    }
    if (info == -9) {
        printf("Starting vector is zero.\n");
    }
    if (info == -10) {
        printf("IPARAM(7) must be 1,2,3,4.\n");
    }
    if (info == -11) {
        printf("IPARAM(7) = 1 and BMAT = 'G' are incompatable.\n");
    }
    if (info == -12) {
        printf("IPARAM(1) must be equal to 0 or 1.\n");
    }
    if (info == -9999) {
        printf("Could not build an Arnoldi factorization.\n");
        printf("IPARAM(5) returns the size of the current Arnoldi\n");
        printf("factorization.\n");
    }
}
static void dneupd_info_string(int info)
{
    if (info == 0) {
        printf("Normal exit.\n");
    }
    if (info == 1) {
        printf("The Schur form computed by LAPACK routine dlahqr\n");
        printf("could not be reordered by LAPACK routine dtrsen .\n");
        printf("Re-enter subroutine dneupd with IPARAM(5)=NCV and\n");
        printf("increase the size of the arrays DR and DI to have\n");
        printf("dimension at least dimension NCV and allocate at least NCV\n");
        printf("columns for Z. NOTE: Not necessary if Z and V share\n");
        printf("the same space. Please notify the authors if this error\n");
        printf("occurs.\n");
    }
    if (info == -1) {
        printf("N must be positive.\n");
    }
    if (info == -2) {
        printf("NEV must be positive.\n");
    }
    if (info == -3) {
        printf("NCV-NEV >= 2 and less than or equal to N.\n");
    }
    if (info == -5) {
        printf("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
    }
    if (info == -6) {
        printf("BMAT must be one of 'I' or 'G'.\n");
    }
    if (info == -7) {
        printf("Length of private work WORKL array is not sufficient.\n");
    }
    if (info == -8) {
        printf("Error return from calculation of a real Schur form.\n");
        printf("Informational error from LAPACK routine dlahqr .\n");
    }
    if (info == -9) {
        printf("Error return from calculation of eigenvectors.\n");
        printf("Informational error from LAPACK routine dtrevc .\n");
    }
    if (info == -10) {
        printf("IPARAM(7) must be 1,2,3,4.\n");
    }
    if (info == -11) {
        printf("IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
    }
    if (info == -12) {
        printf("HOWMNY = 'S' not yet implemented\n");
    }
    if (info == -13) {
        printf("HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
    }
    if (info == -14) {
        printf("DNAUPD did not find any eigenvalues to sufficient\n");
        printf("accuracy.\n");
    }
    if (info == -15) {
        printf("DNEUPD got a different count of the number of converged\n");
        printf("Ritz values than DNAUPD got. This indicates the user\n");
        printf("probably made an error in passing data from DNAUPD to\n");
        printf("DNEUPD or that the data was modified before entering\n");
        printf("DNEUPD\n");
    }
}
static void dsaupd_info_string(int info)
{
    if (info == 0) {
        printf("Normal exit.\n");
    }
    if (info == 1) {
        printf("Maximum number of iterations taken.\n");
        printf("All possible eigenvalues of OP has been found. IPARAM(5)\n");
        printf("returns the number of wanted converged Ritz values.\n");
    }
    if (info == 2) {
        printf("No longer an informational error. Deprecated starting\n");
        printf("with release 2 of ARPACK.\n");
    }
    if (info == 3) {
        printf("No shifts could be applied during a cycle of the\n");
        printf("Implicitly restarted Arnoldi iteration. One possibility\n");
        printf("is to increase the size of NCV relative to NEV.\n");
        printf("See remark 4 below.\n");
    }
    if (info == -1) {
        printf("N must be positive.\n");
    }
    if (info == -2) {
        printf("NEV must be positive.\n");
    }
    if (info == -3) {
        printf("NCV must be greater than NEV and less than or equal to N.\n");
    }
    if (info == -4) {
        printf("The maximum number of Arnoldi update iterations allowed\n");
        printf("must be greater than zero.\n");
    }
    if (info == -5) {
        printf("WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.\n");
    }
    if (info == -6) {
        printf("BMAT must be one of 'I' or 'G'.\n");
    }
    if (info == -7) {
        printf("Length of private work array WORKL is not sufficient.\n");
    }
    if (info == -8) {
        printf("Error return from trid. eigenvalue calculation;\n");
        printf("Informatinal error from LAPACK routine dsteqr .\n");
    }
    if (info == -9) {
        printf("Starting vector is zero.\n");
    }
    if (info == -10) {
        printf("IPARAM(7) must be 1,2,3,4,5.\n");
    }
    if (info == -11) {
        printf("IPARAM(7) = 1 and BMAT = 'G' are incompatable.\n");
    }
    if (info == -12) {
        printf("IPARAM(1) must be equal to 0 or 1.\n");
    }
    if (info == -13) {
        printf("NEV and WHICH = 'BE' are incompatable.\n");
    }
    if (info == -9999) {
        printf("Could not build an Arnoldi factorization.\n");
        printf("IPARAM(5) returns the size of the current Arnoldi\n");
        printf("factorization. The user is advised to check that\n");
        printf("enough workspace and array storage has been allocated.\n");
    }
}
static void dseupd_info_string(int info)
{
    if (info == 0) {
        printf("Normal exit.\n");
    }
    if (info == -1) {
        printf("N must be positive.\n");
    }
    if (info == -2) {
        printf("NEV must be positive.\n");
    }
    if (info == -3) {
        printf("NCV must be greater than NEV and less than or equal to N.\n");
    }
    if (info == -5) {
        printf("WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.\n");
    }
    if (info == -6) {
        printf("BMAT must be one of 'I' or 'G'.\n");
    }
    if (info == -7) {
        printf("Length of private work WORKL array is not sufficient.\n");
    }
    if (info == -8) {
        printf("Error return from trid. eigenvalue calculation;\n");
        printf("Information error from LAPACK routine dsteqr .\n");
    }
    if (info == -9) {
        printf("Starting vector is zero.\n");
    }
    if (info == -10) {
        printf("IPARAM(7) must be 1,2,3,4,5.\n");
    }
    if (info == -11) {
        printf("IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
    }
    if (info == -12) {
        printf("NEV and WHICH = 'BE' are incompatible.\n");
    }
    if (info == -14) {
        printf("DSAUPD did not find any eigenvalues to sufficient\n");
        printf("accuracy.\n");
    }
    if (info == -15) {
        printf("HOWMNY must be one of 'A' or 'S' if RVEC = .true.\n");
    }
    if (info == -16) {
        printf("HOWMNY = 'S' not yet implemented\n");
    }
    if (info == -17) {
        printf("DSEUPD got a different count of the number of converged\n");
        printf("Ritz values than DSAUPD got. This indicates the user\n");
        printf("probably made an error in passing data from DSAUPD to\n");
        printf("DSEUPD or that the data was modified before entering\n");
        printf("DSEUPD .\n");
    }
}
static void znaupd_info_string(int info)
{
    if (info == 0) {
        printf("Normal exit.\n");
    }
    if (info == 1) {
        printf("Maximum number of iterations taken.\n");
        printf("All possible eigenvalues of OP has been found. IPARAM(5)\n");
        printf("returns the number of wanted converged Ritz values.\n");
    }
    if (info == 2) {
        printf("No longer an informational error. Deprecated starting\n");
        printf("with release 2 of ARPACK.\n");
    }
    if (info == 3) {
        printf("No shifts could be applied during a cycle of the\n");
        printf("Implicitly restarted Arnoldi iteration. One possibility\n");
        printf("is to increase the size of NCV relative to NEV.\n");
        printf("See remark 4 below.\n");
    }
    if (info == -1) {
        printf("N must be positive.\n");
    }
    if (info == -2) {
        printf("NEV must be positive.\n");
    }
    if (info == -3) {
        printf("NCV-NEV >= 2 and less than or equal to N.\n");
    }
    if (info == -4) {
        printf("The maximum number of Arnoldi update iteration\n");
        printf("must be greater than zero.\n");
    }
    if (info == -5) {
        printf("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
    }
    if (info == -6) {
        printf("BMAT must be one of 'I' or 'G'.\n");
    }
    if (info == -7) {
        printf("Length of private work array is not sufficient.\n");
    }
    if (info == -8) {
        printf("Error return from LAPACK eigenvalue calculation;\n");
    }
    if (info == -9) {
        printf("Starting vector is zero.\n");
    }
    if (info == -10) {
        printf("IPARAM(7) must be 1,2,3.\n");
    }
    if (info == -11) {
        printf("IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
    }
    if (info == -12) {
        printf("IPARAM(1) must be equal to 0 or 1.\n");
    }
    if (info == -9999) {
        printf("Could not build an Arnoldi factorization.\n");
        printf("User input error highly likely. Please\n");
        printf("check actual array dimensions and layout.\n");
        printf("IPARAM(5) returns the size of the current Arnoldi\n");
        printf("factorization.\n");
    }
}
static void zneupd_info_string(int info)
{
    if (info == 0) {
        printf("Normal exit.\n");
    }
    if (info == 1) {
        printf("The Schur form computed by LAPACK routine csheqr\n");
        printf("could not be reordered by LAPACK routine ztrsen .\n");
        printf("Re-enter subroutine zneupd with IPARAM(5)=NCV and\n");
        printf("increase the size of the array D to have\n");
        printf("dimension at least dimension NCV and allocate at least NCV\n");
        printf("columns for Z. NOTE: Not necessary if Z and V share\n");
        printf("the same space. Please notify the authors if this error\n");
        printf("occurs.\n");
    }
    if (info == -1) {
        printf("N must be positive.\n");
    }
    if (info == -2) {
        printf("NEV must be positive.\n");
    }
    if (info == -3) {
        printf("NCV-NEV >= 2 and less than or equal to N.\n");
    }
    if (info == -5) {
        printf("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
    }
    if (info == -6) {
        printf("BMAT must be one of 'I' or 'G'.\n");
    }
    if (info == -7) {
        printf("Length of private work WORKL array is not sufficient.\n");
    }
    if (info == -8) {
        printf("Error return from LAPACK eigenvalue calculation.\n");
        printf("This should never happened.\n");
    }
    if (info == -9) {
        printf("Error return from calculation of eigenvectors.\n");
        printf("Informational error from LAPACK routine ztrevc .\n");
    }
    if (info == -10) {
        printf("IPARAM(7) must be 1,2,3\n");
    }
    if (info == -11) {
        printf("IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
    }
    if (info == -12) {
        printf("HOWMNY = 'S' not yet implemented\n");
    }
    if (info == -13) {
        printf("HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
    }
    if (info == -14) {
        printf("ZNAUPD did not find any eigenvalues to sufficient\n");
        printf("accuracy.\n");
    }
    if (info == -15) {
        printf("ZNEUPD got a different count of the number of converged\n");
        printf("Ritz values than ZNAUPD got. This indicates the user\n");
        printf("probably made an error in passing data from ZNAUPD to\n");
        printf("ZNEUPD or that the data was modified before entering\n");
        printf("ZNEUPD\n");
    }
}
