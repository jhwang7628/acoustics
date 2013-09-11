/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>
#include <vector>
#include <algorithm>

#include <complex>

#include "areigs.h"

typedef std::complex<double> dcomplex;
using std::vector;
using std::copy;
using std::fill;

typedef int fortran_int;

extern "C" {

int dsaupd_(fortran_int*        IDO,
            const char* BMAT,
            fortran_int*        N,
            const char* WHICH,
            fortran_int*        NEV,
            double*     TOL,
            double*     RESID,
            fortran_int*        NCV,
            double*     V,
            fortran_int*        LDV,
            fortran_int*        IPARAM,
            fortran_int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            fortran_int*        LWORKL,
            fortran_int*        INFO );

int dseupd_(fortran_int*        RVEC,
            const char* HOWMNY,
            fortran_int*        SELECT,
            double*     D,
            double*     Z,
            fortran_int*        LDZ,
            double*     SIGMA,
            const char* BMAT,
            fortran_int*        N,
            const char* WHICH,
            fortran_int*        NEV,
            double*     TOL,
            double*     RESID,
            fortran_int*        NCV,
            double*     V,
            fortran_int*        LDV,
            fortran_int*        IPARAM,
            fortran_int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            fortran_int*        LWORKL,
            fortran_int*        INFO );

int dnaupd_(fortran_int*        IDO,
            const char* BMAT,
            fortran_int*        N,
            const char* WHICH,
            fortran_int*        NEV,
            double*     TOL,
            double*     RESID,
            fortran_int*        NCV,
            double*     V,
            fortran_int*        LDV,
            fortran_int*        IPARAM,
            fortran_int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            fortran_int*        LWORKL,
            fortran_int*        INFO );

int dneupd_(fortran_int*        RVEC,
            const char* HOWMNY,
            fortran_int*        SELECT,
            double*     DR,
            double*     DI,
            double*     Z,
            fortran_int*        LDZ,
            double*     SIGMAR,
            double*     SIGMAI,
            double*     WORKEV,
            const char* BMAT,
            fortran_int*        N,
            const char* WHICH,
            fortran_int*        NEV,
            double*     TOL,
            double*     RESID,
            fortran_int*        NCV,
            double*     V,
            fortran_int*        LDV,
            fortran_int*        IPARAM,
            fortran_int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            fortran_int*        LWORKL,
            fortran_int*        INFO );

int znaupd_(fortran_int*        IDO,
            const char* BMAT,
            fortran_int*        N,
            const char* WHICH,
            fortran_int*        NEV,
            double*     TOL,
            dcomplex*   RESID,
            fortran_int*        NCV,
            dcomplex*   V,
            fortran_int*        LDV,
            fortran_int*        IPARAM,
            fortran_int*        IPNTR,
            dcomplex*   WORKD,
            dcomplex*   WORKL,
            fortran_int*        LWORKL,
            double*     RWORK,
            fortran_int*        INFO );

int zneupd_(fortran_int*        RVEC,
            const char* HOWMNY,
            fortran_int*        SELECT,
            dcomplex*   D,
            dcomplex*   Z,
            fortran_int*        LDZ,
            dcomplex*   SIGMA,
            dcomplex*   WORKEV,
            const char* BMAT,
            fortran_int*        N,
            const char* WHICH,
            fortran_int*        NEV,
            double*     TOL,
            dcomplex*   RESID,
            fortran_int*        NCV,
            dcomplex*   V,
            fortran_int*        LDV,
            fortran_int*        IPARAM,
            fortran_int*        IPNTR,
            dcomplex*   WORKD,
            dcomplex*   WORKL,
            fortran_int*        LWORKL,
            double*     RWORK,
            fortran_int*        INFO );

}


#include "arinfo.cc"


// ---- General settings ----


Arpack::Arpack()
{
    n = 0;
    nev = 6;
    ncv = 20;
    rvec = 1;
    mode = 1;
    maxitr = 30;

    bmat[0] = 'I';
    bmat[1] = '\0';

    which[0] = 'L';
    which[1] = 'M';
    which[2] = '\0';

    tol   = 0.0;
}


Arpack::~Arpack()
{
}


void Arpack::set_n(int n_)
{
    n = n_;
}


void Arpack::set_nev(int nev_)
{
    nev = nev_;
    ncv = 2*nev_;
    if (ncv < 20)
        ncv = 20;
}


void Arpack::set_ncv(int ncv_)
{
    ncv = ncv_;
    if (ncv < nev)
        ncv = nev;
}


void Arpack::set_rvec(int rvec_)
{
    rvec = rvec_;
}


void Arpack::set_maxitr(int maxitr_)
{
    maxitr = maxitr_;
}


void Arpack::set_which(const char* which_)
{
    which[0] = which_[0];
    which[1] = which_[1];
}


void Arpack::set_tol(double tol_)
{
    tol = tol_;
}


// ---- DS settings ----


ArpackDS::ArpackDS()
{
    set_mode(1);
    set_shift(0);
}


ArpackDS::~ArpackDS()
{
}


void ArpackDS::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else {
        bmat[0] = 'G';
    }
}


void ArpackDS::set_shift(double sigma_)
{
    sigma = sigma_;
}


int ArpackDS::compute_eigs(double* d, double* v)
{
    vector<double> resid(n);
    vector<double> workd(3*n);
    vector<double> workl((ncv+8)*ncv);
    vector<fortran_int>    select(ncv);

    fortran_int lworkl = (ncv+8)*ncv;
    fortran_int iparam[11];
    fortran_int ipntr[14];
    fortran_int info   = 0;
    fortran_int ierr   = 0;
    fortran_int status = 0;
    fortran_int iter   = 0;
    fortran_int ido    = 0;

    // -- Clear temporaries
    fill(iparam, iparam+11, 0);
    fill(ipntr,  ipntr+14,  0);

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    ido = 0;
    while (1) {

        dsaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[0],
                &ncv, v, &n, iparam, ipntr, &workd[0], 
                &workl[0], &lworkl, &info);

        if (ido == -1) {
            times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 1) {
            if (mode == 3 || mode == 4 || mode == 5) {
                times_OP2(&workd[ipntr[0]-1], &workd[ipntr[1]-1],
                          &workd[ipntr[2]-1]);
            } else {
                times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
            }
        } else if (ido == 2) {
            times_M(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else {
            break;
        }

        ++iter;
    }

    //info = 0;
    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in dsaupd_ = %d\n", info);
        dsaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        dseupd_(&rvec, "A", &select[0], d, v, &n, &sigma,
                bmat, &n, which, &nev, &tol, &resid[0], &ncv,
                v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl, &ierr);

        if (ierr != 0) {
            printf("Error in dseupd_ = %d\n", ierr);
            dseupd_info_string(ierr);
            status = -2;
        }
    }

    return status;
}


void ArpackDS::times_OP1(double*, double*)
{
}


void ArpackDS::times_OP2(double* x, double* opx, double*)
{
    times_OP1(x, opx);
}


void ArpackDS::times_M(double* x, double* Mx)
{
    copy(Mx, Mx+n, x);
}


// ---- DN settings ----


ArpackDN::ArpackDN()
{
    set_mode(1);
    set_shift(0);
}


ArpackDN::~ArpackDN()
{
}


void ArpackDN::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else if (mode == 2) {
        bmat[0] = 'G';
    } else if (mode == 3) {
        bmat[0] = 'G';
    } else if (mode == 4) {
        bmat[0] = 'G';
    }
}


void ArpackDN::set_shift(double sigmar_, double sigmai_)
{
    sigmar = sigmar_;
    sigmai = sigmai_;
}


int ArpackDN::compute_eigs(double* dr, double* di, double* v)
{
    vector<double> resid(n);
    vector<double> workd(3*n);
    vector<double> workl((3*ncv+6)*ncv);
    vector<double> workev(3*ncv);
    vector<fortran_int>    select(ncv);

    fortran_int lworkl = (3*ncv+6)*ncv;
    fortran_int iparam[11];
    fortran_int ipntr[14];
    fortran_int info   = 0;
    fortran_int ierr   = 0;
    fortran_int status = 0;
    fortran_int iter   = 0;
    fortran_int ido    = 0;

    // -- Clear temporaries
    fill(iparam, iparam+11, 0);
    fill(ipntr,  ipntr+14,  0);

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    ido = 0;
    while (1) {

        dnaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[0],
                &ncv, v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl,
                &info);

        if (ido == -1) {
            times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 1) {
            if (mode == 3 || mode == 4)
                times_OP2(&workd[ipntr[0]-1], &workd[ipntr[1]-1],
                          &workd[ipntr[2]-1]);
            else
                times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 2) {
            times_M(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else {
            break;
        }

        ++iter;
    }

    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in dnaupd_ = %d\n", info);
        dnaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        dneupd_(&rvec, "A", &select[0], dr, di, v, &n, &sigmar, &sigmai,
                &workev[0], bmat, &n, which, &nev, &tol, &resid[0], &ncv,
                v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl, &ierr);

        if (ierr != 0) {
            printf("Error in dneupd_ = %d\n", ierr);
            dneupd_info_string(ierr);
            status = -2;
        }
    }

    return status;
}


void ArpackDN::times_OP1(double*, double*)
{
}


void ArpackDN::times_OP2(double* x, double* opx, double*)
{
    times_OP1(x, opx);
}


void ArpackDN::times_M(double* x, double* Mx)
{
    copy(Mx, Mx+n, x);
}


// ---- ZN settings ----


ArpackZN::ArpackZN()
{
    set_mode(1);
    set_shift(0);
}


ArpackZN::~ArpackZN()
{
}


void ArpackZN::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else if (mode == 2) {
        bmat[0] = 'G';
    } else if (mode == 3) {
        bmat[0] = 'G';
    }
}


void ArpackZN::set_shift(dcomplex sigma_)
{
    sigma = sigma_;
}


int ArpackZN::compute_eigs(dcomplex* d, dcomplex* v)
{
    vector<dcomplex> resid(n);
    vector<dcomplex> workd(3*n);
    vector<dcomplex> workl((3*ncv+5)*ncv);
    vector<dcomplex> workev(3*ncv);
    vector<double>   rwork(ncv);
    vector<fortran_int>      select(ncv);

    fortran_int lworkl = (3*ncv+5)*ncv;
    fortran_int iparam[11];
    fortran_int ipntr[14];
    fortran_int info   = 0;
    fortran_int ierr   = 0;
    fortran_int status = 0;
    fortran_int iter   = 0;
    fortran_int ido    = 0;

    // -- Clear temporaries
    fill(iparam, iparam+11, 0);
    fill(ipntr,  ipntr+14,  0);

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    ido = 0;
    while (1) {

        znaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[0],
                &ncv, v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl,
                &rwork[0], &info);

        if (ido == -1) {
            times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 1) {
            if (mode == 3)
                times_OP2(&workd[ipntr[0]-1], &workd[ipntr[1]-1],
                          &workd[ipntr[2]-1]);
            else
                times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 2) {
            times_M(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else {
            break;
        }

        ++iter;
    }

    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in znaupd_ = %d\n", info);
        znaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        zneupd_(&rvec, "A", &select[0], d, v, &n, &sigma, &workev[0], bmat, &n,
                which, &nev, &tol, &resid[0], &ncv, v, &n, iparam, ipntr,
                &workd[0], &workl[0], &lworkl, &rwork[0], &ierr);
        if (ierr != 0) {
            printf("Error in zneupd_ = %d\n", ierr);
            zneupd_info_string(ierr);
            status = -2;
        }
    }

    return status;
}


void ArpackZN::times_OP1(dcomplex*, dcomplex*)
{
}


void ArpackZN::times_OP2(dcomplex* x, dcomplex* opx, dcomplex*)
{
    times_OP1(x, opx);
}


void ArpackZN::times_M(dcomplex* x, dcomplex* Mx)
{
    copy(Mx, Mx+n, x);
}
