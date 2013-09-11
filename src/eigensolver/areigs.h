//ldoc
/**
% Simple C++ interface to ARPACK
 
[ARPACK](http://www.caam.rice.edu/software/ARPACK/) is a collection of
eigenvalue solver routines written in Fortran 77.  It uses a reverse
communication interface.  This C++ wrapper provides an alternate
interface based on virtual method calls.

To use the interface, one derives from `ArpackDS` (for real symmetric
problems), `ArpackDN` (for real nonsymmetric problems), or `ArpackZN`
(for complex nonsymmetric problems) and implements the virtual methods
`times_OP1`, `times_OP2`, and `times_M`.  The details of what these
methods do depends somewhat on the mode, but in general, `times_OP1`
computes the basic matrix-vector multiplication used to generate the
Krylov subspace, `times_M` multiplies by a mass matrix (in the case of
a generalized eigenvalue problem), and `times_OP2` does a composite
operation that takes advantage of a pre-computed product with a mass
matrix.

*/
//ldoc
#ifndef AREIGS_H
#define AREIGS_H

#include <complex>

//ldoc
/**
# General interface

The basic interface allows one to set or get:

- `n` : the matrix dimension
- `nev` : the number of desired eigenvalues
- `ncv` : the number of Krylov vectors to use
- `rvec` : whether or not eigenvectors should be computed
- `maxitr` : the maximum outer iterations allowed
- `which` : which eigenvalues are desired (`LA`, `SA`, `LM`, `SM`, `BE`)
- `tol` : relative accuracy tolerance
- `mode` : integer describing the mode in which the solver will be used
- `shift` : shift used in shift-invert or Cayley transform modes

The `mode` and `shift` setters and getters vary depend on the specific
mode that is used.
*/
class Arpack {
public:
    Arpack();
    virtual ~Arpack();

    void set_n(int n_);
    void set_nev(int nev_);
    void set_ncv(int ncv_);
    void set_rvec(int rvec_);
    void set_maxitr(int maxitr_);
    void set_which(const char* which_);
    void set_tol(double tol_);

    int    get_n()      { return n;      }
    int    get_nev()    { return nev;    }
    int    get_ncv()    { return ncv;    }
    int    get_rvec()   { return rvec;   }
    int    get_maxitr() { return maxitr; }
    char*  get_which()  { return which;  }
    double get_tol()    { return tol;    }

protected:
    typedef int fortran_int;

    fortran_int n;
    fortran_int nev;
    fortran_int ncv;
    fortran_int rvec;
    fortran_int mode;
    fortran_int maxitr;
    char bmat[2];
    char which[3];
    double tol;
};


/**
# Real symmetric interface

The valid `mode` flags for `ArpackDS` are

1.  Regular mode
2.  Regular inverse mode for GEP.  `OP` is $M^{-1}A$, `B` is $M$.
3.  Shift-invert mode for GEP.  `OP` is $(A-\sigma M)^{-1} M$, `B` is $M$.
4.  Buckling mode.  `OP` is $(A-\sigma M)^{-1} M$, `B` is $A$.
5.  Cayley mode.  `OP` is $(A-\sigma M)^{-1} (A+\sigma M)$, `B` is $M$.

*/
class ArpackDS : public Arpack {
public:
    ArpackDS();
    ~ArpackDS();
    int compute_eigs(double* d, double* v);
    void set_mode(int mode_);
    void set_shift(double sigma_);
protected:
    double sigma;
    virtual void times_OP1(double* x, double* opx);
    virtual void times_OP2(double* x, double* opx, double* Mx);
    virtual void times_M  (double* x, double* Mx);
};


/**
# Real nonsymmetric interface

For `ArpackDN`, the valid `mode` flags are

1.  Regular mode
2.  Regular inverse mode for GEP.  `OP` is $M^{-1}A$, `B` is $M$.
3.  Real shift-invert mode for GEP.  `OP` is $(A-\sigma M)^{-1} M$, `B` is $M$.
3.  Complex shift-invert for GEP.  `OP` is $\Re( (A-\sigma M)^{-1} M )$, `B` is $M$.
4.  Complex shift-invert for GEP.  `OP` is $\Im( (A-\sigma M)^{-1} M )$, `B` is $M$.

*/
class ArpackDN : public Arpack {
public:
    ArpackDN();
    ~ArpackDN();
    int compute_eigs(double* dr, double* di, double* v);
    void set_mode(int mode_);
    void set_shift(double sigmar_, double sigmai_ = 0);
    void set_shift(std::complex<double> sigma_) {
        set_shift(real(sigma_), imag(sigma_));
    }
protected:
    double sigmar, sigmai;
    virtual void times_OP1(double* x, double* opx);
    virtual void times_OP2(double* x, double* opx, double* Mx);
    virtual void times_M  (double* x, double* Mx);
};


/**
# Complex non-Hermitian interface

For `ArpackZN`, the valid `mode` flags are

1.  Regular mode
2.  Regular inverse mode for GEP.
3.  Shift-invert mode for GEP.

*/
class ArpackZN : public Arpack {
public:
    ArpackZN();
    ~ArpackZN();
    int compute_eigs(std::complex<double>* d, std::complex<double>* v);
    void set_mode(int mode_);
    void set_shift(std::complex<double> sigma_);
    void set_shift(double sigmar_, double sigmai_ = 0) {
        set_shift(std::complex<double>(sigmar_, sigmai_));
    }
protected:
    std::complex<double> sigma;
    virtual void times_OP1(std::complex<double>* x, 
                           std::complex<double>* opx);
    virtual void times_OP2(std::complex<double>* x, 
                           std::complex<double>* opx, 
                           std::complex<double>* Mx);
    virtual void times_M  (std::complex<double>* x, 
                           std::complex<double>* Mx);
};

//ldoc
#endif /* AREIGS_H */
