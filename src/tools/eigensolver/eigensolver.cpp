/******************************************************************************
 * Eigen-solver running on local computer
 *
 * eigensolver <K-mat-file> <M-mat-file> <# eig> <outfile>
 * It reads the K matrix and M matrix in given files, solve the generalized 
 * eigen-value problem, Kv = lambda Mv, and store the eigenvector and eigenvalues
 * into the given file.
 *
 *****************************************************************************/
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <boost/program_options.hpp>

#include "utils/nano_timer.h"

#include "logging/logging.h"

#include <arlsmat.h>
#include <arlgsym.h>
#include <arlssym.h>

#if 0
#include <arusmat.h>
#include <arugsym.h>
#endif

using namespace std;

static vector<int> ptrrow[2];
static vector<int> idxcol[2];
static vector<double> data[2];
static ARluSymMatrix<double>* mat[2] = {NULL, NULL};
#if 0
static ARumSymMatrix<double>* mat[2] = {NULL, NULL};
#endif
static string stiffMFile, massMFile, outFile;
static int numEigv = 200;
static double EV_THRESHOLD = 1.;
static bool verbose = false;

static uint8_t read_csc_dmatrix(const char* file, 
        std::vector<int>& ptrrow,
        std::vector<int>& idxcol,
        std::vector<double>& data,
        int& nrow, int& ncol)
{
    if ( verbose ) 
    {
        LOGGING_INFO("Load CSC matrix file: %s", file);
    }

    ifstream fin(file, ios::binary);
    if ( fin.fail() )
    {
        cerr << "read_csc_dmatrix:: Cannot open file " << file << " to read" << endl;
        return 255;
    }

    uint8_t ret;
    fin.read((char *)&ret, sizeof(uint8_t));
    if ( ret != 1 )
    {
        cerr << "read_csc_dmatrix:: The matrix data should be in double format" << endl;
        return 255;
    }

    int n;
    fin.read((char *)&ret, sizeof(uint8_t));
    fin.read((char *)&nrow, sizeof(int));
    fin.read((char *)&ncol, sizeof(int));
    fin.read((char *)&n,    sizeof(int));

    if ( (ret & 1) && (nrow != ncol) ) // symmetric
    {
        cerr << "read_csc_dmatrix:: Symmetric matrix should be square" << endl;
        return 255;
    }

    ptrrow.resize(nrow+1);
    idxcol.resize(n);
    data.resize(n);
    fin.read((char *)&ptrrow[0], sizeof(int)*(nrow+1));
    fin.read((char *)&idxcol[0], sizeof(int)*n);
    fin.read((char *)&data[0],   sizeof(double)*n);

    fin.close();
    return ret;
}

static void parse_cmd(int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description genericOpt("Generic options");
    genericOpt.add_options()
            ("help,h", "display help information");
    po::options_description configOpt("Configuration");
    configOpt.add_options()
            ("neig,n", po::value<int>(&numEigv)->default_value(200), 
                    "Number of smallest eigenvalues to compute")
            ("thresh,t", po::value<double>(&EV_THRESHOLD)->default_value(1.), 
                    "Lower bound of the eigenvalues to store in the file")
            ("stiff,s", po::value<string>(&stiffMFile)->default_value(""),
                    "Name of the stiffness matrix file")
            ("mass,m", po::value<string>(&massMFile)->default_value(""),
                    "Name of the mass matrix file")
            ("out,o", po::value<string>(&outFile)->default_value(""),
                    "Name of the output modes file")
            ("verbose,v", "Display details");
    // use configure file to specify the option
    po::options_description cfileOpt("Configure file");
    cfileOpt.add_options()
            ("cfg-file", po::value<string>(), "configuration file");

    po::options_description cmdOpts;
    cmdOpts.add(genericOpt).add(configOpt).add(cfileOpt);

    po::variables_map vm;
    store(po::parse_command_line(argc, argv, cmdOpts), vm);
    if ( vm.count("cfg-file") )
    {
        ifstream ifs(vm["cfg-file"].as<string>().c_str());
        store(parse_config_file(ifs, configOpt), vm);
    }
    po::notify(vm);

    if ( vm.count("help") )
    {
        printf("Usage: %s [options] \n", argv[0]);
        cout << cmdOpts << endl;
        exit(0);
    }
    verbose = vm.count("verbose") > 0;

    if ( massMFile.empty() )
    {
        cerr << "Specify mass matrix file" << endl;
        exit(1);
    }

    if ( stiffMFile.empty() ) 
    {
        cerr << "Specify stiffness matrix file" << endl;
        exit(1);
    }

    if ( outFile.empty() ) 
    {
        cerr << "Specify the output file" << endl;
        exit(1);
    }

    if ( verbose )
    {
        LOGGING_INFO("=============== Problem Summary ===============");
        LOGGING_INFO("Mass Matrix:                %s", massMFile.c_str());
        LOGGING_INFO("Stiffness Matrix:           %s", stiffMFile.c_str());
        LOGGING_INFO("# of eigenvalues:           %d", numEigv);
        LOGGING_INFO("Lower bound of eigenvalues: %lf", EV_THRESHOLD);
        LOGGING_INFO("===============================================");
    }
}

struct EIGV_CMP_
{
    const double* eigv; // eigenvalue
    EIGV_CMP_(const double* ev):eigv(ev) {}

    inline bool operator ()(int a, int b)
    {   return eigv[a] < eigv[b]; }
};

/*
 * File Format:
 * <int>: size of the eigen problem
 * <int>: # of eigenvalues
 * <eigenvalues>
 * <eigenvec_1>
 * ...
 * <eigenvec_n>
 */
static void write_eigenvalues(int nev, int nsz, 
        const int* ids, const double* eval, const double* evec, 
        const char* file)
{
    ofstream fout(file, std::ios::binary);
    if ( !fout.good() ) 
    {
        cerr << "write_eigenvalues::Cannot open file " << file << " to write" << endl;
        return;
    }

    // size of the eigen-problem. Here square matrix is assumed.
    fout.write((char *)&nsz, sizeof(int));
    int stid = 0;
#if 0
    while (stid < nev && eval[ids[stid]] < EV_THRESHOLD ) ++ stid;
#endif

    int nevout = nev - stid;    // # of eigenvalue to output
    fout.write((char *)&nevout, sizeof(int));

    // output eigenvalues
    for(int vid = stid;vid < nev;++ vid)
    {
        fout.write((const char*)&eval[ids[vid]], sizeof(double));
        printf("ev#%3d:  %lf\n", vid, eval[ids[vid]]);
    }
    // output eigenvectors
    for(int vid = stid;vid < nev;++ vid) {
        fout.write((const char*)&evec[ids[vid]*nsz], sizeof(double)*nsz);
    }

    fout.close();
}

int main(int argc, char* argv[])
{
    LoggingManager::instance().set_logging_level(LOG_INFO);
    parse_cmd(argc, argv);

    int nrow1, ncol1;
    uint8_t c = read_csc_dmatrix(stiffMFile.c_str(), ptrrow[0], 
            idxcol[0], data[0], nrow1, ncol1);
    if ( c == 255 )
    {
        cerr << "Fail to read the matrix file: " << stiffMFile << endl;
        exit(1);
    }
    if ( !(c & 1) )
    {
        cerr << "Symmetric matrix is excepted from file " << stiffMFile << endl;
        exit(1);
    }

#if 0
    // FIXME:
    cout << "Stiffness matrix rows: " << endl;
    for ( int i = 0; i < ptrrow[0].size(); i++) {
        cout << "    " << ptrrow[0][i] << endl;
    }

    cout << endl << "Stiffness matrix columns: " << endl;
    for ( int i = 0; i < idxcol[0].size(); i++) {
        cout << "    " << idxcol[0][i] << endl;
    }

    cout << endl << "Stiffness matrix data: " << endl;
    for ( int i = 0; i < data[0].size(); i++) {
        cout << "    " << data[0][i] << endl;
    }
#endif

    int nrow2, ncol2;
    c = read_csc_dmatrix(massMFile.c_str(), ptrrow[1], idxcol[1],
            data[1], nrow2, ncol2);
    if ( c == 255 )
    {
        cerr << "Fail to read the matrix file: " << massMFile << endl;
        exit(1);
    }
    if ( !(c & 1) || !(c & 2) )
    {
        cerr << "Symmetric positive definite matrix is excepted from file " << massMFile << endl;
        exit(1);
    }
    if ( nrow2 != nrow1 )
    {
        cerr << "Two matrices should have the same size" << endl;
        exit(1);
    }

    if ( numEigv <= 0 || numEigv > nrow1-2 )
    {
        cerr << "number of eigenvalues is out of range: maximum=" << nrow1-2 << endl;
        exit(1);
    }

    //// the readin CSC format is 1-based for all the matrix coordinate
    //// transform them into 0-based in order to be used in Arpack++
    for(int i = 0;i < 2;++ i)
    {
        for(size_t j = 0;j < ptrrow[i].size();++ j) -- ptrrow[i][j];
        for(size_t j = 0;j < idxcol[i].size();++ j) -- idxcol[i][j];
    }

    mat[0] = new ARluSymMatrix<double>(nrow1, (int)data[0].size(), 
            &data[0][0], &idxcol[0][0], &ptrrow[0][0]);
    mat[1] = new ARluSymMatrix<double>(nrow2, (int)data[1].size(),
            &data[1][0], &idxcol[1][0], &ptrrow[1][0]);
#if 0
    mat[0] = new ARumSymMatrix<double>(nrow1, (int)data[0].size(), 
            &data[0][0], &idxcol[0][0], &ptrrow[0][0]);
    mat[1] = new ARumSymMatrix<double>(nrow2, (int)data[1].size(),
            &data[1][0], &idxcol[1][0], &ptrrow[1][0]);
#endif

    cout << "Running eigensolver with shift " << EV_THRESHOLD << endl;

    //// solve eigen problem 
    char which[] = "LM";
    /*
    ARluSymStdEig<double> solver(numEigv, *mat[0], 0, which); //0, 0.1);
    */
    //ARluSymGenEig<double> solver('S', numEigv, *mat[0], *mat[1], 0, which);
    ARluSymGenEig<double> solver('S', numEigv, *mat[0], *mat[1], EV_THRESHOLD, which);
#if 0
    ARluSymStdEig<double> solver(numEigv, *mat[0], EV_THRESHOLD, which);
    //solver.SetShiftInvertMode(EV_THRESHOLD);
#endif
#if 0
    ARluSymGenEig<double> solver('S', numEigv, *mat[0], *mat[1], EV_THRESHOLD, which);
#endif

    LOGGING_INFO("max #iter:  %d", solver.GetMaxit());
    LOGGING_INFO("dimension:  %d", solver.GetN());
    LOGGING_INFO("tolenrence: %f", solver.GetTol());

    if ( verbose ) {
        LOGGING_INFO("Start solving eigen problem %d X %d ...", nrow1, ncol1);
    }
    double st = GetMilliTimed();
    solver.FindEigenvectors();

    LOGGING_INFO("#iter used: %d, time used: %lfs", 
            solver.GetIter(), GetMilliTimed() - st);

    int ncov = solver.ConvergedEigenvalues();
    const double* eval = solver.RawEigenvalues();
    const double* evec = solver.RawEigenvectors();
    vector<int> sortids(ncov);
    for(int i = 0;i < ncov;++ i) sortids[i] = i;

    //// sort the eigen value, discard the very low frequency ones
    sort(sortids.begin(), sortids.end(), EIGV_CMP_(eval));

    if ( verbose ) {
        LOGGING_INFO("Wrote eigenvalues and eigenvectors to file: %s", 
                outFile.c_str());
    }
    write_eigenvalues(ncov, nrow1, &sortids[0], eval, evec, outFile.c_str());

    delete mat[0];
    delete mat[1];
    return 0;
}

