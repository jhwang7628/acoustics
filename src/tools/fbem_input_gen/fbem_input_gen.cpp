/*
 * =====================================================================================
 *
 *       Filename:  fbem_input_gen.cpp
 *
 *    Description:  Generate fastBEM input files for
 *                  FastBEM Acoustics V2.2.0
 *
 *        Version:  1.0
 *        Created:  11/02/10 00:59:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "io/TetMeshReader.hpp"
#include "utils/macros.h"
#include "utils/print_msg.h"

static const char FAST_BEM_VERSION[] = "2.0.8";
const double SOUND_SPEED_AIR         = 343.;
const double AIR_DENSITY             = 1.184;
static const int NRULE               = 3;
static const int NGAUSS              = 1;

using namespace std;

namespace fs = boost::filesystem;

static int          modeId    = -1;
static double       density   = -1.;
static double       cutFreq   = 20000.; // cut-off frequency
static string       modesFile = "";
static string       tetFile   = "";
static string       outFPtn   = "";

static int                      nFixed;
static FixVtxTetMesh<double>    tetMsh;
static TriangleMesh<double>     surfMsh;
static vector<int>              vidS2T;

static int              n3;
static int              nModes;
static vector<double>   eigval;
static vector<double>   eigvec;
static vector<double>   freqs;

static void parse_cmd(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "display help message")
        ("modes,m", po::value<string>(&modesFile), "The eigen mode file (.modes file)")
        ("density,d", po::value<double>(&density), "The material density")
        ("tet,t", po::value<string>(&tetFile), "The tet mesh file")
        ("out,o", po::value<string>(&outFPtn), "The output file pattern (e.g. input-%d.txt)")
        ("cut-freq,c", po::value<double>(&cutFreq), "Cutting-off frequency (default 20KHz)")
        ("mode-id,i", po::value<int>(&modeId), "The mode ID (optional)");

    po::variables_map vm;
    store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if ( vm.count("help") )
    {
        cout << "Usage: " << fs::path(argv[0]).filename() << " [options]" << endl;
        cout << "Create FastBEM input file for FastBEM " << FAST_BEM_VERSION << endl;
        cout << desc << endl;
        cout << "Created by Changxi Zheng (cxzheng@cs.cornell.edu)" << endl << endl;
        exit(0);
    }

    if ( modesFile.empty() || tetFile.empty() || outFPtn.empty() )
    {
        PRINT_ERROR("Please specify all the required file names\n");
        cerr << desc << endl;
        exit(1);
    }

    if ( density <= 0. )
    {
        PRINT_ERROR("Please specify the (positive) density value\n");
        cerr << desc << endl;
        exit(1);
    }
}

static void load_mesh()
{
    PRINT_MSG("Load file: %s\n", tetFile.c_str());
    if ( boost::iends_with(tetFile, ".tet") ) {
        //TetMeshLoader_Double::load_mesh(tetFile.c_str(), tetMsh);
        FV_TetMeshLoader_Double::load_mesh(tetFile.c_str(), tetMsh);
    } else {
        PRINT_ERROR("Unknown tet mesh format!\n");
        SHOULD_NEVER_HAPPEN(2);
    }

    map<int,int> idMap;
    tetMsh.extract_surface(&surfMsh);
    tetMsh.surface_id_map(idMap);
    vidS2T.resize(idMap.size());
    const std::map<int,int>::const_iterator end = idMap.end();
    for(std::map<int,int>::const_iterator it = idMap.begin();it != end;++ it)
        vidS2T[it->second] = it->first;

    nFixed = tetMsh.num_fixed_vertices();
    PRINT_MSG(" %d fixed vertices are detected\n", nFixed);
}

static void load_modes()
{
    PRINT_MSG("Load eigenmodes from %s\n", modesFile.c_str());

    ifstream fin(modesFile.c_str(), ios::binary);
    fin.read((char *)&n3, sizeof(int));
    fin.read((char *)&nModes, sizeof(int));
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file: %s\n", modesFile.c_str());
        SHOULD_NEVER_HAPPEN(2);
    }

    eigval.resize(nModes);
    freqs.resize(nModes);
    fin.read((char *)&eigval[0], sizeof(double)*nModes);
    int nmds = 0;
    for(;nmds < nModes;++ nmds)
    {
        freqs[nmds] = sqrt(eigval[nmds] / density) * 0.5 * M_1_PI;
        if ( freqs[nmds] > cutFreq ) break;
    }
    PRINT_MSG(" %d modes are in audible range\n", nmds);
    nModes = nmds;

    eigvec.resize(nmds*n3);
    fin.read((char *)&eigvec[0], sizeof(double)*eigvec.size());
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file: %s\n", modesFile.c_str());
        SHOULD_NEVER_HAPPEN(2);
    }
    fin.close();
}

/*
 * \param mid mode ID
 * \param vid vertex ID in tet mesh
 */
static inline Vector3d get_mode_disp(int mid, int vid)
{
    if ( vid < nFixed ) return Vector3d(0., 0., 0.);

    int idptr = n3*mid + (vid - nFixed) * 3;
    return Vector3d(eigvec[idptr], eigvec[idptr+1], eigvec[idptr+2]);
}

// Writes a list of all frequencies
static void write_frequencies(const char *fileName)
{
    ofstream fout(fileName);
    if ( fout.fail() ) {
        PRINT_ERROR("Failed to open file %s for writing\n", fileName);
        SHOULD_NEVER_HAPPEN(2);
    }

    fout << setprecision(16);
    fout << nModes << " eigenvalues are read" << endl;

    for ( int i = 0; i < nModes; i++ ) {
        double unscaled = freqs[ i ];
        unscaled *= 2 * M_PI;
        unscaled *= unscaled;
        unscaled *= density;

        fout << "Freq# " << i << " : " << freqs[ i ] << " Hz "
             << unscaled << endl;
    }

    fout.close();
}

static void write_fbem_input(int mid)
{
    char fname[128];
    sprintf(fname, outFPtn.c_str(), mid);

    ofstream fout(fname);
    if ( fout.fail() )
    {
        PRINT_ERROR("Fail to open file %s to write\n", fname);
        SHOULD_NEVER_HAPPEN(2);
    }

    fout << setprecision(16);
    fout << "Acoustic radiation input data generated by fbem_input_gen" << endl;
    fout << "Complete 1" << endl;
    fout << "Full 0 0.d0" << endl;

    // No. of elements, nodes, field points and cells
    fout << surfMsh.num_triangles() << ' ' << surfMsh.num_vertices() << " 0 0" << endl;
    // incident waves
    fout << "0 0" << endl;
    // point sources (monopole and dipole)
    fout << "0 0" << endl;
    // Speed of sound and medium density (the last zero is the wavenumber k ratio)
    fout << SOUND_SPEED_AIR << ' ' << AIR_DENSITY << " 2.d-5 1.d-12 0" << endl;
    // the current frequency
    fout << freqs[mid] << ' ' << freqs[mid] << " 1 0 0" << endl;
    // the beta/nrule/ngauss currently is fixed
    fout << "0 " << NRULE << ' ' << NGAUSS << " 0 One" << endl;

    // write vertices' positions
    fout << "$ Nodes:" << endl;
    const vector<Point3d>& vtx = surfMsh.vertices();
    for(size_t i = 0;i < vtx.size();++ i)
        fout << i+1 << ' ' << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z << endl;

    // Write out elements and boundary conditions
    fout << "$ Elements and Boundary Conditions:" << endl;
    const vector<Tuple3ui>& tgl = surfMsh.triangles();
    for(size_t i = 0;i < tgl.size();++ i)
    {
        fout << i+1 << ' ' 
             << tgl[i][0]+1 << ' '
             << tgl[i][2]+1 << ' '
             << tgl[i][1]+1 << " 2 ";   // Neumann bc always
        Vector3d vsum(0, 0, 0);
        for(int j = 0;j < 3;++ j)
            vsum += get_mode_disp(mid, vidS2T[tgl[i][j]]);
        vsum *= 1./3.;

        Vector3d nml = Triangle<double>::normal(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        nml.normalize();
        // output the normal velocity    (disp*freq => velocity)
        fout << "(0, " << (nml.dotProduct(vsum)*freqs[mid]) << ')' << endl;
    }

    fout << "$ Field Points:" << endl;
    fout << "$ Field Cells:" << endl;
    fout << "$ End of file" << endl;

    fout.close();
}

int main(int argc, char* argv[])
{
    parse_cmd(argc, argv);

    load_mesh();
    load_modes();

    if ( modeId >= 0 ) {
        if ( modeId >= nModes )
        {
            PRINT_WARNING("Given mode ID(%d) is out of range\n", modeId);
            return 0;
        }
        write_fbem_input(modeId);
    } else {
        boost::progress_display progress(nModes, cout, 
                "Generate FastBEM input files:\n    ", 
                "    ", "    ");
        for(int i = 0;i < nModes;++ i, ++ progress) {
            write_fbem_input(i);
        }

        write_frequencies( "freqs.txt" );
    }

    return 0;
}

