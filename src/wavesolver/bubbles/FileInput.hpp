#ifndef _FILE_INPUT_HPP
#define _FILE_INPUT_HPP

#include <cmath>
#include <map>
#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Dense>

#include "Mesh.hpp"

typedef Eigen::Matrix<double, 1, 1> VelocityValue;
typedef std::map<int, std::vector<VelocityValue>> SurfaceVelocityData;

struct FileNames
{
    std::string meshFile;
    std::string datFile; // helmholtz solution file
    std::string freqFile; // frequency info file
};

struct BubbleInputInfo
{
    int bubNum;
    double freq; // Negative for bad/dead bubbles
    std::complex<double> xfrIn; // TODO: unnecessary here
};

static std::vector<BubbleInputInfo>
parseFreqFile(const std::string &fName)
{
    using namespace std;
    vector<BubbleInputInfo> bubInfo;

    ifstream in(fName.c_str());

    string line;

    bubInfo.clear();

    // First line is time
    getline(in, line);

    // Read first real line
    getline(in, line);

    vector<string> data;
    while (!line.empty())
    {
        boost::split(data, line, boost::is_any_of(" \n"), boost::token_compress_on);

        if (data.size() < 2)
        {
            break;
        }

        // TODO: freq file bubble number does not match mesh bubble number
        BubbleInputInfo curInfo;
        curInfo.bubNum = atoi(data[0].c_str()) + 1; // change to 1 based indexing here

        if (data.size() == 2)
        {
            // Bad bubble
            curInfo.freq = -1;
        }
        else
        {
            // Good bubble

            curInfo.freq = atof(data[1].c_str());

            vector<string> xfrData;
            boost::split(xfrData, data[2], boost::is_any_of("(),"), boost::token_compress_on);

            curInfo.xfrIn = complex<double>(atof(xfrData[1].c_str()), atof(xfrData[2].c_str()));
        }

        bubInfo.push_back(curInfo);

        getline(in, line);
    }

    return bubInfo;
}

static std::map<double, FileNames>
parseFileNames(const std::string &dataDir)
{
    // Files we need to load:
    // 1. frequency files
    // 2. mesh files
    // 3. surface vibration data

    using namespace std;
    using namespace boost::filesystem;

    vector<string> meshFiles;
    boost::regex gmsh("gmsh.*\.msh");

    for (directory_iterator i(dataDir), endIt; i != endIt; ++i)
    {
        // Skip if not a file
        if (!is_regular_file( i->status() )) continue;

        boost::smatch what;

        // Skip if not a match
        if (!boost::regex_match( i->path().filename().native(), what, gmsh ))
        {
            //std::cout << i->path().filename().native() << " did not match" << std::endl;
            continue;
        }

        // Match
        meshFiles.push_back( i->path().native() );
    }

    // Sort the mesh files
    std::sort(meshFiles.begin(), meshFiles.end());

    // Now build the output
    map<double, FileNames> output;

    for (const string & m : meshFiles)
    {
        string p = path(m).parent_path().native();

        // parse time out
        // hacky, should really use boost split or something
        string fn = path(m).filename().native();
        string meshTime = fn.substr(5, 10);
        double t = boost::lexical_cast<double>(meshTime);

        FileNames f;
        f.meshFile = m;

        f.datFile = p + string("/bemOutput/helmSolution-") + meshTime + string(".dat");
        f.freqFile = p + string("/bemOutput/info-") + meshTime + string(".txt");

        output[t] = f;
    }

    return output;
}

static SurfaceVelocityData
loadSurfaceDatFile(const std::vector<BubbleInputInfo> &bubInfo,
                   const std::string &fileName,
                   const Mesh &mesh)
{
    std::ifstream in(fileName.c_str());
    SurfaceVelocityData output;

    // Count number of fluid surface triangles
    int surfTris = mesh.m_surfTris.size();
    int airTris = 0;

    for (auto type : mesh.m_triType)
    {
        if (type == Mesh::FLUID_AIR)
            ++airTris;
    }


    double rho = 1.184; // Density of air

    for (int i = 0; i < bubInfo.size(); ++i)
    {
        // skip bad bubbles
        if (bubInfo[i].freq < 0) continue;

        double omega = 2 * M_PI * bubInfo[i].freq;

        std::complex<double> factor(0, -omega * rho);
        std::complex<double> val;
        double r,c;

        // Read the dirichlet data and discard it
        for (int j = 0; j < mesh.m_vertices.size(); ++j)
        {
            in.read((char*)&r, sizeof(r));
            in.read((char*)&c, sizeof(c));
        }

        output[bubInfo[i].bubNum].resize(surfTris, Eigen::Matrix<double, 1, 1>::Zero());

        std::vector<Eigen::Matrix<double, 1, 1>> &curOutput = output[bubInfo[i].bubNum];

        // Read the velocity data and store it
        for (int j = 0; j < airTris; ++j)
        {
            in.read((char*)&r, sizeof(r));
            in.read((char*)&c, sizeof(c));

            val.real(r);
            val.imag(c);

            // Convert from pressure gradient to velocity here
            curOutput.at(mesh.m_fluidToSurf.at(j)) << std::abs((val / factor));
        }

        //double maxVal = 0;
        //for (auto val : curOutput)
        //{
        //    if (std::fabs(val(0)) > maxVal)
        //    {
        //        maxVal = std::fabs(val(0));
        //    }
        //}

        //std::cout << "Max loaded vel: " << maxVal << std::endl;
        //exit(1);
    }

    return output;
}

#endif // _FILE_INPUT_HPP

