#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
#include "igl/read_triangle_mesh.h"
#include "igl/write_triangle_mesh.h"
#include "Eigen/Dense"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"

//##############################################################################
// Function ParseFileID
//##############################################################################
int ParseFileID(const boost::filesystem::path &a)
{
    std::vector<std::string> tokens;
    boost::split(tokens, a.stem().string(), [](char c){return c == '_';});
    return std::stoi(tokens.at(1));
};

//##############################################################################
// Function ReadFilenames
//##############################################################################
int ReadFilenames(const std::string &dir,
                   std::vector<boost::filesystem::path> &filenames)
{
    using namespace boost::filesystem;
    path p(dir.c_str());
    for (directory_iterator it(p); it!=directory_iterator(); ++it)
        if (it->path().extension().string() == ".obj")
            filenames.push_back(it->path());
    std::sort(filenames.begin(), filenames.end(),
            [](const path &a, const path &b){
                return ParseFileID(a) < ParseFileID(b);});
    for (auto f : filenames)
        std::cout << f.string() << std::endl;
    return filenames.size();
}

//##############################################################################
// Function Main
//  Upsample the poser data
//##############################################################################
int main(int argc, char **argv)
{
#ifndef USE_IGL
    std::cerr << "**ERROR** Need libigl\n";
    return 1;
#endif
    if (argc < 6)
    {
        std::cout << "**Usage: " << argv[0]
                                 << " <obj_directory>"
                                 << " <old_frame_rate>"
                                 << " <grid_cell_size>"
                                 << " <output_obj_directory>"
                                 << " <output_obj_prefix>"
                                 << " [new_frame_rate]"
                                 << "\n";
        return 1;
    }

    // read all filenames and sort
    std::cout << "\nReading obj filenames ...\n";
    std::vector<boost::filesystem::path> filenames;
    ReadFilenames(std::string(argv[1]), filenames);
    const double old_framerate = atoi(argv[2]);
    const double dt = 1./(double)old_framerate;
    std::cout << "...Done. Files found: " << filenames.size() << std::endl;

    // read all objs and figure out max vertex speed
    std::cout << "\nReading all objs to find max vertex speed ...\n";
    int N_V, N_F;
    Eigen::MatrixXd V_old;
    Eigen::MatrixXd vel;
    double velmax = std::numeric_limits<double>::lowest();
    for (int ii=0; ii<filenames.size(); ++ii)
    {
        const auto &f = filenames.at(ii);
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;

        // Read mesh and check coherency
        igl::read_triangle_mesh(f.string().c_str(), V, F);
        if (ii==0)
        {
            N_V = V.rows();
            N_F = F.rows();
        }
        if (N_V != V.rows() || N_F != F.rows())
        {
            std::cerr << "**ERROR** mesh is not conherent\n";
            return 2;
        }
        else
        {
            std::cout << "  Read mesh is coherent: " << f.string() << std::endl;
        }

        if (ii>0)
        {
            vel = (V - V_old) / dt;
            for (int jj=0; jj<vel.rows(); ++jj)
                velmax = std::max(velmax, vel.row(jj).norm());
            std::cout << "    Max vel = " << velmax << std::endl;
        }
        V_old = V;
    }
    std::cout << "...Done. Max vertex speed: " << velmax << std::endl;

    // bound the new frame rate
    std::cout << "\nBounding framerate ...\n";
    const double h = atof(argv[3]);
    const double max_dt = h/velmax;
    const double min_framerate = 1./max_dt;
    int new_framerate = (int)std::ceil(min_framerate/5.0)*5;
    std::cout << "...Done. Maximum dt for cell size " << h
              << " is: " << max_dt << std::endl;
    std::cout << "  This corresponds to the minimum frame rate: "
              << 1./max_dt << std::endl;
    std::cout << "  Round it up to the new frame rate: "
              << new_framerate << std::endl;
    if (argc > 6)
    {
        new_framerate = std::max(new_framerate, atoi(argv[6]));
        std::cout << "  Input frame rate detected, taking max and the new frame rate is: "
                  << new_framerate << std::endl;
    }

    // make sure we want to proceed
    {
        std::cout << "\nProceed? ";
        std::string buf;
        std::getline(std::cin, buf);
    }

    // write new objs
    std::cout << "\nUpsampling and write to new directory ...\n";
    const std::string outdir(argv[4]);
    const std::string outpre(argv[5]);
    const int N_files = filenames.size();
    const double T[2] = {(double)ParseFileID(filenames.at(0        ))*dt,
                         (double)ParseFileID(filenames.at(N_files-1))*dt};
    const double new_dt = 1./(double)new_framerate;
    std::cout << "  Write to directory: " << outdir << std::endl;
    std::cout << "  Interpolation time-range: ["
              << T[0] << " " << T[1] << "]\n";
    double t = T[0];
    int count = 0;
    while (t < T[1])
    {
        std::cout << "  Interpolating for mesh at time: " << t << std::endl;
        // compute interpolation coeff
        const int idx_n = (int)((t-T[0])/dt);
        const int idx_p = std::min(idx_n+1, N_files-1);
        double alpha;
        if (idx_n >= N_files-1)
        {
            alpha = 1.0;
        }
        else
        {
            const double t_n = (double)idx_n*dt;
            const double t_p = (double)idx_p*dt;
            alpha = (t - t_n)/dt;
        }

        // interpolate mesh
        Eigen::MatrixXd V_n, V_p;
        Eigen::MatrixXd F;
        igl::read_triangle_mesh(filenames.at(idx_n).c_str(), V_n, F);
        igl::read_triangle_mesh(filenames.at(idx_p).c_str(), V_p, F);
        Eigen::MatrixXd V_interp = V_n*(1.0 - alpha) + V_p*alpha;

        const std::string filename(
                outdir+"/"+outpre+"_"+std::to_string(count)+".obj");
        std::cout << "  Writing new obj to file: " << filename << std::endl;
        igl::write_triangle_mesh(filename.c_str(), V_interp, F);

        t += new_dt;
        count += 1;
    }

    // write processing metadata to output dir
    {
        std::cout << "\nWriting processing metadata to output dir ...\n";
        std::vector<std::pair<std::string, std::string>> metaField;
        metaField.push_back(
                std::make_pair("old_frame_rate", std::to_string(old_framerate)));
        metaField.push_back(
                std::make_pair("new_frame_rate", std::to_string(new_framerate)));
        std::string filename(outdir+"/"+outpre+"_resample-metadata");
        std::ofstream stream(filename.c_str());
        if (stream)
            for (const auto &field : metaField)
                stream << field.first << " " << field.second << std::endl;
        stream.close();
        std::cout << "...Done. Written file: " << filename << std::endl;
    }

    return 0;
}
