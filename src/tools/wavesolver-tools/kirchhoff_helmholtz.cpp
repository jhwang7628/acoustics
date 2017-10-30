#include <iostream> 
#include <fstream> 
#include <Eigen/Dense>
#include "io/FDTD_ListenShell.hpp"

#define DOUBLE_BYTES 8
const static int N_shell = 20480; 
const static int N_point = 2*N_shell; 
const static int N_steps = 100000; 
const static double DT = 1./1188186.;
const static double SOUND_SPEED = 343.; 
const static double DENSITY = 1.225; 

//##############################################################################
void Dp_Dn(const Eigen::VectorXd d_o,
           const Eigen::VectorXd d_i,
           const std::unique_ptr<FDTD_ListenShell<double>> &outShell, 
           const std::unique_ptr<FDTD_ListenShell<double>> &innShell,
           const int i,
           Eigen::VectorXd &dp_dn_i)
{
    const double h = outShell->Radius() - innShell->Radius();
    dp_dn_i = (d_o - d_i) / h;  // finite-difference
}

//##############################################################################
void Dp_Dt(const Eigen::VectorXd d_o_t,
           const Eigen::VectorXd d_o_tp1, 
           const std::unique_ptr<FDTD_ListenShell<double>> &outShell, 
           const int i,
           Eigen::VectorXd &dp_dt_i)
{
    dp_dt_i = (d_o_tp1 - d_o_t)/DT;
}

//##############################################################################
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "**Usage: " << argv[0] << " <read_steps>\n";
        exit(1); 
    }
    const std::string raw_data("/home/jui-hsien/code/acoustics/work/demo/listen_shell/data/test_all_audio.dat"); 
    const std::string out_shell("test_0_out_shell_config.txt"); 
    const std::string inn_shell("test_0_inn_shell_config.txt"); 
    const int read_steps = std::min(N_steps-1, atoi(argv[1]));

    const auto outShell = FDTD_ListenShell<double>::Build(out_shell); 
    const auto innShell = FDTD_ListenShell<double>::Build(inn_shell); 


    std::vector<Eigen::Vector3d> lisPts; 
    lisPts.push_back(Eigen::Vector3d(0., 0.025, 0.)); 
    lisPts.push_back(Eigen::Vector3d(0., 0.25, 0.)); 
    lisPts.push_back(Eigen::Vector3d(0., 2.5, 0.)); 
  
    const double coeff_1 = DENSITY / (4.0*M_PI); 
    const double coeff_2 = 1./(4.*M_PI*SOUND_SPEED);
    int count = 0;
    for (const auto &lisPt : lisPts)
    {
        std::cout << "\nLIS PT: " << lisPt.transpose() << std::endl; 
        const auto &points = outShell->Points(); 
        const auto &normls = outShell->Normls(); 
        Eigen::VectorXd outdata = Eigen::VectorXd::Zero(read_steps);
        std::ifstream data_stream(raw_data, std::ios::in|std::ios::binary); 
        data_stream.seekg(0, data_stream.end); 
        const long long int stream_length = data_stream.tellg(); 
        data_stream.seekg(0, data_stream.beg); 
        std::cout << "stream length = " << stream_length << std::endl;
        Eigen::VectorXd dp_dt, dp_dn;
        Eigen::VectorXd p_t_out(N_shell), p_t_inn(N_shell);
        Eigen::VectorXd p_tp1_out(N_shell), p_tp1_inn(N_shell); 
        const unsigned long long int unit_offset = DOUBLE_BYTES*N_shell; 
        data_stream.read((char*)p_t_out.data(), unit_offset); 
        data_stream.read((char*)p_t_inn.data(), unit_offset); 
        for (int tt=0; tt<read_steps; ++tt)
        {
            //std::cout << "\rStep " << tt << std::flush;
            const long long int here = data_stream.tellg(); 
            std::cout << "step = " << tt << " and here = " << here << std::endl;

            // read step ahead
            data_stream.read((char*)p_tp1_out.data(), unit_offset); 
            data_stream.read((char*)p_tp1_inn.data(), unit_offset); 
            // compute finite difference quantities
            Dp_Dt(p_t_out, p_tp1_out, outShell,           tt, dp_dt); 
            Dp_Dn(p_t_out, p_t_inn  , outShell, innShell, tt, dp_dn); 
            for (int ss=0; ss<N_shell; ++ss)
            {
                const auto   n = normls.at(ss); 
                const auto   r = Eigen::Vector3d(
                        lisPt[0] - points.at(ss).x,
                        lisPt[1] - points.at(ss).y,
                        lisPt[2] - points.at(ss).z); 
                const double A = n.norm(); 
                const double R = r.norm();
                const double R2 = pow(R,2);
                const double coeff_2_ = coeff_2 * (r[0]*n[0]+r[1]*n[1]+r[2]*n[2])/R/A; 
                const int delay = (int)(R/SOUND_SPEED);
                if (tt+delay >= outdata.size()) continue; 
                outdata(tt+delay) += 
                    (coeff_1 *dp_dn   (ss)/R  +
                     coeff_2_*dp_dt   (ss)/R  + 
                     coeff_2_* p_t_out(ss)/R2*SOUND_SPEED)*A; 
            }
        }
        std::cout << std::endl;
    
        std::ofstream ofs("point_"+std::to_string(count)+".dat", std::ios::binary); 
        ofs << lisPt.transpose() << std::endl;
        ofs << outdata.size() << std::endl; 
        ofs.write((char*)outdata.data(), DOUBLE_BYTES*outdata.size()); 
        ofs.close();
        count++; 
        data_stream.close();
    }

    return 0; 
}
