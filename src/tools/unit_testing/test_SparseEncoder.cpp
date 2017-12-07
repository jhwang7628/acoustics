#include <iostream> 
#include <fstream>
#include "utils/SimpleTimer.h"
#include "parser/ImpulseResponseParser.h"
#include "wavesolver/FDTD_Objects.h"
#include "wavesolver/PML_WaveSolver_Settings.h"
#include "wavesolver/FDTD_RigidSoundObject.h"
#include "linearalgebra/SPARSE_MATRIX.h"

//##############################################################################
// Function Main
//##############################################################################
int main(int argc, char **argv)
{
    if (argc < 3) 
    {
        std::cerr << "**Usage: " << argv[0] << " <xml_file> <steps>\n";
        exit(1); 
    }
        
    const std::string xml_file(argv[1]);
    const int steps = atoi(argv[2]);
    auto parser = std::make_shared<ImpulseResponseParser>(xml_file); 

    auto settings = std::make_shared<PML_WaveSolver_Settings>();
    auto objects = std::make_shared<FDTD_Objects>(); 
    parser->GetSolverSettings(settings); 
    parser->GetObjects(settings, objects); 

    auto object = 
        std::dynamic_pointer_cast<FDTD_RigidSoundObject>(objects->GetPtr(0));
    object->InitializeSparseModalEncoder(false); // don't use projected U

    const Eigen::MatrixXd &U_mass_normalized = object->GetEigenVectors(); 
    const int N_modes = U_mass_normalized.cols(); 
    const int N_dofs  = U_mass_normalized.rows(); 
    Eigen::VectorXd a_exact(N_dofs), a_aprox(N_dofs); 
    //auto encoder = object->GetModalEncoder(); 

    //{
    //    std::ofstream stream("U.txt"); 
    //    for (int ii=0; ii<U_mass_normalized.rows(); ii++)
    //        stream << U_mass_normalized.row(ii) << std::endl; 
    //    stream.close();
    //}

    Eigen::MatrixXd U(N_dofs, N_modes);
    Eigen::MatrixXd R(N_modes, N_modes); 
    R.setZero();
    {
        std::ifstream stream("Q.txt"); 
        std::string line; 
        int ii=0;
        while(std::getline(stream, line))
        {
            std::istringstream iss(line); 
            for (int jj=0; jj<N_modes; ++jj)
                iss >> U(ii, jj);
            ++ii; 
        }
    }

    {
        std::ifstream stream("R.txt"); 
        std::string line; 
        int ii=0;
        while(std::getline(stream, line))
        {
            std::istringstream iss(line); 
            for (int jj=0; jj<N_modes; ++jj)
                iss >> R(ii, jj);
            ++ii; 
        }
    }
    //std::cout << "|| U - QR || = " << (U_mass_normalized - U*R).norm() << std::endl;

    Eigen::VectorXd eigenFreqs(N_modes); 
    for (int ii=0; ii<N_modes; ++ii)
        eigenFreqs[ii] = object->GetModeFrequency(ii); 
    auto encoder = std::make_shared<SparseModalEncoder>(U, eigenFreqs);
    Eigen::VectorXd wi(N_modes); 
    wi.setOnes(); 
    wi /= wi.sum();
    //encoder->SetWeights(wi);

    SimpleTimer t_exact, t_aprox; 
    SimpleTimer t_encode, t_decode;
    SimpleTimer *encode_timers = new SimpleTimer[5];

    Eigen::VectorXd force_0 = Eigen::VectorXd::Zero(N_modes); 
    Eigen::VectorXd force_1 = Eigen::VectorXd::Ones(N_modes); 

    std::ofstream f_sparsity("sparsity_pattern.txt"); 
    std::ofstream f_encoderq("encoder_q.txt"); 
    std::ofstream f_Dq("Dq.txt"); 
    std::ofstream f_a("a.txt"); 
    for (int ii=0; ii<steps; ++ii)
    {
        if (ii == 0)
            object->AdvanceModalODESolvers(1, false, &force_1);
        else
            object->AdvanceModalODESolvers(1, false, &force_0);
        object->UpdateQPointers();

        const Eigen::VectorXd &q = object->GetModalAcceleration(); 
        const Eigen::VectorXd q_prime = R*q; 

        // raw acc compute
        t_exact.Start(); 
        const Eigen::VectorXd tmp = U*q_prime; 
        for (int jj=0; jj<U.rows(); ++jj)
            a_exact(jj) = tmp(jj); 
        t_exact.Pause(); 

        // decode
        t_aprox.Start();
        t_encode.Start();
        encoder->Encode(q_prime, encode_timers); 
        t_encode.Pause();
        t_decode.Start();
        for (int jj=0; jj<U.rows(); ++jj)
            a_aprox(jj) = encoder->Decode(jj); 
        t_decode.Pause();
        t_aprox.Pause();

        const Eigen::VectorXd &Dq = encoder->GetDQ(); 
        const Eigen::SparseVector<double> &dq = encoder->GetdQ(); 

        const Eigen::VectorXd &q_e = encoder->GetQExact();
        const Eigen::VectorXd &q_a = encoder->GetQAprox();
        const Eigen::VectorXd dq_dense = dq; 
        const int M = dq_dense.size(); 
        Eigen::VectorXi sparsity_pattern = Eigen::VectorXi::Zero(M); 
        for (int ii=0; ii<M; ++ii)
        {
            if (fabs(dq_dense(ii)) > 1E-12)
                sparsity_pattern(ii) = 1; 
            f_sparsity << sparsity_pattern(ii) << " "; 
            f_encoderq << std::setprecision(16);
            f_encoderq << q_e(ii) << " " << q_a(ii) << " "; 
            f_Dq << Dq(ii) << " "; 
        }
        f_a << a_exact(0) << " " << a_aprox(0) << std::endl;
        f_sparsity << "\n"; 
        f_encoderq << "\n";
        f_Dq << "\n"; 
    } 
    f_sparsity.close();
    f_encoderq.close();
    f_Dq.close();
    f_a.close();

    const double te = t_exact.DurationAverage(); 
    const double ta = t_aprox.DurationAverage(); 
    std::cout << "timing exact = " << te << std::endl; 
    std::cout << "timing aprox = " << ta << std::endl; 
    std::cout << "timing ratio = " << te/ta << std::endl;
    std::cout << "timing encode = " << t_encode.DurationAverage() << std::endl;
    for (int ii=0; ii<5; ++ii)
        std::cout << " - part " << ii << ": " << encode_timers[ii].DurationAverage() << std::endl;
    std::cout << "timing decode = " << t_decode.DurationAverage() << std::endl;

    return 0;
}
