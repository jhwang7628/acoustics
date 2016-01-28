#include <iostream> 
#include <IO/IO.h> 
#include <Eigen/Dense> 

using namespace std;
using namespace Eigen; 

const double sigma_x=0.2, sigma_t=1./10000., c=343.; 
const double sigma_x_sq=sigma_x*sigma_x; 
const double sigma_t_sq=sigma_t*sigma_t; 

double EvaluateGaussian( const Vector3d &x, const Vector3d &y, const double &t )
{
    const double normSq_y = y.squaredNorm();
    const double norm_x_y = (x-y).norm();
    const double delay = t-norm_x_y/c;

    //if (norm_x_y < 1E-4 || delay < 1E-4) 
    //    return 0.; 

    return 1./norm_x_y * exp( -normSq_y/sigma_x_sq/2.0 - pow(delay,2)/sigma_t_sq/2.0 );
}

int main()
{

    MatrixXd X; 
    IO::readMatrixXd(X, "head_listening.txt", IO::ASCII, 2); 

    const double minBound = -1.701; 
    const double maxBound =  1.701; 
    const int N = 100; 

    const double dx = (maxBound - minBound)/(double)N; 
    const double dx3 = pow(dx,3); 

    const double ccStart = minBound+dx/2.; 
    const VectorXd Time = VectorXd::LinSpaced(400,0.0,0.03);

    const int N_listen = X.rows(); 
    const int N_ts     = Time.size(); 

    MatrixXd Value(N_listen, N_ts); 

    int count = 0; 
    #pragma omp parallel for 
    for (int pp=0; pp<N_listen; pp++) 
    {
        #pragma omp critical 
        {
            std::cout << count << " "; 
            std::cout << std::flush; 
            count ++; 
        }
        VectorXd valueListen(N_ts); 
        valueListen.setZero();
        const Vector3d x = X.row(pp); 
        for (int tt=0; tt<N_ts; tt++) 
        {
            const double t = Time(tt); 
            for (int ii=0; ii<N; ii++) 
            {
                for (int jj=0; jj<N; jj++) 
                {
                    for (int kk=0; kk<N; kk++) 
                    {
                        Vector3d y((double)ii*dx, (double)jj*dx, (double)kk*dx); 
                        y.array() += ccStart; 
                        valueListen(tt) += EvaluateGaussian(x, y, t); 
                    }
                }
            }
        }
        #pragma omp critical 
        Value.row(pp) = valueListen; 


    }

    Value = Value*dx3/(4.0*M_PI); 

    IO::writeMatrixXd( Value, "value_x_cpp.txt", IO::ASCII ); 
    std::cout << std::endl;



    return 0; 
}
