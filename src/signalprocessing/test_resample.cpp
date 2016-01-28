#include <fstream> 
#include <iostream> 
#include <Eigen/Dense> 
#include <cstdlib> 
#include <ctime> 

#include "resample.h" 


using namespace std; 

int main()
{

    double pi = 3.1415926;


    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(5000,0,1);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(t.size());

    /// test with sine
    //for (int ii=0; ii<y.size(); ii++) y(ii) = sin(2.0*pi*t(ii)/0.01);
   

    /// test with gaussian white noise 
    /* Setup constants */
    srand(time(0));
    const static int q = 15;
    const static float c1 = (1 << q) - 1;
    const static float c2 = ((int)(c1 / 3)) + 1;
    const static float c3 = 1.f / c1;
    
    /* random number in range 0 - 1 not including 1 */
    float random = 0.f;
    
    /* the white noise */
    float noise = 0.f;
    
    for (int i = 0; i < y.size(); i++)
    {
        random = ((float)rand() / (float)(RAND_MAX));
        y(i) = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
    }


    Eigen::MatrixXd data; 
    data.resize( t.size(), 2 ); 
    data.col(0) = t; 
    data.col(1) = y; 

    Eigen::MatrixXd data_naivedownsample = SIGNAL_PROCESSING::naiveDownsample( data.transpose(), 4 );
    Eigen::MatrixXd data_boxdownsample = SIGNAL_PROCESSING::BoxfilterDownsample( data.transpose(), 4 );



    data_boxdownsample.row(0) = data_naivedownsample.row(0);
    

    



    {
        ofstream of("tmp.txt"); 
        of << data << endl;
    }

    {
        ofstream of("tmp_naivedownsample.txt"); 
        of << data_naivedownsample.transpose() << endl;
    }

    {
        ofstream of("tmp_boxdownsample.txt");
        of << data_boxdownsample.transpose() << endl;
    }





    return 0; 
}
