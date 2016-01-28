#ifndef RESAMPLE_H
#define RESAMPLE_H 

#include "Eigen/Dense" 
#include <iostream>


namespace SIGNAL_PROCESSING
{


/*
 * Downsample a signal columns by simply skipping sample points. 
 * this could cause ringing because of the high frequency part. 
 *
 * now in-place to save memory.
 */ 
void NaiveDownsample( Eigen::MatrixXd & data, const int & skipEvery )
{

    std::cout<<"downsample signal columns by skipping "<<skipEvery<<" steps."<<std::endl;

    const int newNCols = data.cols()/skipEvery; 

    for ( int ii=0; ii<newNCols; ii++ ) 
        data.col( ii ) = data.col( ii*skipEvery );  // should not over-write the useful data

    data.conservativeResize( data.rows(), newNCols );
}


/* 
 * Downsample a signal columns by first box filter it (taking the mean) 
 * skipEvery is essentially the width of the filter. 
 */
Eigen::MatrixXd BoxfilterDownsample( const Eigen::MatrixXd & data, const int & skipEvery ) 
{

    std::cout<<"downsample signal columns after box filter and skip "<<skipEvery<<" steps."<<std::endl;

    Eigen::MatrixXd data_downsampled; 
    data_downsampled.resize( data.rows(), (int)(data.cols()/skipEvery) ); 
    data_downsampled.setZero(); 

    for ( int ii=0; ii<data.cols()/skipEvery; ii++ )
    {
        Eigen::VectorXd colMean = Eigen::VectorXd::Zero( data.rows() ); 
        for ( int jj=0; jj<data.rows(); jj++ ) 
        {
            Eigen::MatrixXd tmp = data.block(jj,ii*skipEvery,1,skipEvery); 

            colMean(jj)=tmp.mean();
        }

        data_downsampled.col( ii ) = colMean; 
    }

    return data_downsampled;

}

/* 
 * Frequency domain downsample.
 * Take the FFT of the signal, chop-off the high frequency component and 
 * downsample it. 
 */ 
Eigen::MatrixXd FFTDownsample()
{
}
 





} // namespace SIGAL_PROCESSING

#endif
