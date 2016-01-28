#ifndef TEXTURE_H
#define TEXTURE_H 


#include <iostream>
#include <fstream>
#include <Eigen/Dense> 
#include "Wav.h"
#include "Distance.h" 

#include "omp.h"

namespace SIGNAL_PROCESSING
{ 



/* 
 * General class for sound texture
 */ 
class Sound
{ 


    enum ChannelType
    {
        MONO,
        STEREO
    }; 

    protected: 
        Eigen::MatrixXd data_; 
        ChannelType     dataType_; 


    public: 

        Sound( const std::string wavname )
        {
            printf( "reading wav file from %s \n", wavname.c_str() );
            try
            {
                data_ = SIGNAL_PROCESSING::readWavStereo( wavname ); 
                dataType_ = STEREO; 
            }
            catch (int e)
            {
                try 
                {
                    data_ = SIGNAL_PROCESSING::readWavMono( wavname ); 
                    dataType_ = MONO;
                }
                catch (int e)
                {
                    std::cerr << "**ERROR** Cannot read file " << wavname << std::endl;
                }
            }
        }

        inline Eigen::MatrixXd L()      const { return data_.col(0); }
        inline Eigen::MatrixXd R()      const { return isMono() ? data_.col(0) : data_.col(1); }
        inline Eigen::MatrixXd data()   const { return data_; } 
        inline int channels()           const { return data_.cols(); }
        inline int length()             const { return data_.rows(); } 
        inline bool isMono()            const { return dataType_==MONO ? true : false; } 

        friend std::ostream & operator << (std::ostream & os, const Sound& texture ) 
        {
            os << texture.data_ << std::endl; 
        }

}; 


class SoundTexture : public Sound 
{
    public : 

        enum DistanceMetric
        {
            L2
        };

    private : 

        Eigen::MatrixXd graph_; 

    public : 

        SoundTexture( const std::string wavname ) : Sound( wavname ) { } 


        /* 
         * Point-wise energy evaluation using window size 2k+1
         */
        void EvaluateL2Energy( const Eigen::VectorXd & L, const Eigen::VectorXd & R, const int & k )
        {
            std::cout << "Evaluating L2 Energy for the signal" << std::endl;

            const int N = L.size(); 
            std::cout << "matrix size : " << N << "-by-" << N << std::endl;
            graph_.setZero(N,N);

#pragma omp parallel for
            for ( int ii=0; ii<N; ii++ )
            {
                if ( (ii%10)==0 ) std::cout<< ii << std::endl;

                
                for ( int jj=ii; jj<N; jj++ )
                {
                    printf( "%u/%u ", ii, jj);
                    double energy_L_window = 0.0; 
                    double energy_R_window = 0.0; 
                    for ( int kk=-k; kk<=k; kk++ )
                    {
                        int index_i = (ii+kk+N) % N; 
                        int index_j = (jj+kk+N) % N; 

                        printf("%u:%u ", index_i, index_j );
                        energy_L_window += pow(L(index_i) - L(index_j), 2);
                        energy_R_window += pow(R(index_i) - R(index_j), 2);
                    }
                    graph_(ii,jj) += sqrt( energy_L_window ); 
                    graph_(ii,jj) += sqrt( energy_R_window ); 
                }
                printf("\n");
            }

        }

        /* 
         * wrapper for building the graph using some distance metric
         */
        Eigen::MatrixXd & BuildGraph()
        {
            std::cout << "build graph for sound texture" << std::endl;

            // debug
            EvaluateL2Energy( L().block(0,0,30,1), R().block(0,0,30,1), 10 );


            return graph_; 
        }




};



} // end namespace SIGNAL_PROCESSING



#endif 












