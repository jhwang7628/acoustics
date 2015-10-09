
#include <Eigen/Dense>
/*
 * A minimal class to store the data from wave solver. 
 */ 
class WaveSolverPointData
{

    private: 
        // N-by-3 matrix that stores positions of N cells
        Eigen::MatrixXd position_; 

        // M-by-1 matrix that stores time stamp of the simulation snapshot 
        Eigen::VectorXd time_; 

        // N-by-M matrix that stores M time steps of pressure computed on N cells.
        // the time stamp is store in "time" and position of the cell is stored in "position"
        Eigen::MatrixXd pressure_; 


    public: 
        WaveSolverPointData(){ }; 


        inline void SetPosition( const Eigen::MatrixXd & position )
        {
            if ( position.cols() != 3 ) throw runtime_error( "**ERROR** columns of position matrix is not 3" );
            position_ = position;
        }

        inline void PushPressure( const double time, const Eigen::MatrixXd pressure )
        {
            if ( pressure.rows() != position_.rows() || pressure.cols() != 1 ) throw runtime_error( "**ERROR** wrong number of cells" );

            pressure_.conservativeResize( Eigen::NoChange, pressure_.cols()+1 ); 
            pressure_.col( pressure_.cols()-1 ) = pressure; 
            time_.conservativeResize( time_.rows()+1, Eigen::NoChange );
        }


};
