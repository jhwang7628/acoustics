#ifndef VTKCONVERTER_H 
#define VTKCONVERTER_H

#include <iostream> 
#include <fstream> 
#include <Eigen/Dense> 

#include <boost/format.hpp> 

#include "utils/IO/IO.h" 


/* 
 * Convenient VTK converter
 */
class VTKConverter
{
    public: 

        enum DataType{
            BINARY,
            ASCII
        };

        static void VTKStructureGridWithScalarFromEigen( const Eigen::MatrixXd & grid, const Eigen::MatrixXd & scalarData, const string & filename, const string & scalarName, const DataType & type, const Eigen::Vector3i & dimension, const string & vtkname="Example" ) 
        {

            /// first write the grid
            VTKStructureGridFromEigen( grid, filename, type, dimension, vtkname );

            /// write the scalar data
            int Npoints = scalarData.rows();

            std::ofstream vtkstream; 
            if (type==ASCII) 
            {
                vtkstream.open(filename.c_str(), std::ios::out | std::ios::app);
                vtkstream<<boost::format( "\nPOINT_DATA %1%\n" ) %Npoints;
                vtkstream<<boost::format( "SCALARS %1% double\n" ) %scalarName;
                vtkstream<<"LOOKUP_TABLE default"<<std::endl;
                for (int ii=0; ii<Npoints; ii++) 
                    vtkstream<<boost::format( "%1% " ) %scalarData(ii,0);

            } 
            else 
            {

                vtkstream.open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
                vtkstream<<boost::format( "\nPOINT_DATA %1%\n" ) %Npoints;
                vtkstream<<boost::format( "SCALARS %1% double\n" ) %scalarName;
                vtkstream<<"LOOKUP_TABLE default"<<std::endl;
                double doubleBuffer;
                for (int ii=0; ii<Npoints; ii++) 
                {
                    doubleBuffer = IO::SwapEndianessDouble( scalarData(ii) );
                    vtkstream.write((char*) &(doubleBuffer), sizeof(double));
                }
            }
            vtkstream.close();

        }

        /* 
         * Convert N-by-3 Eigen matrix to VTK strucutred grid file output
         */
        static void VTKStructureGridFromEigen( const Eigen::MatrixXd & data, const string & filename, const DataType & type, const Eigen::Vector3i & dimension, const string & vtkname="Example" ) 
        {
            //printf( "write vtk structured grid from eigen matrix to file %s\n", filename.c_str() );

            std::ofstream vtkstream(filename.c_str(), std::ios::out | std::ios::trunc);
            int Npoints = data.rows();

            if (vtkstream) 
            {
                vtkstream<<"# vtk DataFile Version 2.0"<<"\n";
                vtkstream<<vtkname<<"\n";
                if (type==ASCII) 
                {
                    vtkstream<<"ASCII"<<"\n";
                    vtkstream.close();
                    vtkstream.clear();
                    vtkstream.open(filename.c_str(), std::ios::out | std::ios::app);
                    vtkstream<<"DATASET STRUCTURED_GRID"<<std::endl;
                    vtkstream<<boost::format( "DIMENSIONS %1% %2% %3%\n" ) %dimension[0] %dimension[1] %dimension[2];
                    vtkstream<<boost::format( "POINTS %1% double\n" ) %Npoints;
                    for (int ii=0; ii<Npoints; ii++) 
                    {
                        vtkstream<<boost::format( "%1% %2% %3% " ) %data(ii,0) %data(ii,1) %data(ii,2);
                    }

                } 
                else 
                {
                    vtkstream<<"BINARY"<<"\n";
                    vtkstream.close();
                    vtkstream.clear();
                    vtkstream.open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
                    vtkstream<<"DATASET STRUCTURED_GRID"<<std::endl;
                    vtkstream<<boost::format( "DIMENSIONS %1% %2% %3%\n" ) %dimension[0] %dimension[1] %dimension[2];
                    vtkstream<<boost::format( "POINTS %1% double\n" ) %Npoints;
                    for (int ii=0; ii<Npoints; ii++) 
                    {

                        for (int jj=0; jj<3; jj++)
                        {
                            double tmp = IO::SwapEndianessDouble( data(ii,jj) );
                            vtkstream.write((char*) &(tmp), sizeof(double));
                        }
                    }
                }
                vtkstream.close();
            } 
            else 
            {
                throw std::runtime_error( "**ERROR** Cannot open file stream");
            }

        }

};


#endif


