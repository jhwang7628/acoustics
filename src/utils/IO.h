//////////////////////////////////////////////////////////////////////
// IO.h: Some helpful I/O routines
//
//////////////////////////////////////////////////////////////////////

#ifndef IO_H
#define IO_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

#include "trace.h"

#include <ctime>
#include <cstring>
#include <time.h>
#include <stdlib.h>// for _MAX_PATH
#include <stdio.h>
#include <sstream>

#include <Eigen/Dense> 
#include <assert.h> 
#include <dirent.h>


#if defined(__unix__) || defined (__LINUX__)
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#endif

#ifdef WIN32
#include <windows.h>
#endif

using namespace std;

#define SDUMP(x)	" " << #x << "=[ " << x << " ] "

typedef enum {
    BINARY, 
    BINARY_V2,
    BINARY_V3, // only implemented for reading
    ASCII
} IOFileType;

class IO {
	public:
#if 0
		// Save an array of vectors
		static void saveVec3Data( const Vector3Array &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data( data.size() * 3 );

			for ( int i = 0; i < data.size(); i++ )
			{
				v_data( i * 3 + 0 ) = data[i][0];
				v_data( i * 3 + 1 ) = data[i][1];
				v_data( i * 3 + 2 ) = data[i][2];
			}

			v_data.write( fileName );
		}

		// Load an array of vectors
		static void loadVec3Data( Vector3Array &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data;

			v_data.read( fileName );

			data.resize( v_data.size() / 3 );

			for ( int i = 0; i < v_data.size() / 3; i++ )
			{
				data[i][0] = v_data( i * 3 + 0 );
				data[i][1] = v_data( i * 3 + 1 );
				data[i][2] = v_data( i * 3 + 2 );
			}
		}
#endif

		// Save an array of scalar data
		static void saveScalarData( const FloatArray &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data( data.size() );

			for ( size_t i = 0; i < data.size(); i++ ) {
				v_data(i) = data[i];
			}

			v_data.write( fileName );
		}

		// Load an array of scalar data
		static void loadScalarData( FloatArray &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data;

			v_data.read( fileName );

			data.resize( v_data.size() );

			for ( int i = 0; i < v_data.size(); i++ ) {
				data[i] = v_data(i);
			}
		}

		static vector<string> split( const string& s, const string& splitters )
		{
			vector<string> splat;
			int start = 0;
			while( 1 )
			{
				size_t occur = s.find_first_of( splitters, start );

				if( occur == string::npos ) {
					// we're done. add the last string
					splat.push_back( s.substr( start, string::npos ) );
					break;
				}
				else {
					splat.push_back( s.substr( start, occur-start ) );
					start = occur + 1;
				}	
			}

			return splat;
		}

        static void readMatrixXd( Eigen::MatrixXd &matout2D, const char *filename, IOFileType x ) 
        {
            switch (x) {

                case BINARY :
                    {

                        cout << "Read BINARY 2D data : " << filename << endl;

                        ifstream file(filename, ios::in | ios::binary | ios::ate); 

                        if (!file)  
                        {
                            cerr << "** Cannot find file " << filename << "." << endl; 
                        }


                        streampos size = file.tellg(); 
                        cout << " - Queried file size : " << size << endl; 

                        int Nrow, Ncol;  

                        file.seekg (0, ios::beg); 
                        file.read((char*)&Nrow, sizeof(int)); 
                        file.read((char*)&Ncol, sizeof(int)); 

                        // streampos size_left = (int)size - 2*sizeof(int); 

                        matout2D.resize(Nrow, Ncol); 

                        cout << " - progress: " << endl; 
                        cout << "     0 %";
                        for (int ii=0; ii<Nrow; ii++)
                        {
                            if ( ii % 50 == 0 )
                            {
                                cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                                     << "     " << (double)ii/(double)Nrow*100 << " % ";
                            }
                            for (int jj=0; jj<Ncol; jj++) 
                                file.read((char*)&matout2D(ii,jj), sizeof(double)); 
                        }
                        cout << endl; 

                        file.close(); 
                        cout << " - Data read in has dimension = (" << Nrow << ", " << Ncol << ")" << endl;

                        break; 
                    }

                case BINARY_V3 : 
                    {
                        cout << "Read BINARY Vector3 data : " << filename << endl;

                        ifstream file(filename, ios::in | ios::binary | ios::ate); 

                        if (!file)  
                        {
                            cerr << "cannot find file " << filename << "." << endl; 
                        }


                        streampos size = file.tellg(); 
                        cout << " - Queried file size : " << size << endl; 

                        int Nrow, Ncol;  

                        file.seekg (0, ios::beg); 
                        file.read((char*)&Nrow, sizeof(int)); 
                        file.read((char*)&Ncol, sizeof(int)); 

                        // streampos size_left = (int)size - 2*sizeof(int); 

                        matout2D.resize(Nrow, Ncol*3); 

                        cout << " - progress: " << endl; 
                        cout << "     0 %";
                        for (int ii=0; ii<Nrow; ii++)
                        {
                            if ( ii % 50 == 0 )
                            {
                                cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" 
                                     << "     " << (double)ii/(double)Nrow*100 << " % ";
                            }
                            for (int jj=0; jj<Ncol*3; jj++) 
                                file.read((char*)&matout2D(ii,jj), sizeof(double)); 
                        }
                        cout << endl; 

                        file.close(); 
                        cout << " - Data read in has dimension = (" << Nrow << ", " << Ncol << "*3)" << endl;

                        break; 
                    }

                case BINARY_V2 : 
                    {
                        cout << "Read BINARY Vector2 data : " << filename << endl;

                        ifstream file(filename, ios::in | ios::binary | ios::ate); 

                        if (!file)  
                        {
                            cerr << "cannot find file " << filename << "." << endl; 
                        }


                        streampos size = file.tellg(); 
                        cout << " - Queried file size : " << size << endl; 

                        int Nrow, Ncol;  

                        file.seekg (0, ios::beg); 
                        file.read((char*)&Nrow, sizeof(int)); 
                        file.read((char*)&Ncol, sizeof(int)); 

                        // streampos size_left = (int)size - 2*sizeof(int); 

                        matout2D.resize(Nrow, Ncol*2); 

                        cout << " - progress: " << endl; 
                        cout << "     0 %";
                        for (int ii=0; ii<Nrow; ii++)
                        {
                            if ( ii % 50 == 0 )
                            {
                                cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" 
                                     << "     " << (double)ii/(double)Nrow*100 << " % ";
                            }
                            for (int jj=0; jj<Ncol*2; jj++) 
                                file.read((char*)&matout2D(ii,jj), sizeof(double)); 
                        }
                        cout << endl; 

                        file.close(); 
                        cout << " - Data read in has dimension = (" << Nrow << ", " << Ncol << "*2)" << endl;

                        break; 
                    }
                case ASCII :
                    {
                        cout << "Read ASCII data : " << filename << endl; 

                        ifstream ifs(filename); 

                        if (!ifs) 
                        {
                            cerr << "Cannot open file " << filename << ". Returning. " << endl; 
                            return; 
                        }

                        string line; 
                        int Nrow=0, Ncol; 
                        while (getline(ifs, line))
                        {
                            if (Nrow % 500 == 0) { cout << Nrow << " lines read. "<< endl;}
                            Nrow ++; 

                            istringstream iss(line); 

                            vector<double> tmp; 
                            double buffer; 

                            Ncol = 0; 
                            while (iss >> buffer)
                            {
                                Ncol ++; 
                                tmp.push_back(buffer); 
                            }

                            // resize and push the vector content into eigen
                            // matrices
                            matout2D.conservativeResize(Nrow, Ncol); 
                            for (unsigned int ii=0; ii<tmp.size(); ii++) 
                            {
                                matout2D(Nrow-1, ii) = tmp[ii]; 
                            }

                        }

                        ifs.close(); 

                        break; 


                    }
            }

            cout << "Read complete." << endl; 

        }

		//----------------------------------------
		//  Tries to create a directory of the given path.
		//  Returns false if the directory already existed.
		//----------------------------------------
		enum CREATE_DIR_RESULT { OK, EXISTS, FAIL };

#if defined( WIN32 )

		static enum CREATE_DIR_RESULT create_dir( const string& name )
		{
			wchar_t* wideStr = NULL;
			size_t numWideChars = 0;
			BOOL ok = false;

			wideStr = new wchar_t[ name.size() ];
			mbstowcs_s( &numWideChars, wideStr, name.size()+1,
									name.c_str(), _TRUNCATE );

			string temp = wchar2string( wideStr );
			ok = CreateDirectory( temp.c_str(), NULL );
			delete wideStr;

			if( ok )
			{
				return OK;
			}
			else
			{
				cerr << "** Creating directory " << SDUMP(name) << " failed!" << endl;
				// TODO - use GetLastError to get more detailed errror info
				return FAIL;
			}
		}

#elif defined( __unix__ )

		static enum CREATE_DIR_RESULT create_dir( const string& name )
		{
			int status = mkdir( name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

			if( status == 0 )
			{
				return OK;
			}
			else
			{
				switch( errno )
				{
					case EEXIST: return EXISTS;
					default:
					{
						cerr << "** Creating directory " << SDUMP(name);
						cerr << " failed!" << endl;
						return FAIL;
					}
				}
			}
		}

#endif

};

#endif
