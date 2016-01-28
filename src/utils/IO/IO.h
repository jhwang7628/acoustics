#ifndef IO_H
#define IO_H

#include <iostream> 
#include <fstream> 
#include <Eigen/Dense> 
#include <assert.h> 
#include <cstdio>
#include <dirent.h>
#include "StringHelper.h"
#include <sys/stat.h>
#ifdef USE_BOOST
#include <boost/filesystem.hpp> 
#endif // USE_BOOST

#define BUFFER_SIZE 500

using namespace std; 


class IO { 

    public: 


        enum FileType{
            BINARY, 
            //BINARY_V3, // only implemented for reading
            ASCII
        };

        static bool is_number( const std::string & s ) 
        {
            return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
        }

        template <typename T>
        static void readMatrixX( Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &matout2D, const char *filename, FileType x, const int VERBOSE=1 ) 
        {
            switch (x) 
            {
                case BINARY :
                {
                    if (VERBOSE >= 1) cout << "reading Eigen matrix ..." << flush; 
                        ifstream file(filename, ios::in | ios::binary | ios::ate); 
                        if (!file)  throw runtime_error("**ERROR** Cannot find file " + std::string(filename)); 

                        int Nrow, Ncol;  

                        file.seekg (0, ios::beg); 
                        file.read((char*)&Nrow, sizeof(int)); 
                        file.read((char*)&Ncol, sizeof(int)); 
                        matout2D.resize(Nrow, Ncol); 

                        for (int ii=0; ii<Nrow; ii++)
                            for (int jj=0; jj<Ncol; jj++) 
                                file.read((char*)&matout2D(ii,jj), sizeof(T)); 

                        file.close(); 

                        if (VERBOSE >= 1) cout << "OK\n data dimension : (" << Nrow << ", " << Ncol << ")" << endl;
                        break; 
                }
                case ASCII :
                {
                        ifstream ifs(filename); 
                        if (!ifs)  throw runtime_error("**ERROR** Cannot find file " + std::string(filename)); 

                        string line; 
                        int Nrow=0, Ncol; 
                        while (getline(ifs, line))
                        {
                            Nrow ++; 

                            if ( Nrow % 100 == 0 ) cout << "Nrow = " << Nrow << endl; 

                            istringstream iss(line); 

                            vector<T> tmp; 
                            T buffer; 

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
                                matout2D(Nrow-1, ii) = tmp[ii]; 
                        }
                        ifs.close(); 
                        break; 
                }
            }

            //cout << "Read complete." << endl; 

        }

        // Write int Eigen matrix data to a binary file. 
        // Format 
        template <typename T>
        static void writeMatrixX( const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &data, const char *filename, FileType x ) 
        { 
            int Nrow = (int) data.rows(); 
            int Ncol = (int) data.cols(); 
            ofstream of1(filename, ios::out | ios::binary);
            if (!of1)  throw runtime_error("**ERROR** Cannot find file " + std::string(filename)); 

            switch (x) {

                case BINARY :
                {
                        cout << "Write data in BINARY format to file : " << filename << endl;
                        cout << " - Data Dimension: " << Nrow << " x " <<  Ncol << endl; 
                        of1.write((char *) &Nrow, sizeof(int));
                        of1.write((char *) &Ncol, sizeof(int)); 

                        cout << " - Progress: " << endl; 
                        cout << "     " << "0 % "; 
                        for (int ii=0; ii<Nrow; ii++) 
                        {
                            if (ii%50 ==0) cout << "     " << (double)ii/(double)Nrow*100.0 << " % \r"; 
                            for (int jj=0; jj<Ncol; jj++) 
                                of1.write((char *) &data(ii,jj), sizeof(T));
                        }
                        cout << endl; 
                        break;
                }
                case ASCII : 
                {
                        cout << "Write data in ASCII format to file : " << filename << endl;
                        cout << " - Progress: " << endl; 
                        cout << "     " << "0 % "; 
                        for (int ii=0; ii<Nrow; ii++) 
                        {
                            if (ii%50 ==0) cout << "     " << (double)ii/(double)Nrow*100.0 << " % \r"; 
                            for (int jj=0; jj<Ncol; jj++)
                                of1 << data(ii,jj) << " "; 
                            of1<<endl;
                        }
                        of1.close();
                        break; 
                }
            }
            of1.close(); 


        }

        // Write int Eigen matrix data to a binary file. 
        // Format 
        static void writeMatrixXi( const Eigen::MatrixXi &data, const char *filename, FileType x ) 
        { 

            int Nrow = (int) data.rows(); 
            int Ncol = (int) data.cols(); 

            switch (x) 
            {

                case BINARY :
                    {
                        cout << "Write data in BINARY format to file : " << filename << endl;
                        ofstream of1(filename, ios::out | ios::binary);

                        if (!of1) 
                        {
                            cerr << "** file I/O (O) problem. Returning" << endl; 
                            return;
                        }

                        cout << " - Data Dimension: " << Nrow << " x " <<  Ncol << endl; 
                        of1.write((char *) &Nrow, sizeof(int));
                        of1.write((char *) &Ncol, sizeof(int)); 


                        cout << " - Progress: " << endl; 
                        cout << "     " << "0 % "; 
                        for (int ii=0; ii<Nrow; ii++) 
                        {
                            if (ii%50 ==0)
                            {
                                cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
                                cout << "     " << (double)ii/(double)Nrow*100.0 << " % "; 
                            }
                            for (int jj=0; jj<Ncol; jj++) 
                                of1.write((char *) &data(ii,jj), sizeof(int));
                        }
                        cout << endl; 

                        of1.close(); 
                        break;

                    }
                case ASCII : 
                    throw runtime_error( "**ERROR** not implemented " ); 
                    break; 
            }
        }

        // Write Eigen matrix data to a binary file. 
        // Format 
        static void writeMatrixXd( const Eigen::MatrixXd &data, const char *filename, FileType x ) 
        { 

            int Nrow = (int) data.rows(); 
            int Ncol = (int) data.cols(); 

            switch (x) {

                case BINARY :
                    {
                        cout << "Write data in BINARY format to file : " << filename << endl;
                        ofstream of1(filename, ios::out | ios::binary);

                        if (!of1) 
                        {
                            cerr << "** file I/O (O) problem. Returning" << endl; 
                            return;
                        }

                        cout << " - Data Dimension: " << Nrow << " x " <<  Ncol << endl; 
                        of1.write((char *) &Nrow, sizeof(int));
                        of1.write((char *) &Ncol, sizeof(int)); 


                        cout << " - Progress: " << endl; 
                        cout << "     " << "0 % "; 
                        for (int ii=0; ii<Nrow; ii++) 
                        {
                            if (ii%50 ==0)
                            {
                                cout << "\r";
                                cout << "     " << (double)ii/(double)Nrow*100.0 << " % "; 
                            }
                            for (int jj=0; jj<Ncol; jj++) 
                                of1.write((char *) &data(ii,jj), sizeof(double));
                        }
                        cout << endl; 

                        of1.close(); 
                        break;

                    }
                case ASCII : 
                    {
                        FILE *fs1;
                        fs1 = fopen(filename ,"w"); 
                        cout << "Write data in ASCII format to file : " << filename << endl;

                        cout << " - Progress: " << endl; 
                        cout << "     " << "0 % "; 
                        for (int ii=0; ii<Nrow; ii++) 
                        {
                            if (ii%50 ==0)
                            {
                                cout << "\r";
                                cout << "     " << (double)ii/(double)Nrow*100.0 << " % "; 
                            }
                            for (int jj=0; jj<Ncol; jj++)
                                fprintf(fs1, "%.16f ", data(ii,jj));

                            fprintf(fs1, "\n");
                        }
                        fclose(fs1);
                        cout << endl; 
                        break; 
                    }
            }
        }

        template <typename MatrixType> 
        static bool writeMatrix_csv(const MatrixType &data, const std::string &filename)
        {
            const int N_rows = data.rows(); 
            const int N_cols = data.cols(); 

            std::ofstream outFile(filename.c_str()); 
            if (!outFile) return false;
            for (int ii=0; ii<N_rows; ii++) 
            {
                for (int jj=0; jj<N_cols; jj++) 
                    outFile << data(ii,jj) << ","; 
                outFile << "\n"; 
            }
            return true; 
        }

        static void readMatrixXd( Eigen::MatrixXd &matout2D, const char *filename, FileType x, const int VERBOSE=1 ) 
        {
            switch (x) 
            {

                case BINARY :
                    {

                        if (VERBOSE >= 2)
                            cout << "Read BINARY 2D data : " << filename << endl;

                        ifstream file(filename, ios::in | ios::binary | ios::ate); 

                        if (!file)  
                        {
                            cerr << "** Cannot find file " << filename << "." << endl; 
                        }

                        int Nrow, Ncol;  

                        file.seekg (0, ios::beg); 
                        file.read((char*)&Nrow, sizeof(int)); 
                        file.read((char*)&Ncol, sizeof(int)); 

                        matout2D.resize(Nrow, Ncol); 

                        if (VERBOSE >= 2) 
                        {
                            cout << " - progress: " << endl; 
                            cout << "     0 %";
                        }
                        for (int ii=0; ii<Nrow; ii++)
                        {
                            if ( ii % 50 == 0 && VERBOSE >= 1 )
                            {
                                cout << "\r"
                                     << "     " << (double)ii/(double)Nrow*100 << " % ";
                            }
                            for (int jj=0; jj<Ncol; jj++) 
                                file.read((char*)&matout2D(ii,jj), sizeof(double)); 
                        }
                        if (VERBOSE >= 2) 
                            cout << endl; 

                        file.close(); 
                        if (VERBOSE >= 1) 
                            cout << " - Data read in has dimension = (" << Nrow << ", " << Ncol << ")" << endl;

                        break; 
                    }

                case ASCII :
                    {
                        //cout << "Read ASCII data : " << filename << endl; 

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
                            Nrow ++; 

                            if ( Nrow % 100 == 0 ) cout << "Nrow = " << Nrow << endl; 

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

            //cout << "Read complete." << endl; 

        }

        static Eigen::MatrixXd readMatrixXd( const int & Nrow, const int & Ncol, ifstream & file, FileType x ) 
        {

            Eigen::MatrixXd matout2D(Nrow, Ncol); 

            switch (x) 
            {
                case BINARY :
                    if (!file) throw runtime_error( "**ERROR** Input file stream corrupted. ");

                    for (int ii=0; ii<Nrow; ii++)
                        for (int jj=0; jj<Ncol; jj++) 
                            file.read((char*)&matout2D(ii,jj), sizeof(double)); 

                    file.close(); 
                    break; 
                default : throw runtime_error( "**ERROR** file type not recognized. ");
            }
            return matout2D; 

        }

        static Eigen::MatrixXi readMatrixXi( const int & Nrow, const int & Ncol, ifstream & file, FileType x ) 
        {

            Eigen::MatrixXi matout2D(Nrow, Ncol); 

            switch (x) 
            {
                case BINARY :
                    if (!file) throw runtime_error( "**ERROR** Input file stream corrupted. ");

                    for (int ii=0; ii<Nrow; ii++)
                        for (int jj=0; jj<Ncol; jj++) 
                            file.read((char*)&matout2D(ii,jj), sizeof(int)); 

                    file.close(); 
                    break; 
                default : throw runtime_error( "**ERROR** file type not recognized. ");
            }
            return matout2D; 

        }
        /**
         * Discrete read the ascii files. 
         * it is usually useful to skip a couple of lines of headers, 
         * and then discretely pick columns for reading. This method implements
         * this functionality. 
         *
         * rows: an integer specifying number of lines that will be SKIPPED. 
         * cols: a std vector that will be PICKED. (don't need to be ordered)
         *       if cols is a zero-vector, then all cols will be read. 
         *
         * Note: 
         *  1. rows is 1-based (specify size); cols is 0-based (specify position). 
         *  2. if cols contains repeated values it will only be read once. 
         */ 
        static void readMatrixXd_discrete_ASCII(Eigen::MatrixXd &matout2D, const char* filename, 
                                      int rows, vector<int> cols)
        {
            assert (rows >= 0); 

            ifstream ifs(filename); 

            if (!ifs) 
            {
                cerr << "** Cannot open file " << filename << ". Returning. " << endl; 
                return; 
            }




            cout << "Read ASCII data discretely : " << filename << endl << endl; 

            cout << "-------------------------------------------------- " << endl; 
            cout << "The following #rows from the file will be SKIPPED : " << endl; 
            cout << rows << endl; 
            cout << "-------------------------------------------------- " << endl; 
            cout << "The following columns from the file will be READ : " << endl; 
            for (unsigned int ii=0; ii<cols.size(); ii++)
                cout << cols[ii] << " "; 
            cout << endl; 
            cout << "-------------------------------------------------- " << endl; 


            string buffer; 
            for (int ii=0; ii<rows; ii++)
            {
                if (!getline(ifs, buffer))
                {
                    cerr << "** Oops, skipping too many rows?" << endl; 
                    return; 
                }
            }

            int Nrow=0, Ncol;  
            string line; 
            while (getline(ifs, line))
            {
                Nrow ++; 

                istringstream iss(line); 

                vector<double> tmp; 
                double buffer; 


                int count = 0; 
                Ncol = 0; 
                while (iss >> buffer)
                {
                    if (cols.size() > 0)
                    {
                        // simple linear search to check
                        for (unsigned int jj=0; jj<cols.size(); jj++) 
                        {
                            if (count == cols[jj]) 
                            {
                                Ncol ++; 
                                tmp.push_back(buffer); 
                                break; 
                            }
                        }
                    }
                    else  // read all
                    {
                        Ncol ++; 
                        tmp.push_back(buffer); 
                    }
                    count ++; 
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


        }

        /**
         * A function that makes life a bit easier. 
         */
        static bool findInitStr(const std::string & str, const std::string & substr)
        {
            return str.substr(0, substr.size()) == substr; 
        }

        static bool findInitStr(const std::string & str, const std::string & substr, std::string & reststr)
        {
            bool found = false; 
            if ( str.substr(0, substr.size()) == substr)
            {
                found = true; 
                reststr = str.substr(substr.size()-1);
            }
            else 
                reststr = str; 

            return found; 
        }

        static bool FindString(const std::string &str, const std::string &substr)
        {
            std::size_t found = str.find(substr); 
            if (found!=std::string::npos) 
                return true; 
            else 
                return false; 
        }


        /**
         * Split the input eigen matrix which was read with BINARY_V3 to three
         * seperate eigen matrices. 
         */ 
        static void splitDataV3(const Eigen::MatrixXd &D, Eigen::MatrixXd &mat_x,
                                                          Eigen::MatrixXd &mat_y, 
                                                          Eigen::MatrixXd &mat_z) 
        {
            assert (D.cols()%3 == 0);  // in case an incomplete read or wrong identifier when reading

            cout << "Spliting data from Vector3" << endl; 
            int Dactual = D.cols()/3; 

            mat_x.resize(D.rows(), Dactual); 
            mat_y.resize(D.rows(), Dactual); 
            mat_z.resize(D.rows(), Dactual); 

            for (int jj=0; jj<Dactual; jj++)
            {
                mat_x.col(jj) = D.col(jj*3); 
                mat_y.col(jj) = D.col(jj*3 + 1); 
                mat_z.col(jj) = D.col(jj*3 + 2); 
            }
        }

       static void ReadChar(std::ifstream &inFile) 
        {
            char charBuffer; 
            inFile.get(charBuffer); 
        }

        /** 
         * Extract file names in "directory" filenames of which starts with "startString", return results 
         * in filenames. Will attempt a sort on string names. 
         */ 
        static void listDirectory(const char *directory, const char *startString, vector<string> &filenames)
        {

            DIR *dp;
            struct dirent *dirp;

            if((dp  = opendir(directory)) == NULL) 
            {
                cerr << "Error(" << errno << ") opening " << directory << endl;
                exit(1);
            }
              
            while ((dirp = readdir(dp)) != NULL) 
            {
                string dname = string(dirp->d_name); 

                if (findInitStr(dname, startString))
                {
                    filenames.push_back(dname);
                }
            }

            // sort the file names
            std::sort(filenames.begin(), filenames.end());

            closedir(dp);

            if (filenames.size() == 0) 
                cout << "**WARNING** zero files with start string " << startString << " in directory " << directory << " read. " << endl; 

        }

        /** 
         * Extract file names in "directory" filenames of which starts with "startString", return results 
         * in filenames. Will attempt a sort on string names. 
         */ 
        static void listDirectory(const char *directory, const char *startString, const char *negateString, vector<string> &filenames)
        {

            DIR *dp;
            struct dirent *dirp;

            if((dp  = opendir(directory)) == NULL) 
            {
                cerr << "Error(" << errno << ") opening " << directory << endl;
                exit(1);
            }
              
            while ((dirp = readdir(dp)) != NULL) 
            {
                string dname = string(dirp->d_name); 

                if (findInitStr(dname, startString) && dname.find(negateString)==std::string::npos)
                {
                    filenames.push_back(dname);
                }
            }

            // sort the file names
            std::sort(filenames.begin(), filenames.end());

            closedir(dp);

            if (filenames.size() == 0) 
                cout << "**WARNING** zero files with start string " << startString << " in directory " << directory << " read. " << endl; 

        }

        // Split the string s using splitters.  
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

        /* 
         * strip / and return filepath
         * its not very robust as I only detect if first and last character is
         * slash. should be enough for internal use though
         */
        static std::string AssembleFilePath( const std::vector<std::string> & tokens )
        {
            if (tokens.size() == 1) 
                return tokens[0]; 

            std::string filepath = tokens[0]; 

            if ( *(tokens[0].end()-1) == '/' )
                filepath = std::string( tokens[0].begin(), tokens[0].end()-1 );

            for (size_t ii=1; ii<tokens.size(); ii++) 
            {
                auto b = ( *(tokens[ii].begin())=='/' ? tokens[ii].begin()+1 : tokens[ii].begin() );
                auto e = ( *(tokens[ii].end()-1)=='/' ? tokens[ii].end()  -1 : tokens[ii].end()   ); 
                filepath = filepath + "/" + string( b,e );
            }

            return filepath;

        }

        /* 
         * strip / and return filepath
         * its not very robust as I only detect if first and last character is
         * slash. should be enough for internal use though
         */
        static std::string AssembleFilePath( const std::string & token1, const std::string & token2 )
        {
            std::vector<std::string> tokens; 
            tokens.push_back(token1); 
            tokens.push_back(token2);

            return AssembleFilePath( tokens );

        }

        static std::string AssembleFilePath( const std::string & token1, const std::string & token2, const std::string &token3 )
        {
            std::vector<std::string> tokens; 
            tokens.push_back(token1); 
            tokens.push_back(token2);
            tokens.push_back(token3);

            return AssembleFilePath( tokens );

        }

        /* 
         * Swap Eidianess for 8-bytes double
         */ 
        static double SwapEndianessDouble(double d)
        {
            double a;
            unsigned char *dst = (unsigned char *)&a;
            unsigned char *src = (unsigned char *)&d;
        
            dst[0] = src[7];
            dst[1] = src[6];
            dst[2] = src[5];
            dst[3] = src[4];
            dst[4] = src[3];
            dst[5] = src[2];
            dst[6] = src[1];
            dst[7] = src[0];
        
            return a;
        }


        /* 
         * check if file exists
         *
         * ref: 
         * http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
         */ 
        static bool ExistFile (const std::string& name) 
        {
            struct stat buffer;   
            return (stat (name.c_str(), &buffer) == 0); 
        }

        static bool ExistFolder(const std::string &pathname)
        {
            struct stat info;

            if( stat( pathname.c_str(), &info ) != 0 )
                    return false; 
            else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
                    return true; 
            else
                    return false; 
        }



        /*
         * don't compile this part if boost is not found. useful when compile
         * lightweight
         */
#ifdef USE_BOOST

            static bool CreateDirectoryForce(const std::string &directoryName)
            {
                boost::filesystem::path dir(directoryName.c_str());
                if (IO::ExistFolder(directoryName))
                {
                    std::cout << "**WARNING** directory " << directoryName << " exists. do nothing." << std::endl;
                    return true; 
                }
                else 
                    return boost::filesystem::create_directory(dir); 
            }
            /** 
             * Extract file names in "directory" filenames of which starts with "startString", return results 
             * in filenames. Will attempt a sort on string names. 
             */ 
            static void listDirectoryMatch(const char *directory, const char *str, vector<string> &filenames)
            {

                DIR *dp;
                struct dirent *dirp;

                if((dp  = opendir(directory)) == NULL) 
                {
                    cerr << "Error(" << errno << ") opening " << directory << endl;
                    exit(1);
                }


                //cout << "Reading filename list from directory : " << directory << endl; 
                //cout << "The substring that is searched for : \"" << str << "\"" << endl; 

                while ((dirp = readdir(dp)) != NULL) 
                {
                    string dname = string(dirp->d_name); 

                    //std::size_t found = dname.find( str );
                    //if ( found != std::string::npos ) filenames.push_back( dname );
                    if ( StringHelper::matchSubstr( dname, str ) ) filenames.push_back( dname );
                }

                // sort the file names
                std::sort(filenames.begin(), filenames.end());

                closedir(dp);

                if (filenames.size() == 0) 
                    cout << "** WARNING ** There might be some problems, zero files read. " << endl; 
                //else 
                //    cout << "Job finished. " << filenames.size() << " files in directory are read. " << endl; 

            }

            /** 
             * Extract file names in "directory" filenames of which starts with "startString", return results 
             * in filenames. Will attempt a sort on string names. 
             */ 
            static std::vector<std::string> listDirectoryMatch( const char * directory, const char * str )
            {

                std::vector<std::string> filenames; 
                listDirectoryMatch( directory, str, filenames );

                return filenames; 
            }

            /** 
             * Extract file names in "directory" filenames of which starts with "startString", return results 
             * in filenames. Will attempt a sort on string names. 
             */ 
            static void listDirectoryG(const char *directory, const char *str, vector<string> &filenames)
            {

                DIR *dp;
                struct dirent *dirp;

                if((dp  = opendir(directory)) == NULL) 
                {
                    cerr << "Error(" << errno << ") opening " << directory << endl;
                    exit(1);
                }


                //cout << "Reading filename list from directory : " << directory << endl; 
                //cout << "The substring that is searched for : \"" << str << "\"" << endl; 

                while ((dirp = readdir(dp)) != NULL) 
                {
                    string dname = string(dirp->d_name); 

                    //std::size_t found = dname.find( str );
                    //if ( found != std::string::npos ) filenames.push_back( dname );
                    if ( StringHelper::searchSubstr( dname, str ) ) filenames.push_back( dname );
                }

                // sort the file names
                std::sort(filenames.begin(), filenames.end());

                closedir(dp);

                if (filenames.size() == 0) 
                    cout << "** WARNING ** There might be some problems, zero files read. " << endl; 
                //else 
                //    cout << "Job finished. " << filenames.size() << " files in directory are read. " << endl; 

            }

            /** 
             * Extract file names in "directory" filenames of which starts with "startString", return results 
             * in filenames. Will attempt a sort on string names. 
             */ 
                static std::vector<std::string> listDirectoryG( const char * directory, const char * str )
                {

                    std::vector<std::string> filenames; 
                    listDirectoryG( directory, str, filenames );

                    return filenames; 
                }
#endif 




}; 

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
#endif
