#ifndef GEOMETRY_SIMPLEPOINTSET_HPP
#define GEOMETRY_SIMPLEPOINTSET_HPP

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "utils/IO/IO.h"

template<typename T> 
class Vertex
{
    private: 
        T _x; 
        T _y; 
        T _z; 

    public: 
             
        Vertex(T x, T y, T z) : _x(x), _y(y), _z(z) {}
        Vertex() : _x(0), _y(0), _z(0) {}

        void getPosition(vector<T> &pos) 
        {
            pos.resize(3); 
            pos[0] = _x; 
            pos[1] = _y; 
            pos[2] = _z; 
        }

        friend std::ostream & operator << (std::ostream & os, const Vertex<T>& obj) 
        {
            os << "(" << obj._x << ", " << obj._y << ", " << obj._z << ")"; 
            return os;
        }

};


template<typename T> 
class SimplePointSet
{
    public: 

        std::vector<Vertex<T> > vlist; 


        SimplePointSet() {}
        SimplePointSet(const char * filename)
        {
            loadPointSet(filename); 
        }


        /* 
         * Assume point set in ascii file without header, 
         * and is 3-dimensional. 
         */
        void loadPointSet(const char * filename)
        {
            Eigen::MatrixXd tmp; 
            IO::readMatrixXd(tmp,filename, IO::BINARY);
            //IO::readMatrixXd(tmp, filename, BINARY_V3); 
            vlist.resize(tmp.rows()); 

            for (int ii=0; ii<tmp.rows(); ii++) 
            {
                Vertex<T> vert(tmp(ii,0),tmp(ii,1),tmp(ii,2)); 
                vlist[ii] = vert; 
            }
        }

        friend std::ostream & operator << (std::ostream & os, const SimplePointSet<T>& sps) 
        {

            for (int ii=0; ii<sps.vlist.size(); ii++)
            {
                os << sps.vlist[ii] << std::endl; 
            }

            return os; 
        }

};


#endif 


