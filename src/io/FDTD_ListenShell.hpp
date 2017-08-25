#ifndef FDTD_LISTENSHELL_HPP
#define FDTD_LISTENSHELL_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <string>
#include "linearalgebra/Vector3.hpp"
//##############################################################################
// Class FDTD_ListenShell
//##############################################################################
template<typename T>
class FDTD_ListenShell
{
    T                       _r; 
    Vector3<T>              _origin;
    std::vector<Vector3<T>> _points; 
    std::vector<Vector3<T>> _normls;
public: 
    static std::unique_ptr<FDTD_ListenShell<T>> Build(const std::string &file);
    FDTD_ListenShell() = default; 
    FDTD_ListenShell(const FDTD_ListenShell *rhs,
                     const Vector3<T> &neworigin, 
                     const T newradius)
    {
        _points = rhs->_points; 
        for (auto &p : _points)
            p = ((p - rhs->_origin)/rhs->_r*newradius) + neworigin; 
        _normls = rhs->_normls; 
        for (auto &n : _normls)
            n *= pow((newradius/_r),2);
    }
    inline int N(){return _points.size();}
    inline std::vector<Vector3<T>> &Points(){return _points;} 
};

//##############################################################################
//##############################################################################
template<typename T>
std::unique_ptr<FDTD_ListenShell<T>> FDTD_ListenShell<T>:: 
Build(const std::string &file)
{
    auto shell = std::make_unique<FDTD_ListenShell<T>>(); 
    auto FindStr = [](const std::string &str, const std::string &substr)
        {return str.find(substr)!=std::string::npos;};
    std::ifstream ifs(file.c_str()); 
    std::string line; 
    int N1=0, N2=0;
    while(std::getline(ifs,line))
    {
        if (FindStr(line, "START"))
        {
            std::getline(ifs,line); 
            std::istringstream iss(line);
            iss >> N1 >> N2; 
            break; 
        }
    }
    std::cout << "Reading " << N1 << " points and " << N2 << " normals "
              << "for constructing listening shell\n";

    // initialize
    shell->_points.resize(N1); 
    shell->_normls.resize(N2); 
    auto &points = shell->_points; 
    auto &normls = shell->_normls; 
    auto &origin = shell->_origin; 

    // read points and normals
    for (int ii=0; ii<N1; ++ii)
    {
        std::getline(ifs,line); 
        std::istringstream iss(line); 
        iss >> points[ii].x >> points[ii].y >> points[ii].z; 
    }
    for (int ii=0; ii<N2; ++ii)
    {
        std::getline(ifs,line); 
        std::istringstream iss(line); 
        iss >> normls[ii].x >> normls[ii].y >> normls[ii].z; 
    }

    // read origin and radius
    std::getline(ifs,line); 
    assert(FindStr(line, "ORIGIN")); 
    {
        std::getline(ifs,line);
        std::istringstream iss(line); 
        iss >> origin.x >> origin.y >> origin.z; 
    }
    std::getline(ifs,line); 
    assert(FindStr(line, "RADIUS")); 
    {
        std::getline(ifs,line); 
        std::istringstream iss(line); 
        iss >> shell->_r; 
    }
    std::getline(ifs,line); 
    assert(FindStr(line, "END"));

    return shell; 
}

#endif
