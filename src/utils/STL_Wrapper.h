#ifndef STL_WRAPPER_H
#define STL_WRAPPER_H

#include <vector> 

namespace STL_Wrapper
{

template <typename T> 
void VectorSortAndTrimInPlace(std::vector<T> &vec)
{
    std::sort(vec.begin(), vec.end());
    auto last = std::unique(vec.begin(), vec.end());
    vec.erase(last, vec.end());
}

template <typename T> 
std::vector<T> VectorSortAndTrim(std::vector<T> vec)
{
    VectorSortAndTrimInPlace(vec);
    return vec; 
}

template <typename T> 
void PrintVectorContent(std::ostream &os, std::vector<T> &vec, const bool &verticalPrint=false, const int cutoff=-1)
{
    const int half_cutoff = cutoff/2;
    const int N = vec.size();

    if (vec.size()<=cutoff || cutoff<0)
    {
        if (verticalPrint)
            for (auto &i : vec) os << i << std::endl; 
        else
        {
            for (auto &i : vec) os << i << " " << std::flush; 
            os << std::endl; 
        }
    }
    else 
    {
        if (verticalPrint)
        {
            for (int ii=0; ii<half_cutoff; ii++) 
                os << vec[ii] << std::endl;
            os << "..." << std::endl; 
            for (int ii=N-half_cutoff; ii<N; ii++) 
                os << vec[ii] << std::endl;
        }
        else
        {
            for (int ii=0; ii<half_cutoff; ii++) 
                os << vec[ii] << " " << std::flush;
            os << "..." << " "; 
            for (int ii=N-half_cutoff; ii<N; ii++) 
                os << vec[ii] << " " << std::flush;
            os << std::endl;
        }
    }
}

template <typename T> 
void RemoveInPlace(std::vector<T> &vec, const T &value) 
{
    vec.erase(std::remove(vec.begin(),vec.end(),value),vec.end());
}

};





#endif
