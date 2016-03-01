#ifndef STL_WRAPPER_H
#define STL_WRAPPER_H

#include <vector> 

namespace STL_Wrapper
{

// return sorted element indicies in ascending order using lambda expression 
//
// after sorting, use the returned indicies to access v 
// to get sorted elements
//
// example: 
//
// v = [0, 3, 2] 
// for (auto i: SortIndicies(v)) 
//   cout << v[i] << endl;
//
// output: 
// 0
// 2
// 3
//
// Ref:
// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
vector<size_t> SortIndicies(const vector<T> &v) 
{

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// find unique element and trim
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
void PrintVectorContent(std::ostream &os, std::vector<T> &vec, const int cutoff=-1, const bool &verticalPrint=false)
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
