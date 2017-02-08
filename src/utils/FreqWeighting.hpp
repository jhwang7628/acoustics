#ifndef _FREQ_WEIGHTING_HPP
#define _FREQ_WEIGHTING_HPP
#include <iostream>
#include "interp/CSpline.hpp"

//##############################################################################
// Namespace FreqWeighting
//   Utilities for weighting frequency curves
//##############################################################################
namespace FreqWeighting
{
    extern const double EqL40PHonDB[30];
    extern const double EqL40PHonFreq[30];
    extern const double Eql0PhonDB[63];
    extern const double Eql0PhonFreq[63];

    // lowest dB for 0 Phon: -9.1552dB @ 3575Hz
    extern const int    ISO226_0phon_ref_index;
    extern const double ISO226_dB_ref;
    extern const double ISO226_freq_ref;

    //###########################################################################
    // Class ISO226_Interpolator
    //   This is a singleton class
    //###########################################################################
    class ISO226_Interpolator
    {
        private: 
            static ISO226_Interpolator *_instance; 
            CSpline<double, false>      _spline; 

            ISO226_Interpolator(); 
            inline double DBToDecimel(const double &db) const
            {
                return pow(10, db/20.);
            }

        public: 
            static inline ISO226_Interpolator *Instance()
            {
                return _instance; 
            }
            double Weight(const double freq) const;
            void PrintAllWeights(const double &start, const double &stop, 
                                 const double &interval, std::ostream &stream
                                )const; 
            ///// debug /////
            void Test() const;
    }; 
};

#endif // _FREQ_WEIGHTING_HPP

