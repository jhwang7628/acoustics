#include <fstream>
#include "utils/FreqWeighting.hpp" 
namespace FreqWeighting
{
const double EqL40PHonDB[30] = {
    99.8539,
    93.9444,
    88.1659,
    82.6287,
    77.7849,
    73.0825,
    68.4779,
    64.3711,
    60.5855,
    56.7022,
    53.4087,
    50.3992,
    47.5775,
    44.9766,
    43.0507,
    41.3392,
    40.0618,
    40.0100,
    41.8195,
    42.5076,
    39.2296,
    36.5090,
    35.6089,
    36.6492,
    40.0077,
    45.8283,
    51.7968,
    54.2841,
    51.4859,
    100
};

const double EqL40PHonFreq[30] = {
    2.0000e+01,
    2.5000e+01,
    3.1500e+01,
    4.0000e+01,
    5.0000e+01,
    6.3000e+01,
    8.0000e+01,
    1.0000e+02,
    1.2500e+02,
    1.6000e+02,
    2.0000e+02,
    2.5000e+02,
    3.1500e+02,
    4.0000e+02,
    5.0000e+02,
    6.3000e+02,
    8.0000e+02,
    1.0000e+03,
    1.2500e+03,
    1.6000e+03,
    2.0000e+03,
    2.5000e+03,
    3.1500e+03,
    4.0000e+03,
    5.0000e+03,
    6.3000e+03,
    8.0000e+03,
    1.0000e+04,
    1.2500e+04,
    2e+04
};

const double Eql0PhonDB[63] = {
    7.6552e+01,
    7.0722e+01,
    6.5619e+01,
    5.9927e+01,
    5.5123e+01,
    4.9882e+01,
    4.5534e+01,
    4.1221e+01,
    3.7632e+01,
    3.3873e+01,
    3.0865e+01,
    2.7640e+01,
    2.5024e+01,
    2.2538e+01,
    2.0510e+01,
    1.8392e+01,
    1.6646e+01,
    1.4700e+01,
    1.3116e+01,
    1.1496e+01,
    1.0088e+01,
    8.6833e+00,
    7.5436e+00,
    6.2358e+00,
    5.1137e+00,
    3.9577e+00,
    3.0589e+00,
    2.1878e+00,
    1.4824e+00,
    7.9040e-01,
    3.0292e-01,
    -1.0709e-01,
    -3.0265e-01,
    -2.9268e-01,
    -1.0258e-02,
    6.3675e-01,
    1.0335e+00,
    3.1790e-01,
    -1.1863e+00,
    -2.7611e+00,
    -4.1116e+00,
    -5.6831e+00,
    -7.0462e+00,
    -8.3209e+00,
    -9.0260e+00,
    -9.1552e+00,
    -8.4944e+00,
    -6.8665e+00,
    -4.4829e+00,
    -6.2586e-01,
    3.2817e+00,
    7.2672e+00,
    9.8291e+00,
    1.1116e+01,
    1.0476e+01,
    7.8876e+00,
    8.3813e+00,
    2.1454e+01,
    4.1800e+01,
    5.2974e+01,
    6.4000e+01,
    7.5951e+01,
    8.9900e+01
};

const double Eql0PhonFreq[63] = {
    2.0000e+01,
    2.2500e+01,
    2.5000e+01,
    2.8250e+01,
    3.1500e+01,
    3.5750e+01,
    4.0000e+01,
    4.5000e+01,
    5.0000e+01,
    5.6500e+01,
    6.3000e+01,
    7.1500e+01,
    8.0000e+01,
    9.0000e+01,
    1.0000e+02,
    1.1250e+02,
    1.2500e+02,
    1.4250e+02,
    1.6000e+02,
    1.8000e+02,
    2.0000e+02,
    2.2500e+02,
    2.5000e+02,
    2.8250e+02,
    3.1500e+02,
    3.5750e+02,
    4.0000e+02,
    4.5000e+02,
    5.0000e+02,
    5.6500e+02,
    6.3000e+02,
    7.1500e+02,
    8.0000e+02,
    9.0000e+02,
    1.0000e+03,
    1.1250e+03,
    1.2500e+03,
    1.4250e+03,
    1.6000e+03,
    1.8000e+03,
    2.0000e+03,
    2.2500e+03,
    2.5000e+03,
    2.8250e+03,
    3.1500e+03,
    3.5750e+03,
    4.0000e+03,
    4.5000e+03,
    5.0000e+03,
    5.6500e+03,
    6.3000e+03,
    7.1500e+03,
    8.0000e+03,
    9.0000e+03,
    1.0000e+04,
    1.1250e+04,
    1.2500e+04,
    1.4250e+04,
    1.6000e+04,
    1.7000e+04,
    1.8000e+04,
    1.9000e+04,
    2.0000e+04
};
const int    ISO226_0phon_sz        = 63; // length
const int    ISO226_0phon_ref_index = 45; 
const double ISO226_dB_ref          = Eql0PhonDB[ISO226_0phon_ref_index]; 
const double ISO226_freq_ref        = Eql0PhonFreq[ISO226_0phon_ref_index]; 

//###################################################################
ISO226_Interpolator *ISO226_Interpolator::_instance = 
    new ISO226_Interpolator(); 

//###################################################################
// Constructor ISO226_Interpolator 
//###################################################################
ISO226_Interpolator::
ISO226_Interpolator()
{
    _spline.init(ISO226_0phon_sz, &(Eql0PhonFreq[0]), &(Eql0PhonDB[0])); 
}

//###################################################################
// Function Weight
//   Compute frequency weighting based on the following formula
//              W = 10^((ISO226(f_ref)-ISO226(freq))/20)
//   f_ref is the reference frequency at which ISO226 is lowest (most 
//   sensitive). C-Spline is used to interpolate values, and 
//   extrapolation is done by clipping it to boundaries. 
//   @param freq input frequency
//   @return W
//###################################################################
double ISO226_Interpolator::
Weight(const double freq) const
{
    const int &sz = ISO226_0phon_sz;
    if (freq <= Eql0PhonFreq[0   ]) 
    {
        return DBToDecimel(ISO226_dB_ref - Eql0PhonDB[0   ]);
    }
    else if (freq >= Eql0PhonFreq[sz-1]) 
    {
        return DBToDecimel(ISO226_dB_ref - Eql0PhonDB[sz-1]);
    }
    return DBToDecimel(ISO226_dB_ref - _spline.eval(freq));
}

//###################################################################
// Function PrintAllWeights
//###################################################################
void ISO226_Interpolator::
PrintAllWeights(const double &start, const double &stop, 
                const double &interval, std::ostream &stream) const
{
    double freq = start; 
    while (freq < stop)
    {
        stream << freq << " " << Weight(freq) << std::endl;
        freq += interval; 
    }
}

//###################################################################
// Function Test
//###################################################################
void ISO226_Interpolator::
Test() const
{
    std::string file;
    std::ofstream stream; 
    {
        file = "test_weights.txt"; 
        stream.open(file.c_str(), std::ios_base::out|std::ios_base::trunc); 
        PrintAllWeights(0., 44100., 20., stream); 
        stream.close(); 
    }
}

}; // namespace FreqWeighting
