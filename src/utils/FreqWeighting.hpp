#ifndef _FREQ_WEIGHTING_HPP
#define _FREQ_WEIGHTING_HPP

#include <iostream>
#include <utils/Exception.h>

/**
 * Utilities for weighting frequency curves
 */

namespace FreqWeighting
{
    const double ITURFreqs[21] = {31.5,
                                  63,
                                  100,
                                  200,
                                  400,
                                  800,
                                  1000,
                                  2000,
                                  3150,
                                  4000,
                                  5000,
                                  6300,
                                  7100,
                                  8000,
                                  9000,
                                  10000,
                                  12500,
                                  14000,
                                  16000,
                                  20000,
                                  31500};

    const double ITURDbGains[21] = {-29.9,
                                    -23.9,
                                    -19.8,
                                    -13.8,
                                    -7.8,
                                    -1.9,
                                    0,
                                    5.6,
                                    9,
                                    10.5,
                                    11.7,
                                    12.2,
                                    12,
                                    11.4,
                                    10.1,
                                    8.1,
                                    0,
                                    -5.3,
                                    -11.7,
                                    -22.2,
                                    -42.7};

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

    /**
     * ITU-R 468 weighting
     * Takes a sound level in decibels at a specific
     * frequency, and returns the ITU-R weighted DB value.
     */
    double
    ITURWeight(double freq,
               double db)
    {
        // Simple linear interpolation for now
        if (freq <= ITURFreqs[0])
        {
            return db + ITURDbGains[0];
        }

        if (freq >= ITURFreqs[20])
        {
            return db + ITURDbGains[20];
        }

        // Find correct range
        int i = 0;
        while (i < 21 && freq >= ITURFreqs[i])
        {
            ++i;
        }

        if (i == 0 || i > 20)
        {
            std::cerr << "invalid range: " << i << std::endl;
            throw EXCEPTION("invalid range");
        }

        double t = (freq - ITURFreqs[i-1]) / (ITURFreqs[i] - ITURFreqs[i-1]);

        return db + ((1-t) * ITURDbGains[i-1] + t * ITURDbGains[i]);
    }

    /**
     * ISO 226 weighting
     * Takes a sound level in decibels at a specific
     * frequency, and returns the (inverse) ISO226 weighted DB value.
     */
    double
    ISO226Weight(double freq,
               double db)
    {
        const int sz = 63;
        // Simple linear interpolation for now
        if (freq <= Eql0PhonFreq[0])
        {
            return db - Eql0PhonDB[0];
        }

        if (freq >= Eql0PhonFreq[sz-1])
        {
            return db - Eql0PhonDB[sz-1];
        }

        // Find correct range
        int i = 0;
        while (i < sz && freq >= Eql0PhonFreq[i])
        {
            ++i;
        }

        if (i == 0 || i > sz-1)
        {
            std::cerr << "invalid range: " << i << std::endl;
            throw EXCEPTION("invalid range");
        }

        double t = (freq - Eql0PhonFreq[i-1]) / (Eql0PhonFreq[i] - Eql0PhonFreq[i-1]);
        std::cout << "Freq: " << freq << ", i: " << i << ", t: " << t << std::endl;

        return db - ((1-t) * Eql0PhonDB[i-1] + t * Eql0PhonDB[i]);
    }
};

#endif // _FREQ_WEIGHTING_HPP

