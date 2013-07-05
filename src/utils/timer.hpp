/*
 * =====================================================================================
 *
 *       Filename:  timer.hpp
 *
 *    Description:  Timer class
 *
 *        Version:  1.0
 *        Created:  12/13/10 00:16:21
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  UTILS_TIMER_INC
#define  UTILS_TIMER_INC

#include "nano_timer.h"

#include <string>

template <bool NanoAccu>
class Timer
{
    public:
        Timer()
            : tStart_(0.),
              elapsed_(0.),
              cycles_(0)
        {
        }

        Timer(const std::string &name)
            : tStart_(0.),
              elapsed_(0.),
              cycles_(0),
              name_(name)
        {
        }

        void start();
        void pause();

        inline void reset()
        {
            elapsed_ = tStart_ = 0.;
            cycles_ = 0;
        }

        inline double elapsed() const
        {   return elapsed_; }

        double getMsPerCycle() const;

    private:
        double       tStart_;
        double       elapsed_;

        unsigned int cycles_;

        std::string  name_;
};

#endif   /* ----- #ifndef UTILS_TIMER_INC  ----- */

