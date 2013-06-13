/*
 * =====================================================================================
 *
 *       Filename:  timer.h
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

template <bool NanoAccu>
class Timer
{
    public:
        Timer():tStart_(0.), elapsed_(0.) { }

        void start();
        void pause();

        inline void reset()
        {   elapsed_ = tStart_ = 0.; }

        inline double elapsed() const
        {   return elapsed_; }

    private:
        double  tStart_;
        double  elapsed_;
};

template<>
void Timer<false>::start()
{   tStart_ = GetMilliTimed(); }

template<>
void Timer<false>::pause()
{   elapsed_ += (GetMilliTimed()-tStart_); }

template<>
void Timer<true>::start()
{   tStart_ = GetNanoTimed(); }

template<>
void Timer<true>::pause()
{   elapsed_ += (GetNanoTimed()-tStart_); }

#endif   /* ----- #ifndef UTILS_TIMER_INC  ----- */

