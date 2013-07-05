/*
 * =====================================================================================
 *
 *       Filename:  timer.cpp
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

#include "timer.hpp"

template<>
void Timer<false>::start()
{
    tStart_ = GetMilliTimed();
}

template<>
void Timer<false>::pause()
{
    elapsed_ += (GetMilliTimed()-tStart_);
    cycles_ += 1;
}

template<>
void Timer<true>::start()
{
    tStart_ = GetNanoTimed();
}

template<>
void Timer<true>::pause()
{
    elapsed_ += (GetNanoTimed()-tStart_);
    cycles_ += 1;
}

template<>
double Timer<false>::getMsPerCycle() const
{
    return ( cycles_ > 0 ) ? ( elapsed_ / (double)cycles_ ) : 0.0;
}

template<>
double Timer<true>::getMsPerCycle() const
{
    return ( cycles_ > 0 ) ? ( elapsed_ / (double)cycles_ / 1000000.0 ) : 0.0;
}
