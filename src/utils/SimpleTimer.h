#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H
#include <iostream>
#include <chrono>

//##############################################################################
// Class SimpleTimer
//##############################################################################
class SimpleTimer
{
    private:
        std::chrono::high_resolution_clock::time_point  _start;
        std::chrono::high_resolution_clock::duration    _elapsed;

    public:
        SimpleTimer() : _elapsed(0) {}

        inline void Start() { _start = std::chrono::high_resolution_clock::now(); }
        inline void Pause() { _elapsed += std::chrono::high_resolution_clock::now() - _start;}
        inline double Duration() { return (std::chrono::duration_cast<std::chrono::duration<double>>(_elapsed)).count(); }

        //// debug methods ////
        inline void DurationDebug()
        {
            std::cout << "elapsed(native) = " << _elapsed.count() << std::endl;
            std::cout << "elapsed(double) = " << (std::chrono::duration_cast<std::chrono::duration<double>>(_elapsed)).count() << std::endl;
        }
};
//##############################################################################
#endif
