#ifndef MATHUTILS_RK4_HPP
#define MATHUTILS_RK4_HPP

namespace MathUtils
{

/**
 * The standard Runge-Kutta 4 integrator.
 */
template<class State>
class RK4
{
    public:
        typedef typename State::Scalar Scalar;

        template<typename Evaluator>
        State
        step(const State &currState,
             Scalar t,
             Scalar dt,
             Evaluator &func)
        {
            State k1 = func(currState,
                            t);

            State k2 = func(currState + dt/2. * k1,
                            t + dt/2.);

            State k3 = func(currState + dt/2. * k2,
                            t + dt/2.);

            State k4 = func(currState + dt * k3,
                            t + dt);

            return currState + dt/6. * (k1 + 2. * k2 + 2. * k3 + k4);
        }
};

} // namespace MathUtils

#endif // MATHUTILS_RK4_HPP

