#ifndef MATHUTILS_MIDPOINT_METHOD_HPP
#define MATHUTILS_MIDPOINT_METHOD_HPP

namespace MathUtils
{

/**
 * Midpoint method integrator.
 */
template<class State>
class MidpointMethod
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
            State currDeriv = func(currState,
                                   t);

            State deriv = func(currState + dt/2. * currDeriv,
                               t + dt/2.);

            return currState + dt * deriv;
        }
};

} // namespace MathUtils

#endif // MATHUTILS_MIDPOINT_METHOD_HPP

