#ifndef MATHUTILS_FORWARD_EULER_HPP
#define MATHUTILS_FORWARD_EULER_HPP

namespace MathUtils
{

/**
 * A forward Euler integrator. This is a simple but
 * bad integrator.
 */
template<typename State>
class ForwardEuler
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
            State deriv = func(currState,
                               t);

            return currState + dt * deriv;
        }
};

} // namespace MathUtils

#endif // MATHUTILS_FORWARD_EULER_HPP

