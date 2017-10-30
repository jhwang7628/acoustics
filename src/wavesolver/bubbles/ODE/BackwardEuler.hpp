#ifndef MATHUTILS_BACKWARD_EULER_HPP
#define MATHUTILS_BACKWARD_EULER_HPP

namespace MathUtils
{

/**
 * A backward Euler integrator.
 * Find state_{n+1} such that
 * state_{n+1} = state_n + dt * f(t_{n+1}, state_{n+1})
 */
template<typename State, typename Jacobian>
class BackwardEuler
{
    public:
        typedef typename State::Scalar Scalar;

        BackwardEuler(Scalar tolerance = 1.e-6)
            : m_tolerance(tolerance)
        {
        }

        template<typename Evaluator>
        State
        step(const State &currState,
             Scalar t,
             Scalar dt,
             const Evaluator &func)
        {
            State deriv = dt * func(currState,
                                    t + dt);

            Jacobian jac = dt * func.jacobian(currState,
                                              t + dt);

            jac -= Jacobian::Identity();

            State deltaState = jac.colPivHouseholderQr().solve(deriv);

            State newState = currState - deltaState;

            // TODO: should tolerance be based on step size (deltaState)
            // or on the residual newState - currState - dt * deriv?
            while ( deltaState.norm() > m_tolerance )
            {
                // Take another step
                deriv = dt * func(newState,
                                  t + dt);

                jac = dt * func.jacobian(newState,
                                         t + dt);

                jac -= Jacobian::Identity();

                deltaState = jac.colPivHouseholderQr().solve(currState - newState + deriv);

                newState -= deltaState;
            }

            return newState;
        }

    private:
        Scalar m_tolerance;
};

} // namespace MathUtils

#endif // MATHUTILS_BACKWARD_EULER_HPP

