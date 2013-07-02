#ifndef RIGID_BODY_HPP
#   define RIGID_BODY_HPP

#include "config.h"
#include <stdlib.h>
#include "geometry/Point3.hpp"
#include "linearalgebra/Vector3.hpp"
#include "linearalgebra/Quaternion.hpp"

#include <boost/function.hpp>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

typedef boost::function<double (double t)>   StateEvaluator;

struct RigidState {
  StateEvaluator             positionSignal[ 3 ];
  StateEvaluator             velocitySignal[ 3 ];
};

/*!
 * The basic rigid body
 *
 * \param TMesh is the tetrahedron mesh
 */
template <typename T, class TMesh>
class RigidBody
{
    public:
        RigidBody(TMesh* mesh, T density):m_fixed(false), 
                m_density(density), mp_mesh(mesh), mp_stateSignal(NULL)
        {
            init();
        }

        /*
         * This method should only be called when an object 
         * is unbreakable, in which case we don't need the 
         * tet mesh for stress analysis
         */
        TMesh* clean_tet_mesh()
        {
            TMesh* ret = mp_mesh;
            mp_mesh = NULL;
            return ret;
        }

        /*
         * given the initial position of a point on this rigid body,
         * return the current position of that point due to the movement 
         * of this rigid body
         */
        inline Point3<T> current_position(const Point3<T>& pt) const
        {   return m_x + m_q.rotate(pt - m_x0);  }
        
        void advance_velocity(T dt, T t = 0.0)
        {
            if ( specified_state() )
            {
              m_v[ 0 ] = mp_stateSignal->velocitySignal[ 0 ]( t );
              m_v[ 1 ] = mp_stateSignal->velocitySignal[ 1 ]( t );
              m_v[ 2 ] = mp_stateSignal->velocitySignal[ 2 ]( t );
              
              m_omega.zero();
            }
            else
            {
              m_v.scaleAdd(dt, m_acc);
              m_omega.scaleAdd(dt, m_angAcc);     // advance angular velocity
            }
        }

        void advance_state(T dt, double t = 0.0)
        {
            if ( specified_state() )
            {
              m_x[ 0 ] = mp_stateSignal->positionSignal[ 0 ]( t );
              m_x[ 1 ] = mp_stateSignal->positionSignal[ 1 ]( t );
              m_x[ 2 ] = mp_stateSignal->positionSignal[ 2 ]( t );
            }
            else
            {
              m_x.scaleAdd(dt, m_v);

              // the rotation in dt is [cos(|omega|*dt*0.5),
              // sin(|omega|*dt*0.5)*omega/|omega|)]
              T omeganorm = m_omega.length();
              if ( omeganorm > EPS )
              {
                  const T dang = omeganorm * dt * 0.5;
                  Quaternion<T> rot(cos(dang), m_omega * (sin(dang) / omeganorm));
                  m_q = rot * m_q;
                  m_q.normalize();
              }
              m_invq = m_q.conjugate();

              // m_invI = R * m_invI0 * R^T
              const Matrix3<T> R = m_q.toMatrix3();
              m_invI = R * m_invI0 * R.transpose();
            }
        }

        inline void set_velocity(const Vector3<T>& v, const Vector3<T>& angV)
        {
            m_v = v;
            m_omega = angV;
        }

        inline void add_impulses(const Vector3<T>& J, const Vector3<T>& L)
        {   
            m_v += (J*m_invMass);
            m_omega += (m_invI*L);
        }

        void set_position(const Point3<T>& p, const Quaternion<T>& rot)
        {
            m_x = p;
            m_q = rot;
            m_invq = m_q.conjugate();
            const Matrix3<T> R = m_q.toMatrix3();
            m_invI = R * m_invI0 * R.transpose();
        }

        void translate(T dx, T dy, T dz)
        {
            m_x.x += dx;
            m_x.y += dy;
            m_x.z += dz;
        }

        void init_origin_rotation(T rw, T rx, T ry, T rz)
        {
            m_q = Quaternion<T>(rw, rx, ry, rz);
            m_q.normalize();
            m_invq = m_q.conjugate();

            m_x -= m_x0;
            Vector3<T> xx = m_q.rotate(m_x0);
            m_x += xx;
        }

        void rotate(const Quaternion<T>& rot)
        {
            m_q = rot * m_q;
            m_q.normalize();
            m_invq = m_q.conjugate();
            const Matrix3<T> R = m_q.toMatrix3();
            m_invI = R * m_invI0 * R.transpose();
        }

        /*!
         * Initialize velocity
         *
         * This method simply sets the current velocity. It is often used for 
         * setting the initial velocity at the beginning of the simulation.
         */
        inline void init_velocity(T vx, T vy, T vz)
        {   m_v.set(vx, vy, vz); }

        /*
         * given a vector in the object's initial space, transform it
         * into the current space
         */
        inline void to_current_space(Vector3<T>& vec) const
        {   m_q.rotate_vector(vec); }

        /*
         * Return the velocity of a given point on this rigid body
         * pt is in the current coordinate space
         */
        inline Vector3<T> current_velocity(const Point3<T>& pt) const
        {
            return m_v + m_omega.crossProduct(pt - m_x);
        }

        /*
         * given a vector in the object's current space, transform it
         * back to the initial space
         */
        inline void to_init_space(Vector3<T>& vec) const
        {   m_invq.rotate_vector(vec); }

        inline T mass_inverse() const
        {   return m_invMass; }
        const Matrix3<T>& current_inertia_inverse() const
        {   return m_invI; }
        const Matrix3<T>& init_inertia_inverse() const
        {   return m_invI0; }
        const Point3<T>& initial_mass_center() const
        {  return m_x0; }
        const Point3<T>& mass_center() const
        {  return m_x; }
        const Quaternion<T>& rotation() const
        {  return m_q; }
        T density() const
        {  return m_density; }
        T mass() const
        {  return m_mass; }
        const Vector3<T>& velocity() const
        {  return m_v; }
        const Vector3<T>& angular_velocity() const
        {  return m_omega; }
        bool is_fixed() const
        {  return m_fixed; }
        bool specified_state() const
        {  return ( mp_stateSignal != NULL ); }
        void set_fixed(bool f)
        {  m_fixed = f; }
        const TMesh* mesh() const
        {  return mp_mesh; }
        TMesh* mesh() 
        {  return mp_mesh; }

        void specify_state( RigidState *stateSignal )
        {  mp_stateSignal = stateSignal; }

        T trans_energy() const
        {  return 0.5*m_mass*m_v.lengthSqr(); }
        T rot_energy() const
        {  return 0.5*m_omega.dotProduct(m_invI.inverse()*m_omega); }
        T kinetic_energy() const
        {  return trans_energy() + rot_energy(); }

    protected:
        void init();

    protected:
        bool            m_fixed;
        T               m_density;

        /*======= constant quantities =======*/
        T               m_mass;
        T               m_invMass;
        Point3<T>       m_x0;       // center of mass
        Matrix3<T>      m_invI0;    // inverse of inertia at initial configuration

        /*======= State Variables ======*/
        Point3<T>       m_x;        // current position
        Quaternion<T>   m_q;        // rotation
        Quaternion<T>   m_invq;     // inverse rotation: rotate from current state to initial state
        Matrix3<T>      m_invI;     // inverse of current inertia

        /*======= Derived quantities =======*/
        Vector3<T>      m_v;        // current velocity
        Vector3<T>      m_omega;    // current angular velocity

        /*======= Computed quantities =======*/
        Vector3<T>      m_acc;      // acceleration
        Vector3<T>      m_angAcc;   // angular acceleration

        TMesh*          mp_mesh;

        RigidState     *mp_stateSignal;
};

//////////////////////////////////////////////////////////////////////////////////////

template <typename T, class TMesh>
void RigidBody<T, TMesh>::init()
{
    const std::vector<T>&           ms  = mp_mesh->masses();
    const std::vector< Point3<T> >& vtx = mp_mesh->rest_positions();

    Matrix3<T> I0;  // inertia
    m_mass = 0;
    m_x0.zero();
    // TODO: parallel reduction
    for(size_t i = 0;i < ms.size();++ i)
    {
        m_mass += ms[i];
        m_x0 += ms[i] * vtx[i];
    }

    if ( m_mass < EPS )
    {
        fprintf(stderr, "ERROR: object has zero mass\n");
        exit(1);
    }

    m_x0 /= m_mass;
    m_mass *= m_density;
    m_invMass = 1. / m_mass;

    m_x = m_x0;

    // update inertia
    for(size_t i = 0;i < ms.size();++ i)
    {
        Vector3<T> ri = vtx[i] - m_x0;
        I0 += Matrix3<T>(
                ms[i]*(M_SQR(ri.y)+M_SQR(ri.z)), -ms[i]*ri.x*ri.y, -ms[i]*ri.x*ri.z,
                -ms[i]*ri.y*ri.x, ms[i]*(M_SQR(ri.x)+M_SQR(ri.z)), -ms[i]*ri.y*ri.z,
                -ms[i]*ri.z*ri.x, -ms[i]*ri.z*ri.y, ms[i]*(M_SQR(ri.x)+M_SQR(ri.y)));
    }
    I0 *= m_density;
    m_invI0 = I0.inverse();
    m_invI  = m_invI0;
}

#ifdef USE_NAMESPACE
}
#endif
#endif
