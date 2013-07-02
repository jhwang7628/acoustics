#include "Simulator.h"
#include <algorithm>
#include "rigid/LSCollisionDetect.hpp"
#include "logging/logging.h"

#include <linearalgebra/Vector3.hpp>

RigidBodySimulator::RigidBodySimulator():m_ts(0)
{
    m_constraints.push_back(new GroundCollisionConstraint<REAL, TCollProc> );

#ifdef USE_RECORDER
    //// initialize recorder
    printf("Initializing recorders\n" );
    m_impRec.init("impulses.txt", "modalImpulses.txt");
    printf("Initialized impulse recorder\n");
    m_dispRec.init("displace.bin");
    printf("Done\n" );
#endif
}

#ifdef USE_RECORDER
RigidBodySimulator::RigidBodySimulator(
        const char* recFile, const char *modalRecFile,
        const char* dispFile):m_ts(0)
{
    m_constraints.push_back(new GroundCollisionConstraint<REAL, TCollProc>);
    //// initialize recorder
    m_impRec.init(recFile, modalRecFile);
    m_dispRec.init(dispFile);
}
#endif

RigidBodySimulator::~RigidBodySimulator()
{
    delete m_constraints[0];
}

/*
 * - update external force
 * - collision detection and response, may apply some impulse on rigid body
 * - advance velocity
 * - contact resolution
 * - advance position
 */
void RigidBodySimulator::advance(REAL dt)
{
    const int NRB = m_rigidObjs.size();
#ifdef USE_RECORDER
    //// write the displacement of all the rigid objects into file
    m_dispRec.begin_record(m_ts);
    for(size_t i = 0;i < NRB;++ i)
        m_dispRec.record_displacement(m_rigidObjs[i]);
    m_dispRec.end_record();
#endif
    //// update external force
    for(size_t i = 0;i < NRB;++ i)
    {
        if ( m_rigidObjs[i]->is_fixed() )
        {
          continue;
        }
        m_rigidObjs[i]->apply_forces();
    }

    //// collision response
    resolve_collision(dt);

    //// advance velocity
    for(size_t i = 0;i < NRB;++ i)
        if ( !m_rigidObjs[i]->is_fixed() ) 
            m_rigidObjs[i]->advance_velocity(dt, m_ts);

    //// contact response
    resolve_contact(dt);

    //// advance position
    for(size_t i = 0;i < NRB;++ i)
        if ( !m_rigidObjs[i]->is_fixed() ) 
            m_rigidObjs[i]->advance_state(dt, m_ts);

    m_ts += dt;
}

void RigidBodySimulator::init()
{
    m_shuffledObjIds.resize(m_rigidObjs.size());
    for(size_t i = 0;i < m_shuffledObjIds.size();++ i)
        m_shuffledObjIds[i] = i;

    m_contactGraph.init(m_rigidObjs, m_constraints);
}

/*
 * compute the collision matrix, i.e. v_new = v_old + K*j
 */
void RigidBodySimulator::collision_matrix(const TRigidBody* b, 
        const Vector3<REAL>& r, Matrix3<REAL>& K)
{
    // K = R*I^{-1}*R
    Matrix3<REAL> R(0., -r.z,  r.y,
                   r.z,    0, -r.x,
                  -r.y,  r.x,    0);
    K = R*b->m_predinvI*R;
    // K = -R*I^{-1}*R
    K *= -1;
    // K = 1/M - R*I^{-1}*R
    K.cols[0][0] += b->mass_inverse();
    K.cols[1][1] += b->mass_inverse();
    K.cols[2][2] += b->mass_inverse();
}

/* $ COLLISION
 */
void RigidBodySimulator::resolve_collision(REAL dt)
{
    // increase the timestamp such that m_predx/m_predinvq get updated
    // using the current m_acc
    for(size_t i = 0;i < m_rigidObjs.size();++ i)
        ++ m_rigidObjs[i]->collision_processor()->m_predTimestamp;

    // shuffle the obj id by swap two in the m_shuffledObjIds list
    int i = rand() % m_rigidObjs.size();
    int j = rand() % m_rigidObjs.size();
    if ( i != j ) std::swap(m_shuffledObjIds[i], m_shuffledObjIds[j]);

    int isColliding = 1;
    int collIts;
    for(collIts = 0;isColliding && collIts < COLLISION_DETECT_TOT_ITS;++ collIts)
    {
        isColliding = 0;

        //// collision constraints
        for(i = 0;i < (int)m_shuffledObjIds.size();++ i)
        {
            if ( m_rigidObjs[m_shuffledObjIds[i]]->is_fixed()
              || m_rigidObjs[m_shuffledObjIds[i]]->specified_state())
            {
              continue;
            }

            isColliding |= collision_constraint_resp(
                    m_rigidObjs[m_shuffledObjIds[i]],
                    dt);
        }

        // check each pair of obj to see if they collide
        for(i = 0;i < (int)m_shuffledObjIds.size();++ i)
        for(j = i+1;j < (int)m_shuffledObjIds.size();++ j)
        {
            if ( ( m_rigidObjs[m_shuffledObjIds[i]]->is_fixed() &&
                   m_rigidObjs[m_shuffledObjIds[j]]->is_fixed() )
              || ( m_rigidObjs[m_shuffledObjIds[i]]->specified_state() &&
                   m_rigidObjs[m_shuffledObjIds[j]]->specified_state() ) )
            {
                continue;
            }

            isColliding |= collision_pair_resp(
                    m_rigidObjs[m_shuffledObjIds[i]], 
                    m_rigidObjs[m_shuffledObjIds[j]],
                    dt);
        }
        // if no collision detected, isColliding will be 0
    }
    LOGGING_INFO("%d iterations for collision detection [ts=%lf]", collIts, m_ts);
}

/* $ COLLISION
 * resolve collisions for the given pair of objects
 */
int RigidBodySimulator::collision_pair_resp(TRigidBody* b1, TRigidBody* b2, REAL dt)
{
    CollisionRec<REAL> cRecb12, cRecb21;
    TCollProc* cp1 = b1->collision_processor();
    TCollProc* cp2 = b2->collision_processor();

    for(int collIts = 0;collIts < COLLISION_DETECT_PAIR_ITS;++ collIts)
    {
        //// move the object to the predicted position
        b1->update_force_predicted_position(dt);
        b2->update_force_predicted_position(dt);

        cp1->m_collidedVtx.clear();
        cp2->m_collidedVtx.clear();

        //// get the candidate vertices for collision detection
        detect_tree_node_collisions(
                cp1->bounding_tree_root(),
                cp2->bounding_tree_root(),
                cp1, cp2);

        //// find deepest vertex of b1 that is inside b2, and the deepest vertex
        //   of b2 that is inside b1
        bool cb12 = cp1->deepest_penetrating_vertex(b2, cRecb12);
        bool cb21 = cp2->deepest_penetrating_vertex(b1, cRecb21);
        if ( !cb12 && !cb21 )
            return collIts > 0;
        else // collision detected 
        {
            if ( cb12 && (!cb21 || cRecb12.depth < cRecb21.depth) )
                apply_collision_impulse(b1, b2, cRecb12);   // cRecb12.vtxId refers to the vertex in b1
            else
                apply_collision_impulse(b2, b1, cRecb21);   // cRecb21.vtxId refers to the vertex in b2
        }
    }
    return 1;
}

/* $ COLLISION
 * apply collision constraints to the given body.
 * collision constraints is the constraints imposed by collisions like wall, ground, etc
 */
int RigidBodySimulator::collision_constraint_resp(TRigidBody* b, REAL dt)
{
    CollisionRec<REAL> cRec;
    TCollProc* cp = b->collision_processor();

    for(int collIts = 0;collIts < COLLISION_DETECT_PAIR_ITS;++ collIts)
    {
        b->update_force_predicted_position(dt);

        if ( !cp->detect_collision_constraint(m_constraints, cRec) )
            return collIts > 0;
        else
            apply_collision_impulse(b, cRec);
    }
    return 1;
}

/* $ COLLISION
 * The collison impulse of a object with a fixed object, like the ground, walls, etc
 * Collision response with friction impulse considered.
 */
void RigidBodySimulator::apply_collision_impulse(
        TRigidBody* b, const CollisionRec<REAL>& cRec)
{
    //// assuming stick friction, compute the impulse
    const Vector3<REAL> r = cRec.pt - b->m_predx;
    Matrix3<REAL> K;
    collision_matrix(b, r, K);   // K should be a 3x3 matrix
    Vector3<REAL> Jt = K.inverse() * (cRec.impulseDir * 
            (-cRec.eps*cRec.vnrel) - cRec.vrel);
    const REAL JnSqr = M_SQR(Jt.dotProduct(cRec.impulseDir));  // square of normal component of the impulse
    const REAL normJSqr = Jt.lengthSqr();

#ifdef USE_RECORDER
    //m_impRec.record_impulse(m_ts, Jt, cRec.vtxId, b, true);
    m_impRec.record_constraint_impulse(m_ts, Jt, cRec, b);
#endif

    //// check if the impulse is in the friction cone.
    if ( normJSqr - JnSqr > M_SQR(cRec.mu)*JnSqr )     // |J_t|^2 <= \mu^2|J_n|^2 (sticking friction)
    {       // sliding friction
        Vector3<REAL> T = cRec.vrel - cRec.vnrel*cRec.impulseDir;   // tangential dir
        T.normalize();
        Vector3<REAL> td = cRec.impulseDir;
        td.scaleAdd(-cRec.mu, T);                       // td = N - \mu*T

        REAL jn = -(1. + cRec.eps) * cRec.vnrel;
        jn /= cRec.impulseDir.dotProduct(K*td);

        // J = jn*N - \mu*jn*T
        Jt = cRec.impulseDir * jn;
        Jt.scaleAdd(-cRec.mu*jn, T);
    }

    b->apply_impulse_to_prediction(Jt, r);
    ++ b->collision_processor()->m_predTimestamp;
}

/* $ COLLISION
 * NOTE: the cRec.impulseDir is pointing outward of bb object
 *       and cRec.vrel is the relative velocity (v_a - v_b) 
 * NOTE: cRec.vtxId always refers to the vertex where the impulse 
 *       is applied in the rigid body \ba
 */
void RigidBodySimulator::apply_collision_impulse(
        TRigidBody* ba, TRigidBody* bb, const CollisionRec<REAL>& cRec)
{
    const Vector3<REAL> ra = cRec.pt - ba->m_predx;
    const Vector3<REAL> rb = cRec.pt - bb->m_predx;
    Matrix3<REAL> Ka, Kb, K;

    if ( !ba->is_fixed() )
    {
        collision_matrix(ba, ra, Ka);
        K += Ka;
    }

    if ( !bb->is_fixed() )
    {
        collision_matrix(bb, rb, Kb);
        K += Kb;
    }

    Vector3<REAL> Jt = K.inverse() * (cRec.impulseDir * 
            (-cRec.eps*cRec.vnrel) - cRec.vrel);
    const REAL JnSqr = M_SQR(Jt.dotProduct(cRec.impulseDir));
    const REAL normJSqr = Jt.lengthSqr();

#ifdef USE_RECORDER
    //m_impRec.record_impulse(m_ts,  Jt, cRec.vtxId, ba, true);
    //m_impRec.record_impulse(m_ts, -Jt, cRec.pt, bb);
    m_impRec.record_inter_obj_impulse(m_ts, Jt, cRec, ba, bb);
#endif

    if ( normJSqr - JnSqr > M_SQR(cRec.mu)*JnSqr )
    {   // sliding friction
        Vector3<REAL> T = cRec.vrel - cRec.vnrel*cRec.impulseDir;   // tangential dir
        T.normalize();
        Vector3<REAL> td = cRec.impulseDir;
        td.scaleAdd(-cRec.mu, T);                       // td = N - \mu*T

        REAL jn = -(1. + cRec.eps) * cRec.vnrel;
        jn /= cRec.impulseDir.dotProduct(K*td);

        Jt = cRec.impulseDir * jn;
        Jt.scaleAdd(-cRec.mu*jn, T);
    }

    ba->apply_impulse_to_prediction( Jt, ra);
    bb->apply_impulse_to_prediction(-Jt, rb);

    ++ ba->collision_processor()->m_predTimestamp;
    ++ bb->collision_processor()->m_predTimestamp;
}

// ------------------------------------------------------------------------------------
/* $ CONTACT
 */
void RigidBodySimulator::resolve_contact(REAL dt)
{
    int i, j;
    //// create the contact graph
    // label the m_predTimestamp because m_v has been changed by advance_velocity
    for(i = 0;i < (int)m_rigidObjs.size();++ i)
    {
        ++ m_rigidObjs[i]->collision_processor()->m_predTimestamp;
        m_rigidObjs[i]->update_velocity_predicted_position(dt);
    }
    m_contactGraph.update();
    //m_contactGraph.print_graph();

    int isColliding = 1;
    int collIts;
    const std::vector<int>& sortedObj = m_contactGraph.sorted_objs();

    for(collIts = 0;isColliding && collIts < CONTACT_DETECT_TOT_ITS;++ collIts)
    {
        isColliding = 0;

        //// collision constraints
        for(i = 0;i < (int)sortedObj.size();++ i)
        {
            // fixed obj is assumed not to collide with any constraints
            if ( m_rigidObjs[sortedObj[i]]->is_fixed() 
              || m_rigidObjs[sortedObj[i]]->specified_state() ) continue;

            isColliding |= contact_constraint_resp(
                    m_rigidObjs[sortedObj[i]], dt, collIts);
        }

        // check each pair of obj to see if they collide
        for(i = 1;i < (int)sortedObj.size();++ i)
        for(j = 0;j < i;++ j)
        {
            if ( ( m_rigidObjs[sortedObj[i]]->is_fixed() && 
                   m_rigidObjs[sortedObj[j]]->is_fixed() ) 
              || ( m_rigidObjs[sortedObj[i]]->specified_state() &&
                   m_rigidObjs[sortedObj[j]]->specified_state() ) ) continue;

            isColliding |= contact_pair_resp(
                    m_rigidObjs[sortedObj[i]], m_rigidObjs[sortedObj[j]],
                    dt, collIts);
        }
        // if no collision detected, isColliding will be 0
    }

    LOGGING_INFO("%d iterations for contact detection [ts=%lf]", collIts, m_ts);

    if ( isColliding )  // last iteration using shock propagation
    {
        const std::vector<int>& levels = m_contactGraph.levels();
        for(i = 0;i < (int)sortedObj.size();++ i)
        {
            if ( levels[sortedObj[i]] < 0 )         // fixed obj
            {
                m_rigidObjs[sortedObj[i]]->m_pinned = true;
                continue;
            }

            isColliding = 1;
            for(int ccnt2 = 0;isColliding && ccnt2 < CONTACT_SHOCK_PROP_ITS;++ ccnt2)
            {
                isColliding = 0;
                for(j = 0;j < i;++ j)
                    isColliding |= contact_pair_resp(
                        m_rigidObjs[sortedObj[i]], m_rigidObjs[sortedObj[j]],
                        dt, collIts + ccnt2);

                isColliding |= contact_constraint_resp(
                        m_rigidObjs[sortedObj[i]], dt, collIts + ccnt2);
            }
            m_rigidObjs[sortedObj[i]]->m_pinned = true;
        }

        for(i = 0;i < (int)m_rigidObjs.size();++ i) m_rigidObjs[i]->m_pinned = false;
    }
}

/* $ CONTACT
 * apply contact constraints to the given body.
 * contact constraints is the constraints imposed by collisions like wall, ground, etc
 */
int RigidBodySimulator::contact_constraint_resp(TRigidBody* b, REAL dt, int numIt)
{
    CollisionRec<REAL> cRec;
    TCollProc* cp = b->collision_processor();

    for(int collIts = 0;collIts < CONTACT_DETECT_PAIR_ITS;++ collIts)
    {
        //// predict the position based on current velocity
        b->update_velocity_predicted_position(dt);

        //// check if there is still collisions with those constrains
        if ( !cp->detect_collision_constraint(m_constraints, cRec) )
            return collIts > 0;
        else
            apply_contact_impulse(b, cRec, contact_rest_coeff(numIt));
    }
    return 1;
}

/* $ CONTACT
 * resolve contact for each pair
 */
int RigidBodySimulator::contact_pair_resp(TRigidBody* b1, TRigidBody* b2, REAL dt, int numIt)
{
    CollisionRec<REAL> cRecb12, cRecb21;
    TCollProc* cp1 = b1->collision_processor();
    TCollProc* cp2 = b2->collision_processor();

    for(int collIts = 0;collIts < CONTACT_DETECT_PAIR_ITS;++ collIts)
    {
        //// move the object to the predicted position
        b1->update_velocity_predicted_position(dt);
        b2->update_velocity_predicted_position(dt);

        cp1->m_collidedVtx.clear();
        cp2->m_collidedVtx.clear();

        //// get the candidate vertices for collision detection
        detect_tree_node_collisions(
                cp1->bounding_tree_root(),
                cp2->bounding_tree_root(),
                cp1, cp2);

        //// find deepest vertex of b1 that is inside b2, and the deepest vertex
        //   of b2 that is inside b1
        bool cb12 = cp1->deepest_penetrating_vertex(b2, cRecb12);
        bool cb21 = cp2->deepest_penetrating_vertex(b1, cRecb21);
        if ( !cb12 && !cb21 )
            return collIts > 0;
        else // collision detected 
        {
            if ( cb12 && (!cb21 || cRecb12.depth < cRecb21.depth) )
                apply_contact_impulse(b1, b2, cRecb12, contact_rest_coeff(numIt));
            else
                apply_contact_impulse(b2, b1, cRecb21, contact_rest_coeff(numIt));
        }
    }
    return 1;
}

/* 
 * This method should never be called if object is pinned or fixed
 * No friction for contact impulse 
 */
void RigidBodySimulator::apply_contact_impulse(
        TRigidBody* b, const CollisionRec<REAL>& cRec,
        const REAL restC)
{
    REAL j = -(1. + restC) * cRec.vnrel;
    const Vector3<REAL> r = cRec.pt - b->m_predx;

    REAL term3 = cRec.impulseDir.dotProduct((
                b->m_predinvI * (r.crossProduct(cRec.impulseDir))).crossProduct(r));
    j /= (b->mass_inverse() + term3);

#ifdef USE_RECORDER
    const Vector3<REAL> Jt = j*cRec.impulseDir;
    b->apply_impulse_to_prediction(Jt, r);
#warning ------ NOTE HERE ------
    //m_impRec.record_impulse(m_ts, Jt, cRec.vtxId, b, true);
#else
    b->apply_impulse_to_prediction(j*cRec.impulseDir, r);
#endif
    ++ b->collision_processor()->m_predTimestamp;
}

/* 
 * Apply contact impulse to objects b1 and b2
 *
 * No friction for contact impulse 
 */
void RigidBodySimulator::apply_contact_impulse(
        TRigidBody* b1, TRigidBody* b2, 
        const CollisionRec<REAL>& cRec, const REAL restC)
{
    const Vector3<REAL> ra = cRec.pt - b1->m_predx;
    const Vector3<REAL> rb = cRec.pt - b2->m_predx;
    const bool dob1 = !b1->m_pinned && !b1->is_fixed()
                        && !b1->specified_state(); // if obj1 is fixed
    const bool dob2 = !b2->m_pinned && !b2->is_fixed()
                        && !b2->specified_state(); // if obj2 is fixed

    REAL j = -(1. + restC) * cRec.vnrel;
    REAL denom = 0.;
    if ( dob1 )
    {
        denom += b1->mass_inverse();
        denom += cRec.impulseDir.dotProduct((
                b1->m_predinvI * (ra.crossProduct(cRec.impulseDir))).crossProduct(ra));
    }

    if ( dob2 )
    {
        denom += b2->mass_inverse();
        denom += cRec.impulseDir.dotProduct((
                b2->m_predinvI * (rb.crossProduct(cRec.impulseDir))).crossProduct(rb));
    }

    assert(denom > 1E-9);
    j /= denom;
    Vector3<REAL> impulse = j * cRec.impulseDir;

    if ( dob1 )
    {
        b1->apply_impulse_to_prediction( impulse, ra);
        ++ b1->collision_processor()->m_predTimestamp;
    }
    if ( dob2 )
    {
        b2->apply_impulse_to_prediction(-impulse, rb);
        ++ b2->collision_processor()->m_predTimestamp;
    }
}

