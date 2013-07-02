#ifndef RIGID_COLLISION_REC_HPP
#   define RIGID_COLLISION_REC_HPP

#include "geometry/Point3.hpp"
#include "linearalgebra/Vector3.hpp"

/*!
 * CollisionRec is associated with each rigid body to record the information about collisions.
 * The information includes:
 * - the penetration depth at each vertex
 * - direction of collision impulse at the vertex
 */
template <typename T>
struct CollisionRec
{
    T           depth;      // negative value
    Point3<T>   pt;         // the collision point in predicted configuration
    Vector3<T>  impulseDir;
    Vector3<T>  vrel;
    T           vnrel;
    T           eps;        // restitution coefficient
    T           mu;         // friction coefficient
#ifdef USE_RECORDER
    int         vtxId;      // the vtx id in surface triangle mesh. the referred 
                            // vertex is the current deepest vertx that penetrates
                            // into the other objects
#endif
};

#endif
