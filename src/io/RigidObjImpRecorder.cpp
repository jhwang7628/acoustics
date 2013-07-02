#include "RigidObjImpRecorder.h"
#include <iomanip>

using namespace std;

void RigidObjImpRecorder::init(const char* file, const char *modalFile,
                               int precision)
{
    printf("In RigidObjImpRecorder::init\n");

    if ( m_fout.is_open() ) {
        m_fout.close();
    }
    m_fout.open(file);
    m_fout << setprecision(precision);

    if ( m_modalImpulseOut.is_open() ) {
        m_modalImpulseOut.close();
    }
    m_modalImpulseOut.open(modalFile);
    m_modalImpulseOut << setprecision(precision);

    printf("Done with RigidObjImpRecorder::init\n");
}

/*
 * - Transform the impulse vector from the current prediction to
 *   the object's initial(rest) configuration
 *
 * the last field in each line is a char, either T or S, indicating
 * if the given vertex id is the surface vtx id or tet mesh vtx id.
 * 'S' means surface vtx id, 'T' means tet id
 * The recording in the file is gonna to be
 * <time>    <obj id>    <applied vtx id>    <impulse>  <T/S>
 * ...
 */
/*
void RigidObjImpRecorder::record_impulse(
        REAL ts, const Vector3<REAL>& imp, int vtxId,
        const TRigidBody* body, bool surfVtx)
{
    // map the impulse vector to object's rest configuration
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);

    m_fout << ts << ' ' << body->id() << ' ' 
           << vtxId << ' '      // vtxId is 0-based
           << impVec.x << ' '
           << impVec.y << ' '
           << impVec.z << ' ' 
           << (surfVtx ? 'S' : 'T') << std::endl;
}
*/

void RigidObjImpRecorder::record_constraint_impulse(
                                      REAL ts, const Vector3<REAL>& imp,
                                      const CollisionRec<REAL>& cRec,
                                      const TRigidBody* body)
{
  {
    // map the impulse vector to object's rest configuration
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);
    /*
    const Vector3<REAL> relVel = body->predicted_inverse_rotation().rotate(cRec.vrel);

    m_modalImpulseOut
            << ts << ' ' << body->id() << ' ' << cRec.vtxId  << ' '
            << impVec.x << ' ' << impVec.y << ' ' << impVec.z << ' '
            <<-relVel.x << ' ' <<-relVel.y << ' ' <<-relVel.z << " S" << endl;
    */
    const Point3d& pos = body->predicted_mass_center();
    const Quat4d&    q = body->predicted_rotation();
#if 0
    m_modalImpulseOut
            << ts << " C " << body->id() << ' ' << cRec.vtxId << ' '
            << impVec.x << ' ' << impVec.y << ' ' << impVec.z << ' '
            << cRec.vrel.x << ' ' << cRec.vrel.y << ' ' << cRec.vrel.z << ' '
            << pos.x << ' ' << pos.y << ' ' << pos.z << ' ' 
            << q.w << ' ' << q.v.x << ' ' << q.v.y << ' ' << q.v.z << endl;
#endif
    m_modalImpulseOut
            << ts << ' ' << body->id() << ' ' << cRec.vtxId << ' '
            << impVec.x << ' ' << impVec.y << ' ' << impVec.z << " S C" << endl;
  }

  // Get rest pose positions of the impact
  const Point3<REAL>     pointA = body->initial_predicted_position( cRec.pt );

  // Center of mass and rotation
  const Point3d         &posA = body->predicted_mass_center();

  const Quat4d          &qA = body->predicted_inverse_rotation();

  const Vector3d        &contactN = cRec.impulseDir;

  REAL                   impulseMagnitude = imp.length();
  REAL                   relativeSpeed = cRec.vnrel;

  // Write the impulse type (C for constraint), the time of the impulse and
  // an (object id, body vertex position) tuple for the object
  //
  // For the body, also write the current inverse rotation quaternion and
  // current center of mass position
  m_fout << "C " << ts << ' ';
  m_fout << body->id() << ' ';
  m_fout << pointA.x << ' ' << pointA.y << ' ' << pointA.z << ' ';
  m_fout << posA.x << ' ' << posA.y << ' ' << posA.z << ' ';
  m_fout << qA.w << ' ' << qA.v.x << ' ' << qA.v.y << ' ' << qA.v.z << ' ';

  // Write out impulse information: relative speed along the contact
  // direction, impulse magnitude, and impulse direction
  m_fout << relativeSpeed << ' ' << impulseMagnitude << ' ';
  m_fout << contactN.x << ' ' << contactN.y << ' ' << contactN.z << endl;
}

void RigidObjImpRecorder::record_inter_obj_impulse(
                                      REAL ts, const Vector3<REAL>& imp,
                                      const CollisionRec<REAL>& cRec,
                                      const TRigidBody* ba,
                                      const TRigidBody* bb)
{
  {
    Vector3<REAL> impVec = ba->predicted_inverse_rotation().rotate(imp);
    Vector3<REAL> relVel = ba->predicted_inverse_rotation().rotate(cRec.vrel);
    m_modalImpulseOut
            << ts << ' ' << ba->id() << ' ' << cRec.vtxId << ' '
            << impVec.x << ' ' << impVec.y << ' ' << impVec.z << " S P" << endl;
#if 0
    m_modalImpulseOut
            << ts << ' ' << ba->id() << ' ' << cRec.vtxId << ' '
            << impVec.x << ' ' << impVec.y << ' ' << impVec.z << ' '
            <<-relVel.x << ' ' <<-relVel.y << ' ' <<-relVel.z << " S" << endl;
#endif
    
    int vtxIdB = bb->kdtree().find_nearest(cRec.pt);
    impVec = bb->predicted_inverse_rotation().rotate(imp);
    relVel = bb->predicted_inverse_rotation().rotate(cRec.vrel);
    m_modalImpulseOut
            << ts << ' ' << bb->id() << ' ' << vtxIdB << ' '
            <<-impVec.x << ' ' <<-impVec.y << ' ' <<-impVec.z << " T P" << endl;
#if 0
    m_modalImpulseOut
            << ts << ' ' << bb->id() << ' ' << vtxIdB << ' '
            <<-impVec.x << ' ' <<-impVec.y << ' ' <<-impVec.z << ' '
            << relVel.x << ' ' << relVel.y << ' ' << relVel.z << " T" << endl;
#endif
  }

  // Get rest pose positions of the impact in both meshes
  const Point3<REAL>     pointA = ba->initial_predicted_position( cRec.pt );
  const Point3<REAL>     pointB = bb->initial_predicted_position( cRec.pt );

  // Center of mass and rotation for each object at collision time
  const Point3d         &posA = ba->predicted_mass_center();
  const Point3d         &posB = bb->predicted_mass_center();

  const Quat4d          &qA = ba->predicted_inverse_rotation();
  const Quat4d          &qB = bb->predicted_inverse_rotation();

  const Vector3d        &contactN = cRec.impulseDir;

  REAL                   impulseMagnitude = imp.length();
  REAL                   relativeSpeed = cRec.vnrel;

  // Write the type of impulse (P for object pair), the time of the
  // impulse, and (object id, body vertex position) tuples
  // for both objects
  //
  // For each body, also write the current inverse rotation quaternion, and
  // current center of mass position
  m_fout << "P " << ts << ' ';
  m_fout << ba->id() << ' ';
  m_fout << pointA.x << ' ' << pointA.y << ' ' << pointA.z << ' ';
  m_fout << posA.x << ' ' << posA.y << ' ' << posA.z << ' ';
  m_fout << qA.w << ' ' << qA.v.x << ' ' << qA.v.y << ' ' << qA.v.z << ' ';
  m_fout << bb->id() << ' ';
  m_fout << pointB.x << ' ' << pointB.y << ' ' << pointB.z << ' ';
  m_fout << posB.x << ' ' << posB.y << ' ' << posB.z << ' ';
  m_fout << qB.w << ' ' << qB.v.x << ' ' << qB.v.y << ' ' << qB.v.z << ' ';

  // Write out impulse information: relative speed along the contact
  // direction, impulse magnitude, and impulse direction
  m_fout << relativeSpeed << ' ' << impulseMagnitude << ' ';
  m_fout << contactN.x << ' ' << contactN.y << ' ' << contactN.z << endl;

#if 0
  {
    Vector3<REAL> impVec = ba->predicted_inverse_rotation().rotate(imp);
    const Point3d& posA = ba->predicted_mass_center();
    const Quat4d&    qA = ba->predicted_rotation();
    m_modalImpulseOut
            << ts << " P " << ba->id() << ' ' << cRec.vtxId << ' ' 
            << impVec.x    << ' ' << impVec.y    << ' ' << impVec.z    << ' '
            << cRec.vrel.x << ' ' << cRec.vrel.y << ' ' << cRec.vrel.z << ' '
               // relative velocity
            << posA.x      << ' ' << posA.y      << ' ' << posA.z      << ' ' 
            << qA.w        << ' ' << qA.v.x      << ' ' << qA.v.y      << ' '
            << qA.v.z << endl;

    impVec = bb->predicted_inverse_rotation().rotate(imp);
    const Point3<REAL> ipt = bb->initial_predicted_position(cRec.pt);
    const Point3d& posB = bb->predicted_mass_center();
    const Quat4d&    qB = bb->predicted_rotation();
    m_modalImpulseOut
            << ts << " P " << bb->id() << ' '
            << ipt.x       << ' ' << ipt.y       << ' ' << ipt.z       << ' '
            << impVec.x    << ' ' << impVec.y    << ' ' << impVec.z    << ' '
            << posB.x      << ' ' << posB.y      << ' ' << posB.z      << ' '
            << qB.w        << ' ' << qB.v.x      << ' ' << qB.v.y      << ' '
            << qB.v.z << endl;
  }
#endif
}

/*
 * \pt is the point where the impulse is applied in the current "predicted" 
 * configuration. This method looks up the nearest vertex to \pt in the
 * vertices of the given rigid body \body, and records the impulse at that 
 * vertex
 */
/*
void RigidObjImpRecorder::record_impulse(REAL ts, 
        const Vector3<REAL>& imp, const Point3<REAL>& pt,
        const TRigidBody* body)
{
#ifdef USE_RECORDER
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);
    //// transform \pt from "predicted" configuration into rest configuration
    const Point3<REAL> ipt = body->initial_predicted_position(pt);
    //// find the nearest vertex from \ipt
    //   vid is the id in tet mesh because the kd-tree is built from the entire 
    //   tet vertices
    int vtxId = body->kdtree().find_nearest(ipt);

    m_fout << ts << ' ' << body->id() << ' '
           << vtxId << ' '
           << impVec.x << ' '
           << impVec.y << ' '
           << impVec.z << " S" << std::endl;
#endif
}
*/

