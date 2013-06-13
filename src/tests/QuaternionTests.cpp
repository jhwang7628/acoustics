#include <gtest/gtest.h>
#include "linearalgebra/Quaternion.hpp"

using namespace std;

TEST(TestQuaternion, rotate1)
{
    Point3d a(0, 0, 1);
    Quat4d rot = Quat4d::fromAxisRotD(Vector3d(0,0,1), 90) * Quat4d::fromAxisRotD(Vector3d(0,1,0),90);
    cout << "ROT RET-1: " << rot.rotate(a) << endl;
}

TEST(TestQuaternion, composite1)
{
    double alpha = M_PI/3, beta = M_PI/4, gamma = M_PI/5;

    Matrix3<double> Rx(1., 0., 0.,
                       0., cos(gamma), -sin(gamma),
                       0., sin(gamma),  cos(gamma));
    Matrix3<double> Ry(cos(beta), 0., sin(beta),
                       0.,        1., 0.,
                      -sin(beta), 0., cos(beta));
    Matrix3<double> Rz(cos(alpha), -sin(alpha), 0.,
                       sin(alpha),  cos(alpha), 0.,
                        0.,         0.,         1.);
    Vector3d vec(0.1, 0.3, 0.4);
    Vector3d out1 = (Rz*(Ry*(Rx*vec)));

    Quat4d rot = Quat4d::fromAxisRotR(Vector3d(0.,0.,1.), alpha) *
                 Quat4d::fromAxisRotR(Vector3d(0.,1.,0.), beta) *
                 Quat4d::fromAxisRotR(Vector3d(1.,0.,0.), gamma);
    Vector3d out2 = rot.rotate(vec);
    EXPECT_NEAR(out1.x, out2.x, 1E-10);
    EXPECT_NEAR(out1.y, out2.y, 1E-10);
    EXPECT_NEAR(out1.z, out2.z, 1E-10);

    Matrix3<double> mmat = rot.toMatrix3();
    cout << mmat << endl;
    cout << "=================================" << endl;
    cout << Rz*Ry*Rx << endl;
}
