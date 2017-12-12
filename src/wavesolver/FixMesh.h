#ifndef _FIX_MESH_H
#define _FIX_MESH_H

#include <Eigen/Dense>

bool myMeshFix(const Eigen::MatrixXd & V,
               const Eigen::MatrixXi & F,
               Eigen::MatrixXd & W,
               Eigen::MatrixXi & G,
               bool keepAllComponents=false);

#endif // _FIX_MESH_H

