#include "FixMesh.h"
#include "meshfix.h"

bool myMeshFix(const Eigen::MatrixXd & V,
               const Eigen::MatrixXi & F,
               Eigen::MatrixXd & W,
               Eigen::MatrixXi & G,
               bool keepAllComponents)
{
    return meshfix(V, F, W, G, keepAllComponents);
}

