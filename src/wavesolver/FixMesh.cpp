#include "FixMesh.h"
#ifdef MESHFIX_WITH_EIGEN
#include "meshfix.h"
#endif

bool myMeshFix(const Eigen::MatrixXd & V,
               const Eigen::MatrixXi & F,
               Eigen::MatrixXd & W,
               Eigen::MatrixXi & G,
               bool keepAllComponents)
{
#ifdef MESHFIX_WITH_EIGEN
    return meshfix(V, F, W, G, keepAllComponents);
#else
    return false;
#endif
}

