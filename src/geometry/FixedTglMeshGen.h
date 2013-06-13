/*
 * =====================================================================================
 *
 *       Filename:  FixedTglMeshGen.h
 *
 *    Description:  Generated a triangle mesh with given vertices are fixed
 *                  All the fixed vertices in the new mesh are put at the beginning
 *                  of the vertice list
 *
 *        Version:  1.0
 *        Created:  10/27/10 10:15:14
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  FIXED_TGL_MESH_GEN_INC
#define  FIXED_TGL_MESH_GEN_INC

#include "TriangleMesh.hpp"

void gen_fixed_tgl_mesh(const TriangleMesh<double>* inMsh, const bool* isFixed,
                        TriangleMesh<double>* outMsh, int* vidMap);

#endif   /* ----- #ifndef FIXED_TGL_MESH_GEN_INC  ----- */

