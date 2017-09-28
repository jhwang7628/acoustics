set(GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidMesh.h                ${CMAKE_CURRENT_SOURCE_DIR}/RigidMesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tribox3.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tribox3.cpp
    CACHE FILE "Geometry sources")

set(GTS_GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/GTS_TriMesh.h              ${CMAKE_CURRENT_SOURCE_DIR}/GTS_TriMesh.cpp
    CACHE FILE "GTS Geometry sources")

set(CGAL_GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/MeshTree.h                 ${CMAKE_CURRENT_SOURCE_DIR}/MeshTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ClosestPointMesh.h         ${CMAKE_CURRENT_SOURCE_DIR}/ClosestPointMesh.cpp
    CACHE FILE "CGAL Geometry sources")

set(MESH_UTILS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/TetMeshIndexToSurfaceMesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/TetMeshIndexToSurfaceMesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TriangleMeshGraph.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TriangleMeshGraph.cpp
    CACHE FILE "Mesh Utils sources")

if (USE_VLFEAT)
    set(MESH_UTILS_SRC
        ${MESH_UTILS_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/TriangleMeshKDTree.hpp
        CACHE FILE "Mesh Utils sources")
endif (USE_VLFEAT)

set(EIGENIO_SRC
    /home/jui-hsien/Research/code/TurbulentHead/tools/IO/IO.h
    CACHE FILE "Eigen IO sources")
