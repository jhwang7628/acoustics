set(GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidMesh.h                ${CMAKE_CURRENT_SOURCE_DIR}/RigidMesh.cpp
    CACHE FILE "Geometry sources")

set(GTS_GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/GTS_TriMesh.h              ${CMAKE_CURRENT_SOURCE_DIR}/GTS_TriMesh.cpp
    CACHE FILE "GTS Geometry sources")

set(CGAL_GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/MeshTree.h                 ${CMAKE_CURRENT_SOURCE_DIR}/MeshTree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ClosestPointMesh.h         ${CMAKE_CURRENT_SOURCE_DIR}/ClosestPointMesh.cpp
    CACHE FILE "CGAL Geometry sources")