set(DISTANCEFIELD_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/boundingBox.h         ${CMAKE_CURRENT_SOURCE_DIR}/boundingBox.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/closestPointField.h   ${CMAKE_CURRENT_SOURCE_DIR}/closestPointField.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/distanceField.h       ${CMAKE_CURRENT_SOURCE_DIR}/distanceField.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mat3d.h               ${CMAKE_CURRENT_SOURCE_DIR}/mat3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/minivector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/objfile.h             ${CMAKE_CURRENT_SOURCE_DIR}/objfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/objfileOctree.h       ${CMAKE_CURRENT_SOURCE_DIR}/objfileOctree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/obj.h
    ${CMAKE_CURRENT_SOURCE_DIR}/octree.h              ${CMAKE_CURRENT_SOURCE_DIR}/octree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sphere.h              ${CMAKE_CURRENT_SOURCE_DIR}/sphere.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/triangle.h            ${CMAKE_CURRENT_SOURCE_DIR}/triangle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tribox2.h             ${CMAKE_CURRENT_SOURCE_DIR}/tribox2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/triple.h
    ${CMAKE_CURRENT_SOURCE_DIR}/vec2d.h               ${CMAKE_CURRENT_SOURCE_DIR}/vec2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vec3d.h               ${CMAKE_CURRENT_SOURCE_DIR}/vec3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vec_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FieldBuilder.h        ${CMAKE_CURRENT_SOURCE_DIR}/FieldBuilder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AdaptiveDistanceField.h        
    ${CMAKE_CURRENT_SOURCE_DIR}/AdaptiveDistanceField.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lodepng.h        
    ${CMAKE_CURRENT_SOURCE_DIR}/lodepng.cpp

    CACHE FILE "Distance Field sources")
