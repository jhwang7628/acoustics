set(DEFORMABLE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/linear.cpp         ${CMAKE_CURRENT_SOURCE_DIR}/ModeData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/stvk.cpp           ${CMAKE_CURRENT_SOURCE_DIR}/StVKMesh.cpp

    CACHE FILE "Deformable source")

add_subdirectory(python)