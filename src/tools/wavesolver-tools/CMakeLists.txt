# =============================================================================
#   Builds tools for processing FDTD wavesolver data
# =============================================================================

if ( BUILD_KLEINPAT_TOOLS )
    add_executable(extract_surface_mesh
        ${GENERIC_SRC}
        ${MATH_SRC}
        ${IO_SRC}
        ${LINEARALGEBRA_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}
        ${GEOMETRY_SRC})
    message("FIXME debug BLAS_LAPACK_LIBS = ${BLAS_LAPACK_LIBS}")
    target_link_libraries(extract_surface_mesh PUBLIC
        ${BLAS_LAPACK_LIBS}
        #mkl::mkl
        igl::core)

    #add_executable(extract_surface_mesh
    #    ${DISTANCEFIELD_SRC}
    #    ${WAVESOLVER_SRC}
    #    ${MESH_UTILS_SRC}
    #    ${GENERAL_TESTING_SRCS}
    #    ${FDTD_ACOUSTIC_SIMULATOR_TESTING_SRCS}
    #    ${GENERIC_SRC}
    #    ${IO_SRC}
    #    ${LOGGING_SRC}
    #    ${LINEARALGEBRA_SRC}
    #    ${MATH_SRC}
    #    ${LINEARALGEBRA_BLAS_LAPACK_SRC}
    #    ${GEOMETRY_SRC}
    #    ${FIELD_SRC}
    #    ${UTILS_SRC}
    #    ${PARSER_SRC}
    #    ${GTS_GEOMETRY_SRC}
    #    ${CGAL_GEOMETRY_SRC}
    #    ${TRANSFER_SRC}
    #    ${WAVELET_TRANSFER_SRC}
    #    ${SPARSE_LINEAR_ALGEBRA_SRC}
    #    ${MODAL_MODEL_SRC}
    #    ${MODAL_MATERIAL_SRC}
    #    ${DEFORMABLE_SRC}
    #    ${RIGID_SRC}
    #    ${RIGID_IO_SRC}
    #    ${IMPULSE_SERIES_SRC}
    #    ${RIGID_PARSER_SRC}
    #    ${SND_IO_SRC}
    #    extract_surface_mesh.cpp)

    #target_link_libraries(extract_surface_mesh
    #    ${MKL_LIBRARY}
    #    ${BLAS_LAPACK_LIBS}
    #    ${GTS_LIBRARIES}
    #    ${GLIB_LIBRARIES}
    #    ${CGAL_LIBRARIES}
    #    ${GSL_LIBRARY}
    #    ${LibConfig_LIBRARY}
    #    ${VLFEAT_LIBRARY}
    #    sndfile
    #    boost_thread
    #    boost_system
    #    boost_timer
    #    boost_filesystem
    #    boost_regex)
endif ()


if ( BUILD_WAVESOLVER_TOOLS )

    # ===================
    # Definitions
    # ===================
    add_definitions(-DUSE_RECORDER) # needed for testing recorder

    # ===================
    # Target
    # ===================
    add_executable(convert-dat-vtk
        ${GRID_SOURCES}
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
        ${FIELD_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${RIGID_PARSER_SRC}               ${LIBCONFIG_PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${WAVESOLVER_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
        ${SPARSE_LINEAR_ALGEBRA_SRC}
        ${LEASTSQUARESURFACE_SOURCES}
        ${RIGID_SRC}
        ${RIGID_IO_SRC}
        ${IMPULSE_SERIES_SRC}
        ${MODAL_MODEL_SRC}
        ${DEFORMABLE_SRC}
        ${MESH_UTILS_SRC}
        ${SND_IO_SRC}
        ${DISTANCEFIELD_SRC}
        ${WAVESOLVER_SRC}
        ${MESH_UTILS_SRC}
        ${GENERAL_TESTING_SRCS}
        ${FDTD_ACOUSTIC_SIMULATOR_TESTING_SRCS}
        ${GENERIC_SRC}
        ${IO_SRC}
        ${LOGGING_SRC}
        ${LINEARALGEBRA_SRC}
        ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}
        ${GEOMETRY_SRC}
        ${FIELD_SRC}
        ${UTILS_SRC}
        ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}
        ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}
        ${WAVELET_TRANSFER_SRC}
        ${SPARSE_LINEAR_ALGEBRA_SRC}
        ${MODAL_MODEL_SRC}
        ${MODAL_MATERIAL_SRC}
        ${DEFORMABLE_SRC}
        ${RIGID_SRC}
        ${RIGID_IO_SRC}
        ${IMPULSE_SERIES_SRC}
        ${RIGID_PARSER_SRC}
        ${SND_IO_SRC}
        convert_dat_vtk.cpp)

    add_executable(convert-dat-vtk-nonuniform
        convert_dat_vtk_nonuniform.cpp)

    add_executable(kirchhoff_helmholtz
        kirchhoff_helmholtz.cpp)

    add_executable(find_speaker_vertices
        find_speaker_vertices.cpp)

    add_executable(triangulate_obj
        triangulate_obj.cpp)

    add_executable(remove_duplicates
        remove_duplicates.cpp)

    add_executable(remove_unused_vertices_obj
        remove_unused_vertices_obj.cpp)

    add_executable(obj_centroid
        obj_centroid.cpp)

    add_executable(read_write_test
        read_write_test.cpp)

    add_executable(upsample_poser
        upsample_poser.cpp)

    add_executable(extract_surface_mesh
        ${DISTANCEFIELD_SRC}
        ${WAVESOLVER_SRC}
        ${MESH_UTILS_SRC}
        ${GENERAL_TESTING_SRCS}
        ${FDTD_ACOUSTIC_SIMULATOR_TESTING_SRCS}
        ${GENERIC_SRC}
        ${IO_SRC}
        ${LOGGING_SRC}
        ${LINEARALGEBRA_SRC}
        ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}
        ${GEOMETRY_SRC}
        ${FIELD_SRC}
        ${UTILS_SRC}
        ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}
        ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}
        ${WAVELET_TRANSFER_SRC}
        ${SPARSE_LINEAR_ALGEBRA_SRC}
        ${MODAL_MODEL_SRC}
        ${MODAL_MATERIAL_SRC}
        ${DEFORMABLE_SRC}
        ${RIGID_SRC}
        ${RIGID_IO_SRC}
        ${IMPULSE_SERIES_SRC}
        ${RIGID_PARSER_SRC}
        ${SND_IO_SRC}
        extract_surface_mesh.cpp)

    add_executable(map_list_of_points_to_triangle_mesh
        map_list_of_points_to_triangle_mesh.cpp)

    add_executable(extract_surf_modes
        ${CMAKE_CURRENT_SOURCE_DIR}/../../geometry/TetMeshIndexToSurfaceMesh.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../../geometry/TetMeshIndexToSurfaceMesh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../../deformable/ModeData.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../../deformable/ModeData.cpp
        extract_surf_modes.cpp)

    target_link_libraries(convert-dat-vtk
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}
        ${LibConfig_LIBRARY}
        ${VLFEAT_LIBRARY}
        sndfile
        boost_thread
        boost_timer
        boost_system
        boost_filesystem
        boost_regex)

    target_link_libraries(convert-dat-vtk-nonuniform
        boost_thread
        boost_timer
        boost_system
        boost_filesystem
        boost_regex)

    target_link_libraries(upsample_poser
        boost_thread
        boost_timer
        boost_system
        boost_filesystem
        boost_regex)

    target_link_libraries(extract_surface_mesh
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}
        ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}
        ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}
        ${LibConfig_LIBRARY}
        ${VLFEAT_LIBRARY}
        sndfile
        boost_thread
        boost_system
        boost_timer
        boost_filesystem
        boost_regex)

endif ( BUILD_WAVESOLVER_TOOLS )
