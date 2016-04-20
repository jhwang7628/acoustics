# =============================================================================
#   Builds the unit testing tools
# =============================================================================
if (BUILD_UNIT_TESTING)
    # ===================
    # Dependencies
    # ===================

    set(GENERAL_TESTING_SRCS
        testing.h) 

    set(GRID_TESTING_SRCS 
        test_GridWithObject.cpp)

    set(LEASTSQUARESURFACE_TESTING_SRCS 
        test_LeastSquareSurface.cpp)

    set(GRAPH_TESTING_SRCS
        test_Graph.cpp) 

    set(SPARSE_LINEAR_SYSTEM_SOLVER_TESTING_SRCS
        test_SparseLinearSystem.cpp) 

    set(TRILINEAR_TESTING_SRCS
        test_Trilinear.cpp)

    set(FDTD_RIGID_OBJECT_TESTING_SRCS
        test_FDTD_RigidObject.cpp)

    set(TEST_INCLUDE
        /home/jui-hsien/code/acoustics/src/TYPES.h
        )


    # set(TURBSOUND_INCLUDE_DIR
    #     /home/jui-hsien/code/turbsound_postfluid/src
    #     /home/jui-hsien/code/turbsound_postfluid/tools)


    # ===================
    # Settings
    # ===================
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


    # ===================
    # Target
    # ===================


    add_executable(unit_testing_leastSquareSurface
        ${GENERAL_TESTING_SRCS}
        ${LEASTSQUARESURFACE_TESTING_SRCS}
        ${LEASTSQUARESURFACE_SOURCES}
        )
    add_executable(unit_testing_grid
        ${GENERAL_TESTING_SRCS}
        ${GRID_TESTING_SRCS}
        ${GRID_SOURCES}
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
        )

    add_executable(unit_testing_graph
        ${GENERAL_TESTING_SRCS}
        ${GRAPH_TESTING_SRCS}
        ${GRAPH_SRC} 
        ) 

    add_executable(unit_testing_sparseLinearSystemSolver
        ${GENERAL_TESTING_SRCS}
        ${SPARSE_LINEAR_SYSTEM_SOLVER_TESTING_SRCS}
        ${SPARSE_LINEAR_ALGEBRA_SRC} 
        ) 

    add_executable(unit_testing_trilinear
        ${TRILINEAR_TESTING_SRCS}
        ) 

    # message(${GENERAL_TESTING_SRCS}           "\n")
    # message(${FDTD_RIGID_OBJECT_TESTING_SRCS} "\n")
    # message(${GENERIC_SRC}                    "\n")
    # message(${IO_SRC}                         "\n")
    # message(${LOGGING_SRC}                    "\n")
    # message(${LINEARALGEBRA_SRC}              "\n")
    # message(${DISTANCEFIELD_SRC}              "\n")
    # message(${MATH_SRC}                       "\n")
    # message(${LINEARALGEBRA_BLAS_LAPACK_SRC}  "\n")
    # message(${GEOMETRY_SRC}                   "\n")
    # message(${FIELD_SRC}                      "\n")
    # message(${WAVESOLVER_SRC}                 "\n")
    # message(${UTILS_SRC}                      "\n")
    # message(${PARSER_SRC}                     "\n")
    # message(${GTS_GEOMETRY_SRC}               "\n")
    # message(${CGAL_GEOMETRY_SRC}              "\n")
    # message(${TRANSFER_SRC}                   "\n")
    # message(${WAVELET_TRANSFER_SRC}           "\n")
    # message(${SPARSE_LINEAR_ALGEBRA_SRC}      "\n")

    add_executable(unit_testing_FDTD_RigidObject
        ${DISTANCEFIELD_SRC}              
        ${WAVESOLVER_SRC} 
        ${GENERAL_TESTING_SRCS} 
        ${FDTD_RIGID_OBJECT_TESTING_SRCS} 
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
        ) 

    target_link_libraries(unit_testing_sparseLinearSystemSolver
        boost_timer
        boost_thread
        boost_system                    
        boost_regex) 

    target_link_libraries(unit_testing_grid
        ${MKL_LIBRARY}
        ${LEASTSQUARESURFACE_SOURCES}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        boost_thread
        boost_system                    
        boost_regex) 

    target_link_libraries(unit_testing_FDTD_RigidObject
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        boost_thread
        boost_system                    
        boost_timer
        boost_regex) 

    # include_directories(${isostuffer_SOURCE_DIR}/src) 
    # target_include_directories(
    #     unit_testing_FDTD_RigidObject2
    #     PUBLIC ${isostuffer_SOURCE_DIR}/src) 

endif ( BUILD_UNIT_TESTING )
