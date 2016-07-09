# =============================================================================
#   Builds the unit testing tools
# =============================================================================
if (BUILD_UNIT_TESTING)
    # ===================
    # Dependencies
    # ===================
    set(GENERAL_TESTING_SRCS
        testing.h) 

    # ===================
    # Settings
    # ===================
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    
    # ===================
    # Definitions
    # ===================
    add_definitions(-DUSE_RECORDER) # needed for testing recorder

    # ===================
    # Target
    # ===================
    add_executable(unit_testing_leastSquareSurface
        ${GENERAL_TESTING_SRCS}
        ${LEASTSQUARESURFACE_SOURCES}
        test_LeastSquareSurface.cpp
        )
    add_executable(unit_testing_grid
        ${GENERAL_TESTING_SRCS}
        ${GRID_SOURCES}
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
        test_GridWithObject.cpp
        )

    add_executable(unit_testing_graph
        ${GENERAL_TESTING_SRCS}
        ${GRAPH_TESTING_SRCS}
        test_Graph.cpp
        ) 

    add_executable(unit_testing_sparseLinearSystemSolver
        ${GENERAL_TESTING_SRCS}
        ${SPARSE_LINEAR_ALGEBRA_SRC} 
        test_SparseLinearSystem.cpp
        ) 

    add_executable(unit_testing_trilinear
        test_Trilinear.cpp
        ) 

    add_executable(unit_testing_FDTD_AcousticSimulator
        ${DISTANCEFIELD_SRC}              
        ${WAVESOLVER_SRC} 
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
        test_FDTD_AcousticSimulator.cpp
        ) 

    add_executable(unit_testing_Modal
        test_Modal.cpp
        ${GENERAL_TESTING_SRCS} 
        ${LINEARALGEBRA_SRC} 
        ${DEFORMABLE_SRC} 
        ${RIGID_SRC}
        ${RIGID_IO_SRC} 
        ${GEOMETRY_SRC}
        ${LOGGING_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  
        ) 

    # ===================
    # Linking 
    # ===================
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

    target_link_libraries(unit_testing_FDTD_AcousticSimulator
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        boost_thread
        boost_system                    
        boost_timer
        boost_filesystem
        boost_regex) 

    target_link_libraries(unit_testing_Modal
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             
        ${GTS_LIBRARIES}
        boost_thread
        boost_system                    
        boost_timer
        boost_filesystem
        boost_regex) 


endif ( BUILD_UNIT_TESTING )
