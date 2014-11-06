# =============================================================================
#   Builds a GUI for the PAN Wavesolver
# =============================================================================
if (BUILD_PAN_PRECOMPUTE_GUI)
    # ===================
    # Dependencies
    # ===================
    
        # package for opengl and glut
        find_package(OpenGL REQUIRED)
        find_package(GLUT REQUIRED)

        if(USE_QT5)
            find_package(Qt5Widgets REQUIRED)
        else(USE_QT5)
            find_package(Qt4 REQUIRED)
        endif(USE_QT5)

        if (USE_CUDA)
		    find_package(CUDA)
		    add_definitions(-DUSE_CUDA)
		endif (USE_CUDA)

        # check the existence of libQGLViewer
        find_package(QGLViewer 1.0 REQUIRED)

        set(PAN_PRECOMPUTE_SRCS
            precompute_acceleration_pulse.cpp
            ../../ui/WaveViewer.h 			../../ui/WaveViewer.cpp
            )

    # ===================
    # Settings
    # ===================
        set(QT_USE_QTXML 1)
        set(QT_USE_QTOPENGL 1)

    # ===================
    # Target
    # ===================

        if(USE_QT5)
            include_directories(${QT_USE_FILE})
        else(USE_QT5)
            include(${QT_USE_FILE})
        endif(USE_QT5)

        set(CMAKE_INCLUDE_CURRENT_DIR ON)

        if(USE_CUDA)
	        set(CUDA_64_BIT_DEVICE_CODE ON)
	        set(CUDA_NVCC_FLAGS -arch=sm_12 --maxrregcount=32 -O4 -use_fast_math -Xptxas -dlcm=ca)
	        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \"-DBOOST_NOINLINE=__attribute__((noinline))\" ")

	        cuda_add_executable(precompute-pan-gui
	            ${PAN_PRECOMPUTE_SRCS}
	            ${GENERIC_SRC}                    ${IO_SRC}
	            ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
	            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
	            ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
	            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
	            ${WAVESOLVER_CUDA_PAN_SRC} 		  ${WAVESOLVER_CUDA_PAN_WRAPPER_SRC}
	            ${UTILS_SRC}                      ${PARSER_SRC}
	            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
	            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
	            ${DEFORMABLE_SRC}                 
	            ${MULTIPOLE_SRC})

	        target_link_libraries(precompute-pan-gui
	            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
	            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
	            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
	            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
	            ${GSL_LIBRARY}                  boost_thread
	            boost_system                    boost_serialization
	            ${MKL_LIBRARY}
	            ${CUDA_LIBRARIES}               )
	    else(USE_CUDA)
	    	add_executable(precompute-pan-gui
	            ${PAN_PRECOMPUTE_SRCS}
	            ${GENERIC_SRC}                    ${IO_SRC}
	            ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
	            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
	            ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
	            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
	            ${UTILS_SRC}                      ${PARSER_SRC}
	            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
	            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
	            )

	        target_link_libraries(precompute-pan-gui
	            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
	            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
	            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
	            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
	            ${GSL_LIBRARY}                  boost_thread
	            boost_system                    ${MKL_LIBRARY})
	    endif(USE_CUDA)

	    target_include_directories(precompute-pan-gui
            PUBLIC ${GLUT_INCLUDE_DIR}
            PUBLIC ${OPENGL_INCLUDE_DIR}
            PUBLIC ${QGLViewer_INCLUDE_DIR}
            )

        set_target_properties(precompute-pan-gui PROPERTIES AUTOMOC TRUE)
        set_target_properties(precompute-pan-gui PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
        set_target_properties(precompute-pan-gui PROPERTIES AUTOUIC TRUE)
        set_target_properties(precompute-pan-gui PROPERTIES AUTORCC TRUE)

        if(USE_QT5)
            qt5_use_modules(precompute-pan-gui Xml OpenGL)
        endif(USE_QT5)

        set_target_properties(precompute-pan-gui PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_PAN_PRECOMPUTE_GUI )