# =============================================================================
#   Builds a GUI for the Impulse Response Wavesolver
# =============================================================================
if (BUILD_PRECOMPUTE_IMPULSE_RESPONSE_GUI)
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

        set(IR_PRECOMPUTE_GUI_SRCS
            precompute_impulse_response.cpp
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

        include_directories("/home/jui-hsien/code/turbsound_postfluid/tools")

        set(CMAKE_INCLUDE_CURRENT_DIR ON)

        if(USE_CUDA)
	        set(CUDA_64_BIT_DEVICE_CODE ON)
	        set(CUDA_NVCC_FLAGS -arch=sm_12 --maxrregcount=32 -O4 -use_fast_math -Xptxas -dlcm=ca)
	        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \"-DBOOST_NOINLINE=__attribute__((noinline))\" ")

	        cuda_add_executable(precompute-impulse-response
	            ${IR_PRECOMPUTE_GUI_SRCS}
	            ${GENERIC_SRC}                    ${IO_SRC}
	            ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
	            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
	            ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
	            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
	            ${WAVESOLVER_CUDA_PAT_SRC} 		  ${WAVESOLVER_CUDA_PAT_WRAPPER_SRC}
	            ${UTILS_SRC}                      ${PARSER_SRC}
	            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
	            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
                ${DEFORMABLE_SRC}                 ${EIGENIO_SRC}
	            ${MULTIPOLE_SRC})

	        target_link_libraries(precompute-impulse-response
	            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
	            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
	            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
	            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
	            ${GSL_LIBRARY}                  boost_thread
	            boost_system                    boost_serialization
	            ${MKL_LIBRARY}
	            ${CUDA_LIBRARIES}               )
	    else(USE_CUDA)
	    	add_executable(precompute-impulse-response
	            ${IR_PRECOMPUTE_GUI_SRCS}
	            ${GENERIC_SRC}                    ${IO_SRC}
	            ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
	            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
	            ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
	            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
	            ${UTILS_SRC}                      ${PARSER_SRC}
	            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
	            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
	            )

	        target_link_libraries(precompute-impulse-response
	            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
	            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
	            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
	            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
	            ${GSL_LIBRARY}                  boost_thread
	            boost_system                    ${MKL_LIBRARY})
	    endif(USE_CUDA)

	    target_include_directories(precompute-impulse-response
            PUBLIC ${GLUT_INCLUDE_DIR}
            PUBLIC ${OPENGL_INCLUDE_DIR}
            PUBLIC ${QGLViewer_INCLUDE_DIR}
            )

        set_target_properties(precompute-impulse-response PROPERTIES AUTOMOC TRUE)
        set_target_properties(precompute-impulse-response PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
        set_target_properties(precompute-impulse-response PROPERTIES AUTOUIC TRUE)
        set_target_properties(precompute-impulse-response PROPERTIES AUTORCC TRUE)

        if(USE_QT5)
            qt5_use_modules(precompute-impulse-response Xml OpenGL)
        endif(USE_QT5)

        set_target_properties(precompute-impulse-response PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_PRECOMPUTE_IMPULSE_RESPONSE_GUI )

# =============================================================================
#   Builds the Impulse Response Wavesolver
# =============================================================================
if (BUILD_PRECOMPUTE_IMPULSE_RESPONSE)

    set(CMAKE_CXX_FLAGS "-std=c++11")
    # ===================
    # Dependencies
    # ===================

        if (USE_CUDA)
		    find_package(CUDA)
		    add_definitions(-DUSE_CUDA)
		endif (USE_CUDA)

        set(IR_PRECOMPUTE_SRCS
            precompute_impulse_response_text.cpp
            )

    # ===================
    # Settings
    # ===================


    # ===================
    # Target
    # ===================

        if(USE_CUDA)
	        set(CUDA_64_BIT_DEVICE_CODE ON)
	        set(CUDA_NVCC_FLAGS -arch=sm_12 --maxrregcount=32 -O4 -use_fast_math -Xptxas -dlcm=ca)
	        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \"-DBOOST_NOINLINE=__attribute__((noinline))\" ")

	        cuda_add_executable(precompute-impulse-response-text
	            ${IR_PRECOMPUTE_SRCS}
	            ${GENERIC_SRC}                    ${IO_SRC}
	            ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
	            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
	            ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
	            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
	            ${WAVESOLVER_CUDA_PAT_SRC} 		  ${WAVESOLVER_CUDA_PAT_WRAPPER_SRC}
	            ${UTILS_SRC}                      ${PARSER_SRC}
	            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
	            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
	            ${DEFORMABLE_SRC}                 
	            ${MULTIPOLE_SRC})

	        target_link_libraries(precompute-impulse-response-text
	            ${MKL_LIBRARY}
	            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
	            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
	            ${GSL_LIBRARY}                  boost_thread
	            boost_system                    boost_serialization
	            ${CUDA_LIBRARIES})
	    else(USE_CUDA)
	    	add_executable(precompute-impulse-response-text
	            ${IR_PRECOMPUTE_SRCS}
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

	        target_link_libraries(precompute-impulse-response-text
                ${MKL_LIBRARY}
                ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
	            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
	            ${GSL_LIBRARY}                  boost_thread
	            boost_system                    boost_regex)
	    endif(USE_CUDA)

        set_target_properties(precompute-impulse-response-text PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_PRECOMPUTE_IMPULSE_RESPONSE )
