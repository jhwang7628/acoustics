if (BUILD_PYTHON_EXTENSIONS)
    find_package(PythonLibs REQUIRED)
    FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
        
        if(USE_CUDA)
            cuda_add_library(pat_solver SHARED
                    pat_wrapper.cpp

                    ${GENERIC_SRC}                    ${IO_SRC}
                    ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
                    ${DISTANCEFIELD_SRC}              ${MATH_SRC}
                    ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
                    ${FIELD_SRC}                      ${WAVESOLVER_SRC}
                    ${WAVESOLVER_CUDA_PAT_SRC}        ${WAVESOLVER_CUDA_PAT_WRAPPER_SRC}
                    ${UTILS_SRC}                      ${PARSER_SRC}
                    ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
                    ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
                    ${DEFORMABLE_SRC}                 
                    ${MULTIPOLE_SRC})

            set_target_properties(pat_solver PROPERTIES PREFIX "_" OUTPUT_NAME solver)
            set_target_properties(pat_solver PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}/solver)

            target_link_libraries(pat_solver
                    ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
                    ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
                    ${GSL_LIBRARY}                  boost_thread
                    boost_system                    boost_serialization
                    ${MKL_LIBRARY}
                    ${CUDA_LIBRARIES}               
                    ${PYTHON_LIBRARIES}             ${Boost_LIBRARIES})

            file(WRITE "${LIBRARY_OUTPUT_PATH}/solver/__init__.py" "from _solver import *\n")

            add_dependencies(acoustic_python mode_data)
        else(USE_CUDA)
            message(SEND_ERROR "The Python Wave Solver requires CUDA")
        endif(USE_CUDA)

endif(BUILD_PYTHON_EXTENSIONS)