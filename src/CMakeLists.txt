## CMake specification for src directory

set(CONFIG_SOURCES config.h)

set(ISOSTUFFER_SOURCES
    isostuffer/tables.h            isostuffer/IsoStuffer.hpp
    isostuffer/OctTree.h           isostuffer/OctTree.cpp
    isostuffer/DistValProc.h       isostuffer/DistValProc.cpp)

# utils/
set(UTILS_SRC
    utils/macros.h      utils/math.hpp
    utils/nano_timer.h  utils/nano_timer.c
    utils/timer.hpp     utils/timer.cpp
    utils/MathUtil.cpp  utils/MathUtil.cpp)

# linearalgebra
set(LINEAR_SRC
    linearalgebra/eig3.h          linearalgebra/eig3.cpp
    linearalgebra/Quaternion.hpp  linearalgebra/Quaternion.cpp
    linearalgebra/mat4inv.cpp)

# Linear algebra sources that depend on BLAS/LAPACK implementations
set(LINEAR_BLAS_LAPACK_SRC
    linearalgebra/MATRIX.h        linearalgebra/MATRIX_FAST.cpp
    linearalgebra/VECTOR.h        linearalgebra/VECTOR_FAST.cpp
    linearalgebra/SPARSE_MATRIX.h linearalgebra/SPARSE_MATRIX.cpp)

set(GENERIC_SRC
    generic/precision_type.cpp)

set(LOGGING_SRC
    logging/logging.cpp)

set(DEFORMABLE_SRC
    deformable/linear.cpp         deformable/ModeData.cpp
    deformable/stvk.cpp           deformable/StVKMesh.cpp)

set(IO_SRC
    io/MatrixIO.cpp               io/TglMeshReader.cpp
    io/ImpulseIO.cpp)

set(RIGID_SRC
    rigid/Simulator.cpp           rigid/Simulator_QS.cpp
    rigid/ContactGraph.cpp)

set(RIGID_IO_SRC
    io/RigidObjImpRecorder.cpp    io/RigidObjRecorder.cpp)

set(FIELD_SRC
    field/ScalarField.h           field/ScalarField.cpp)

set(WAVESOLVER_SRC
    wavesolver/Laplacian.cpp      wavesolver/Laplacian.h
    wavesolver/MAC_Grid.cpp       wavesolver/MAC_Grid.h
    wavesolver/PML_WaveSolver.cpp wavesolver/PML_WaveSolver.h
    wavesolver/WaveSolver.cpp     wavesolver/WaveSolver.h)

set(WAVESOLVER_CUDA_WRAPPER_SRC
    wavesolver/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.h
    wavesolver/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.cpp
    wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h
    wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.cpp)

set(WAVESOLVER_CUDA_SRC
    wavesolver/gpusolver/cuda/cuda_PAN_wave_3d.h
    wavesolver/gpusolver/cuda/cuda_PAN_wave_3d.cu
    wavesolver/gpusolver/cuda/cuda_PAN_wave_3d_kernel.cu
    wavesolver/gpusolver/cuda/cuda_PAN_wave_3d_kernel_math.cu
    wavesolver/gpusolver/cuda/cuda_PAT_wave_3d.h
    wavesolver/gpusolver/cuda/cuda_PAT_wave_3d.cu
    wavesolver/gpusolver/cuda/cuda_PAT_wave_3d_kernel.cu
    wavesolver/gpusolver/cuda/cuda_PAT_wave_3d_kernel_math.cu)

set(SUPERLU_INTERFACE_SRC
    superlu-interface/SuperLU_Interface.h
    superlu-interface/SuperLU_Interface.cpp)

set(EIGENSOLVER_SRC
    eigensolver/areigs.cpp        eigensolver/areigs.h
    eigensolver/arinfo.cc
    eigensolver/Eigensolver.cpp   eigensolver/Eigensolver.h)

# Sources for Jernej Barbic's flat signed distance function classes
set(DISTANCEFIELD_SRC
    distancefield/boundingBox.h         distancefield/boundingBox.cpp
    distancefield/closestPointField.h   distancefield/closestPointField.cpp
    distancefield/distanceField.h       distancefield/distanceField.cpp
    distancefield/mat3d.h               distancefield/mat3d.cpp
    distancefield/minivector.h
    distancefield/objfile.h             distancefield/objfile.cpp
    distancefield/objfileOctree.h       distancefield/objfileOctree.cpp
    distancefield/obj.h
    distancefield/octree.h              distancefield/octree.cpp
    distancefield/sphere.h              distancefield/sphere.cpp
    distancefield/triangle.h            distancefield/triangle.cpp
    distancefield/tribox2.h             distancefield/tribox2.cpp
    distancefield/triple.h
    distancefield/vec2d.h               distancefield/vec2d.cpp
    distancefield/vec3d.h               distancefield/vec3d.cpp
    distancefield/vec_types.h
    distancefield/FieldBuilder.h        distancefield/FieldBuilder.cpp)

set(MATH_SRC
    math/InterpolationFunction.h        math/InterpolationFunction.cpp
    math/Function.h                     math/Function.cpp
    math/SampledFunction.h              math/SampledFunction.cpp)

set(GEOMETRY_SRC
    geometry/RigidMesh.h                geometry/RigidMesh.cpp)

set(GTS_GEOMETRY_SRC
    geometry/GTS_TriMesh.h              geometry/GTS_TriMesh.cpp)

set(CGAL_GEOMETRY_SRC
    geometry/MeshTree.h                 geometry/MeshTree.cpp
    geometry/ClosestPointMesh.h         geometry/ClosestPointMesh.cpp)

set(TRANSFER_SRC
    transfer/AccelerationNoiseModel.h   transfer/AccelerationNoiseModel.cpp
    transfer/MultiTermApproximation.h   transfer/MultiTermApproximation.cpp
    transfer/PulseApproximation.h       transfer/PulseApproximation.cpp
    transfer/ProxyManager.h             transfer/ProxyManager.cpp
    transfer/RadialApproximation.h      transfer/RadialApproximation.cpp)

set(WAVELET_TRANSFER_SRC
    transfer/CompressedMultiTermApproximation.h
    transfer/CompressedMultiTermApproximation.cpp)

# tinyXML-based Parser
set(PARSER_SRC
    parser/Parser.cpp                   parser/Parser.h
    tinyxml/tinystr.cpp                 tinyxml/tinystr.h
    tinyxml/tinyxml.cpp                 tinyxml/tinyxmlerror.cpp
    tinyxml/tinyxml.h                   tinyxml/tinyxmlparser.cpp)

# Sound generation sources
set(SNDGEN_SRC
    sndgen/FBemTransferInterp.h         sndgen/FBemTransferInterp.cpp
    sndgen/ImpulseInterp.hpp
    sndgen/RigidModal.h                 sndgen/RigidModal.cpp
    sndgen/RigidSoundObj.hpp
    sndgen/SoundBuffer.hpp
    sndgen/SoundWriter.hpp)

set(SNDPLAYER_SRC
    sndgen/player/SoundPlayer.h         sndgen/player/SoundPlayer.cpp
    sndgen/player/SineWavePlayer.h
    )

set(MULTIPOLE_SRC
    multipole/MultipoleData.h           multipole/MultipoleData.cpp
    multipole/MultipoleUtil.h           multipole/MultipoleUtil.cpp
    multipole/MultipoleMath.h           multipole/MultipoleMath.cpp
    multipole/MultipolePlayer.h         multipole/MultipolePlayer.cpp)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Boost_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

if ( USE_GTS )
    include_directories(
        ${GTS_INCLUDE_DIR}
        ${GLIB_INCLUDE_DIRS})
endif ( USE_GTS )

if ( USE_CGAL )
    include_directories(
        ${CGAL_INCLUDE_DIRS})
endif ( USE_CGAL )

if ( USE_GSL )
    include_directories(
        ${GSL_INCLUDE_DIR})
endif ( USE_GSL )

# =============================================================================
#   Isostuffer GUI
# =============================================================================
if ( BUILD_ISOSTUFFER_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    #include(${QT_USE_FILE})
    include_directories(${QT_USE_FILE})

    include_directories(
        ${GLUT_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
        ${QGLViewer_INCLUDE_DIR})

    # setup the qrc files
    set(UI_QRC_FILES ui/images.qrc)
    #QT4_ADD_RESOURCES(UI_QRC_SRCS ${UI_QRC_FILES})
    QT5_ADD_RESOURCES(UI_QRC_SRCS ${UI_QRC_FILES})

    # headers to be moc-ified
    set(UI_MOC_HEADERS 
        ui/IsoStufferFrame.h    ui/IsoStufferCanvas.h
        ui/ParamsDialog.h)
    # the ui files
    set(UI_UI_FILES ui/IsoStufferFrame.ui ui/Params.ui)

    # Modification made here to fix compiler error with Boost 1.48.
    # Original version of this command is:
    #       QT4_WRAP_CPP(UI_MOC_SRCS ${UI_MOC_HEADERS})
    # -Jeff Chadwick
    #QT4_WRAP_CPP(UI_MOC_SRCS ${UI_MOC_HEADERS} 
                 #OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
    #QT4_WRAP_UI(UI_UI_SRCS ${UI_UI_FILES})
    QT5_WRAP_CPP(UI_MOC_SRCS ${UI_MOC_HEADERS} 
                 OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
    QT5_WRAP_UI(UI_UI_SRCS ${UI_UI_FILES})

    set(GUIAPP_SRCS tools/isostuffer/GuiMain.cpp
        ui/IsoStufferFrame.h    ui/IsoStufferFrame.cpp
        io/TglMeshReader.hpp    io/TglMeshReader.cpp
        ui/IsoStufferCanvas.h   ui/IsoStufferCanvas.cpp
        ${UI_UI_SRCS} ${UI_MOC_SRCS} ${UI_QRC_SRCS} )
    
    # =========================================================================
    # set for the final executable for isostuffer
    add_executable(isostuffer-gui
        ${CONFIG_SOURCES}       ${GUIAPP_SRCS}          ${ISOSTUFFER_SOURCES}
        ${UTILS_SRC}            ${LINEAR_SRC}           ${GENERIC_SRC}
        ${IO_SRC})

    target_link_libraries(isostuffer-gui
        ${GLUT_LIBRARIES}       ${QT_LIBRARIES}         ${QGLViewer_LIBRARY}
        ${ICC_LIBS}             ${OPENGL_LIBRARIES}     rt
        ${MKL_LIBRARY})

    qt5_use_modules(isostuffer-gui Xml OpenGL)

    set_target_properties(isostuffer-gui PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_ISOSTUFFER_GUI ) 


# =============================================================================
#   Command-line isostuffer
# =============================================================================
if ( BUILD_ISOSTUFFER )
  add_executable(isostuffer
      io/TglMeshReader.hpp    io/TglMeshReader.cpp
      tools/isostuffer/main.cpp
      ${CONFIG_SOURCES}           ${ISOSTUFFER_SOURCES}       ${UTILS_SRC}
      ${LINEAR_SRC}               ${GENERIC_SRC}              ${IO_SRC})

  target_link_libraries(isostuffer ${ICC_LIBS}        boost_program_options
                        rt)

  set_target_properties(isostuffer PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_ISOSTUFFER )

# =============================================================================
#   Tetmesh viewer
# =============================================================================
if ( BUILD_TETVIEWER )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    #include(${QT_USE_FILE})
    include_directories(${QT_USE_FILE})
    include_directories(
        ${GLUT_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
        ${QGLViewer_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR})
    
    # =================== TETVIEWER ==================
    # the qrc files
    set(TETVIEWER_QRC_FILES tools/tetviewer/images.qrc)
    #QT4_ADD_RESOURCES(TETVIEWER_QRC_SRCS ${TETVIEWER_QRC_FILES})

    #QT4_WRAP_UI(TETVIEWER_UI_SRC tools/tetviewer/tetviewer.ui)
    #QT4_WRAP_CPP(TETVIEWER_MOC_SRC 
            #tools/tetviewer/TetViewerFrame.h tools/tetviewer/TetViewerCanvas.h)
    QT5_ADD_RESOURCES(TETVIEWER_QRC_SRCS ${TETVIEWER_QRC_FILES})

    QT5_WRAP_UI(TETVIEWER_UI_SRC tools/tetviewer/tetviewer.ui)
    QT5_WRAP_CPP(TETVIEWER_MOC_SRC 
            tools/tetviewer/TetViewerFrame.h tools/tetviewer/TetViewerCanvas.h)

    add_executable(tetviewer
        tools/tetviewer/TetViewerFrame.h      tools/tetviewer/TetViewerFrame.cpp
        tools/tetviewer/TetViewerCanvas.h     tools/tetviewer/TetViewerCanvas.cpp
        generic/precision_type.hpp  generic/precision_type.cpp
        ${TETVIEWER_QRC_SRCS}
        ${TETVIEWER_UI_SRC}         ${TETVIEWER_MOC_SRC}
        ${IO_SRC}                   ${DEFORMABLE_SRC}
        ${LINEAR_SRC})
    target_link_libraries(tetviewer
        #${GLUT_LIBRARIES}    ${OPENGL_LIBRARIES}
        #${QT_LIBRARIES}      ${QGLViewer_LIBRARY})
        ${GLUT_LIBRARIES}       ${QT_LIBRARIES}         ${QGLViewer_LIBRARY}
        ${ICC_LIBS}             ${OPENGL_LIBRARIES}     #${BLAS_LAPACK_LIBS}
        rt)

    qt5_use_modules(tetviewer Xml OpenGL)
endif ( BUILD_TETVIEWER )

# =============================================================================
#   Elasticity solver (builds mass and stiffness matrices)
# =============================================================================
if ( BUILD_ELASTICITY_SOLVER )
    link_directories(${MKL_LIBRARY_DIR})

    add_executable(elasticity_solver
        tools/elasticity/elasticity_solver.cpp
        ${LOGGING_SRC}        ${UTILS_SRC}
        ${DEFORMABLE_SRC}     ${LINEAR_SRC}
        ${GENERIC_SRC}        ${IO_SRC})
    target_link_libraries(elasticity_solver
        ${BLAS_LAPACK_LIBS}   boost_program_options
        rt                    ${MKL_LIBRARY})
endif ( BUILD_ELASTICITY_SOLVER )

# =============================================================================
#   Modal eigensolver
# =============================================================================
if ( BUILD_EIGENSOLVER )
    include_directories(
        #${ARPACK_INCLUDE_DIR})
        #/home/chadwick/packages/arpack++/include)
        arpack++/include)

    include_directories(
        ${SUPERLU_INCLUDES})
        
    message(STATUS "YOO: " ${SUPERLU_INCLUDES})

    link_directories(${MKL_LIBRARY_DIR})

    add_executable(eigensolver
        tools/eigensolver/eigensolver.cpp
        ${LOGGING_SRC}        ${UTILS_SRC}
        ${IO_SRC}             ${SUPERLU_INTERFACE_SRC}
        ${LINEAR_BLAS_LAPACK_SRC})
    target_link_libraries(eigensolver
        ${SUPERLU_LIBRARIES}  ${ARPACK_LIBRARY}
        ${BLAS_LAPACK_LIBS}   ${MKL_LIBRARY}
        boost_program_options rt)

    add_executable(arpack-eigensolver
        tools/eigensolver/arpack-eigensolver.cpp
        ${LOGGING_SRC}        ${UTILS_SRC}
        ${IO_SRC}             ${SUPERLU_INTERFACE_SRC}
        ${EIGENSOLVER_SRC}
        ${LINEAR_BLAS_LAPACK_SRC})
    target_link_libraries(arpack-eigensolver
        ${SUPERLU_LIBRARIES}  ${ARPACK_LIBRARY}
        ${BLAS_LAPACK_LIBS}   ${MKL_LIBRARY}
        boost_program_options rt)

    add_executable(test-superlu
        tools/test-superlu.cpp
        ${LOGGING_SRC}        ${UTILS_SRC}
        ${IO_SRC}             ${SUPERLU_INTERFACE_SRC}
        ${EIGENSOLVER_SRC}    ${LINEAR_BLAS_LAPACK_SRC})
    target_link_libraries(test-superlu
        ${SUPERLU_LIBRARIES}  ${ARPACK_LIBRARY}
        ${BLAS_LAPACK_LIBS}   ${MKL_LIBRARY}
        boost_program_options rt)
endif ( BUILD_EIGENSOLVER )

# =============================================================================
#   Simple utility for generating a tet mesh with one tet
# =============================================================================
if ( BUILD_SINGLE_TET )
    add_executable(single_tet
        tools/single_tet.cpp
        ${UTILS_SRC}          ${IO_SRC}
        ${LINEAR_SRC}         ${GENERIC_SRC})

    target_link_libraries(single_tet
        rt)
endif ( BUILD_SINGLE_TET )

# =============================================================================
#   Builds input files for the FastBEM acoustics solver
# =============================================================================
if ( BUILD_FBEM_INPUT_GEN )
    add_executable(fbem_input_gen
        tools/fbem_input_gen/fbem_input_gen.cpp
        ${UTILS_SRC}        ${IO_SRC}
        ${GENERIC_SRC})

    target_link_libraries(fbem_input_gen
        boost_program_options   boost_filesystem  boost_system  rt)

    add_executable(fieldpoints_gen
        tools/fbem_input_gen/fieldpoints_gen.cpp
        ${UTILS_SRC}        ${IO_SRC}
        ${GENERIC_SRC}      ${LOGGING_SRC})

    target_link_libraries(fieldpoints_gen
        boost_program_options   boost_filesystem  boost_system  rt)
endif ( BUILD_FBEM_INPUT_GEN )

# =============================================================================
#   Builds a rigid-body solver
# =============================================================================
if ( BUILD_RIGIDSIM )
    add_executable(rigidsim
        tools/RigidSimDrop.cpp
        ${GENERIC_SRC}                    ${DEFORMABLE_SRC}
        ${IO_SRC}                         ${RIGID_SRC}
        ${LOGGING_SRC}                    ${LINEAR_SRC}
        ${RIGID_IO_SRC})

    # Turn on impact and displacement recording
    add_definitions(-DUSE_RECORDER)

    target_link_libraries(rigidsim
        ${LibConfig_LIBRARY}              ${BLAS_LAPACK_LIBS}
        ${MKL_LIBRARY})
endif ( BUILD_RIGIDSIM )

# =============================================================================
#   Builds a rigid-body solver GUI
# =============================================================================
if ( BUILD_RIGIDSIM_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    include_directories(${QT_USE_FILE})

    # Set up the QRC files
    set(RIGIDSIM_UI_QRC_FILES tools/rigidsim-gui/images.qrc)
    #QT4_ADD_RESOURCES(RIGIDSIM_UI_QRC_SRCS ${RIGIDSIM_UI_QRC_FILES})
    QT5_ADD_RESOURCES(RIGIDSIM_UI_QRC_SRCS ${RIGIDSIM_UI_QRC_FILES})

    # Headers to be moc-ified
    set(RIGIDSIM_UI_MOC_HEADERS tools/rigidsim-gui/RigidSim.h)
    # The ui files
    set(RIGIDSIM_UI_UI_FILES tools/rigidsim-gui/RigidSim.ui)

    add_definitions(-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)

    #QT4_WRAP_CPP(RIGIDSIM_UI_MOC_SRCS ${RIGIDSIM_UI_MOC_HEADERS})
    #QT4_WRAP_UI(RIGIDSIM_UI_UI_SRCS ${RIGIDSIM_UI_UI_FILES})
    QT5_WRAP_CPP(RIGIDSIM_UI_MOC_SRCS ${RIGIDSIM_UI_MOC_HEADERS})
    QT5_WRAP_UI(RIGIDSIM_UI_UI_SRCS ${RIGIDSIM_UI_UI_FILES})

    # Turn on impact and displacement recording
    add_definitions(-DUSE_RECORDER)

    add_executable(rigidsim-gui
        tools/rigidsim-gui/RigidSim.cpp   tools/rigidsim-gui/RigidCanvas.cpp
        tools/rigidsim-gui/demo/DemoDropObjs.cpp
        tools/rigidsim-gui/demo/DemoDropObjsWithFixed.cpp
        tools/rigidsim-gui/demo/DemoShakingPiggy.cpp
        ${GENERIC_SRC}                    ${DEFORMABLE_SRC}
        ${IO_SRC}                         ${RIGID_SRC}
        ${LOGGING_SRC}                    ${LINEAR_SRC}
        ${RIGID_IO_SRC}                   ${RIGIDSIM_UI_UI_SRCS}
        ${RIGIDSIM_UI_MOC_SRCS}           ${RIGIDSIM_UI_QRC_SRCS})

    target_link_libraries(rigidsim-gui
        ${GLUT_LIBRARIES}                 ${QT_LIBRARIES}
        ${QGLViewer_LIBRARY}              ${OPENGL_LIBRARIES}
        ${LibConfig_LIBRARY}              ${BLAS_LAPACK_LIBS}
        ${MKL_LIBRARY})

    qt5_use_modules(rigidsim-gui Xml OpenGL)

    # ========================================================================
    #   Also build a simple "replayer" that we can use to make quick OpenGL
    #   renderings
    # ========================================================================
    #QT4_WRAP_CPP(REPLAYER_MOC_SRC tools/rigidsim-gui/replayer_multiobj.h)
    QT5_WRAP_CPP(REPLAYER_MOC_SRC tools/rigidsim-gui/replayer_multiobj.h)

    add_executable(replayer
        tools/rigidsim-gui/replayer_multiobj.cpp
        ${REPLAYER_MOC_SRC}               ${GENERIC_SRC}
        ${RIGID_SRC}                      ${LOGGING_SRC}
        ${LINEAR_SRC}                     ${IO_SRC}
        ${RIGID_IO_SRC}                   ${DEFORMABLE_SRC})
    target_link_libraries(replayer
        ${GLUT_LIBRARIES}                 ${QT_LIBRARIES}
        ${QGLViewer_LIBRARY}              ${OPENGL_LIBRARIES}
        ${LibConfig_LIBRARY})

    qt5_use_modules(replayer Xml OpenGL)

endif ( BUILD_RIGIDSIM_GUI )

# =============================================================================
#   Builds a GUI for the PAN Wavesolver
# =============================================================================
if ( BUILD_PAN_PRECOMPUTE_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    include_directories(${QT_USE_FILE})

    # Headers to be moc-ified
    set(PAN_PRECOMPUTE_GUI_MOC_HEADERS ui/WaveViewer.h)

    #QT4_WRAP_CPP(PAN_PRECOMPUTE_GUI_MOC_SRCS ${PAN_PRECOMPUTE_GUI_MOC_HEADERS})
    QT5_WRAP_CPP(PAN_PRECOMPUTE_GUI_MOC_SRCS ${PAN_PRECOMPUTE_GUI_MOC_HEADERS})

    if(USE_CUDA)
        set(CUDA_64_BIT_DEVICE_CODE ON)
        set(CUDA_NVCC_FLAGS -arch=sm_12 --maxrregcount=32 -O4 -use_fast_math -Xptxas -dlcm=ca)
        set(CUDA_NVCC_FLAGS
      "${CUDA_NVCC_FLAGS} \"-DBOOST_NOINLINE=__attribute__((noinline))\" ")
        find_library(CUDA_CUDART_LIBRARY cudart PATHS /usr/local/cuda/lib64/)

        cuda_add_executable(precompute-pan-gui
            tools/precompute_acceleration_pulse.cpp
            ui/WaveViewer.cpp
            ${PAN_PRECOMPUTE_GUI_MOC_SRCS}
            ${GENERIC_SRC}                    ${IO_SRC}
            ${LOGGING_SRC}                    ${LINEAR_SRC}
            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
            ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
            wavesolver/gpusolver/cuda/cuda_PAN_wave_3d.cu ${WAVESOLVER_CUDA_WRAPPER_SRC}
            wavesolver/gpusolver/cuda/cuda_PAT_wave_3d.cu
            ${UTILS_SRC}                      ${PARSER_SRC}
            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
            ${DEFORMABLE_SRC}                 ${SNDPLAYER_SRC}
            ${MULTIPOLE_SRC})

        target_link_libraries(precompute-pan-gui
            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
            ${GSL_LIBRARY}                  boost_thread
            boost_system                    boost_serialization
            ${MKL_LIBRARY}
            ${CUDA_LIBRARIES}               ${PORTAUDIO_LIBRARY})
    else(USE_CUDA)
        add_executable(precompute-pan-gui
            tools/precompute_acceleration_pulse.cpp
            ui/WaveViewer.cpp
            ${PAN_PRECOMPUTE_GUI_MOC_SRCS}
            ${GENERIC_SRC}                    ${IO_SRC}
            ${LOGGING_SRC}                    ${LINEAR_SRC}
            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
            ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
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
    qt5_use_modules(precompute-pan-gui Xml OpenGL)
endif ( BUILD_PAN_PRECOMPUTE_GUI )

# =============================================================================
#   Builds a GUI for the PAN Wavesolver
# =============================================================================
if ( BUILD_PAT_PRECOMPUTE_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    include_directories(${QT_USE_FILE})

    # Headers to be moc-ified
    #set(PAT_PRECOMPUTE_GUI_MOC_HEADERS ui/WaveViewer.h)

    # QT4_WRAP_CPP(PAT_PRECOMPUTE_GUI_MOC_SRCS ${PAT_PRECOMPUTE_GUI_MOC_HEADERS})

    if(USE_CUDA)
        set(CUDA_64_BIT_DEVICE_CODE ON)
        set(CUDA_NVCC_FLAGS -arch=sm_12 --maxrregcount=32 -O4 -use_fast_math -Xptxas -dlcm=ca)
        set(CUDA_NVCC_FLAGS
      "${CUDA_NVCC_FLAGS} \"-DBOOST_NOINLINE=__attribute__((noinline))\" ")
        find_library(CUDA_CUDART_LIBRARY cudart PATHS /usr/local/cuda/lib64/)

            message(STATUS "YOO1: " ${PORTAUDIO_INCLUDE_DIRS})
     message(STATUS "YOO2: " ${PORTAUDIO_LIBRARIES})
   message(STATUS "YOO3: " ${PORTAUDIO_LIBRARY_DIRS})

        cuda_add_executable(precompute-pat-gui
            tools/precompute_acoustic_transfer.cpp
            ui/WaveViewer.cpp
            ${PAN_PRECOMPUTE_GUI_MOC_SRCS}
            ${GENERIC_SRC}                    ${IO_SRC}
            ${LOGGING_SRC}                    ${LINEAR_SRC}
            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
            ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
            wavesolver/gpusolver/cuda/cuda_PAN_wave_3d.cu ${WAVESOLVER_CUDA_WRAPPER_SRC}
            wavesolver/gpusolver/cuda/cuda_PAT_wave_3d.cu
            ${UTILS_SRC}                      ${PARSER_SRC}
            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
            ${DEFORMABLE_SRC}                 ${SNDPLAYER_SRC}
            ${MULTIPOLE_SRC})

        target_link_libraries(precompute-pat-gui
            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
            ${GSL_LIBRARY}                  boost_thread
            boost_system                    boost_serialization
            ${MKL_LIBRARY}
            ${CUDA_LIBRARIES}               ${PORTAUDIO_LIBRARY})
    else(USE_CUDA)
        add_executable(precompute-pat-gui
            tools/precompute_acoustic_transfer.cpp
            ui/WaveViewer.cpp
            ${PAN_PRECOMPUTE_GUI_MOC_SRCS}
            ${GENERIC_SRC}                    ${IO_SRC}
            ${LOGGING_SRC}                    ${LINEAR_SRC}
            ${DISTANCEFIELD_SRC}              ${MATH_SRC}
            ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
            ${FIELD_SRC}                      ${WAVESOLVER_SRC}
            ${UTILS_SRC}                      ${PARSER_SRC}
            ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
            ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
            ${DEFORMABLE_SRC})

        target_link_libraries(precompute-pat-gui
            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
            ${GSL_LIBRARY}                  boost_thread
            boost_system                    ${MKL_LIBRARY})
    endif(USE_CUDA)
    qt5_use_modules(precompute-pat-gui Xml OpenGL)
endif ( BUILD_PAT_PRECOMPUTE_GUI )

# =============================================================================
#   Builds a command-line PAN wavesolver
# =============================================================================
if ( BUILD_RUN_PULSE_SOLVE )
    add_executable(run-pulse-solve
        tools/run_pulse_solve.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEAR_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC})

    target_link_libraries(run-pulse-solve
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  boost_thread
        boost_system                    ${MKL_LIBRARY})
endif ( BUILD_RUN_PULSE_SOLVE )

# =============================================================================
#   Builds a program which reads the output of run-pulse-solve or
#   precompute-pan-gui and converts it in to a multi-term PAN approximation
# =============================================================================
if ( BUILD_MULTITERM_PAN )
    add_executable(build-multiterm-pan
        tools/build_pulse_approximation_system.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEAR_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC})

    target_link_libraries(build-multiterm-pan
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  boost_thread
        boost_system                    ${MKL_LIBRARY})
endif ( BUILD_MULTITERM_PAN )

# =============================================================================
#   Builds a program which reads the output of run-pulse-solve or
#   precompute-pan-gui and converts it in to a multi-term PAN approximation
#   with wavelet-based compression
# =============================================================================
if ( BUILD_COMPRESSED_PAN )
    add_executable(build-compressed-pan
        tools/build_compressed_pulse_approximation.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEAR_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC})

    target_link_libraries(build-compressed-pan
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  boost_thread
        boost_system                    ${MKL_LIBRARY})
endif ( BUILD_COMPRESSED_PAN )

# =============================================================================
#   Builds tools for writing rigid-body data to disk (inertia and mass)
# =============================================================================
if ( BUILD_INIT_RIGID_TOOLS )
    add_executable(write-inertia
        tools/init-rigid-tools/write_inertia.cpp
        ${GENERIC_SRC}                    ${IO_SRC})

    add_executable(write-mass-center
        tools/init-rigid-tools/write_mass_center.cpp
        ${GENERIC_SRC}                    ${IO_SRC})

    add_executable(init-mass-centers
        tools/init-rigid-tools/init_mass_centers.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${GEOMETRY_SRC}                   ${LINEAR_BLAS_LAPACK_SRC})

    target_link_libraries(init-mass-centers
        ${BLAS_LAPACK_LIBS}               ${MKL_LIBRARY})
endif ( BUILD_INIT_RIGID_TOOLS )

# =============================================================================
#   Builds a tool for generating rigid-body sound, based on existing transfer
#   and impact data
# =============================================================================
if ( BUILD_RIGID_SNDGEN )
    add_executable(sound-generator
        tools/soundgen/RigidFBemSndGen.cpp
        ${SNDGEN_SRC}                     ${MATH_SRC}
        ${GENERIC_SRC}                    ${IO_SRC})

    target_link_libraries(sound-generator
        ${LibConfig_LIBRARY}              ${BLAS_LAPACK_LIBS}
        sndfile                           ${MKL_LIBRARY})
endif ( BUILD_RIGID_SNDGEN )

# =============================================================================
#   Builds a tool for estimating Hertz time scales for each impact in a
#   simulation and writing this data out to a set of files which can later
#   be used by sound-generator
# =============================================================================
if ( BUILD_ESTIMATE_TIMESCALE )
    add_executable(estimate-timescales
        tools/soundgen/estimate_timescales.cpp
        ${IO_SRC}                         ${TRANSFER_SRC}
        ${UTILS_SRC}                      ${GEOMETRY_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${PARSER_SRC}                     ${LINEAR_SRC}
        ${LINEAR_BLAS_LAPACK_SRC}         ${MATH_SRC}
        ${WAVELET_TRANSFER_SRC})

    target_link_libraries(estimate-timescales
        ${CGAL_LIBRARIES}                 ${GTS_LIBRARIES}
        ${BLAS_LAPACK_LIBS}               ${GSL_LIBRARY}
        boost_thread                      boost_system
        ${MKL_LIBRARY})
endif ( BUILD_ESTIMATE_TIMESCALE )

if (BUILD_PYTHON_EXTENSIONS)
    find_package(PythonLibs REQUIRED)
    FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)
    add_custom_target(acoustic_python)

    # Multipole
        #MultipoleMath
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(multipole_math SHARED multipole/python/math.cpp
                            multipole/MultipoleMath.h multipole/MultipoleMath.cpp)

            set_target_properties(multipole_math PROPERTIES PREFIX "_" OUTPUT_NAME multipole_math)
            TARGET_LINK_LIBRARIES(multipole_math ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)
            install(TARGETS multipole_math DESTINATION ${PROJECT_BINARY_DIR}/lib/multipole/math)
            add_dependencies(acoustic_python multipole_math)
        #MultipoleData
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(multipole_data SHARED multipole/python/data.cpp
                            multipole/MultipoleMath.h multipole/MultipoleMath.cpp
                            multipole/MultipoleData.h multipole/MultipoleData.cpp)

            set_target_properties(multipole_data PROPERTIES PREFIX "_" OUTPUT_NAME multipole_data)
            TARGET_LINK_LIBRARIES(multipole_data ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)
            install(TARGETS multipole_data DESTINATION ${PROJECT_BINARY_DIR}/lib/multipole/)
            add_dependencies(acoustic_python multipole_data)

        #MultipoleUtil
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(multipole_util SHARED multipole/python/util.cpp
                            multipole/MultipoleUtil.h multipole/MultipoleUtil.cpp
                            multipole/MultipoleMath.h multipole/MultipoleMath.cpp
                            multipole/MultipoleData.h multipole/MultipoleData.cpp)

            set_target_properties(multipole_util PROPERTIES PREFIX "_" OUTPUT_NAME multipole_util)
            TARGET_LINK_LIBRARIES(multipole_util ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)
            install(TARGETS multipole_util DESTINATION ${PROJECT_BINARY_DIR}/lib/multipole/)
            add_dependencies(acoustic_python multipole_util)
    # Mode
        #MultipoleMath
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(mode_data SHARED deformable/python/data.cpp
                            deformable/ModeData.h deformable/ModeData.cpp
                            ${CONFIG_SOURCES})

            set_target_properties(mode_data PROPERTIES PREFIX "_" OUTPUT_NAME mode_data)
            TARGET_LINK_LIBRARIES(mode_data ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)
            install(TARGETS mode_data DESTINATION ${PROJECT_BINARY_DIR}/lib/mode/)
            add_dependencies(acoustic_python mode_data)

    # Solver
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})

        cuda_add_library(pat_solver SHARED
                tools/precompute-pat/python/pat_wrapper.cpp
                ui/WaveViewer.cpp
                ${PAN_PRECOMPUTE_GUI_MOC_SRCS}
                ${GENERIC_SRC}                    ${IO_SRC}
                ${LOGGING_SRC}                    ${LINEAR_SRC}
                ${DISTANCEFIELD_SRC}              ${MATH_SRC}
                ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
                ${FIELD_SRC}                      ${WAVESOLVER_SRC}
                wavesolver/gpusolver/cuda/cuda_PAN_wave_3d.cu ${WAVESOLVER_CUDA_WRAPPER_SRC}
                wavesolver/gpusolver/cuda/cuda_PAT_wave_3d.cu
                ${UTILS_SRC}                      ${PARSER_SRC}
                ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
                ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
                ${DEFORMABLE_SRC}                 ${SNDPLAYER_SRC}
                ${MULTIPOLE_SRC})

        set_target_properties(pat_solver PROPERTIES PREFIX "_" OUTPUT_NAME solver)
        target_link_libraries(pat_solver
            ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
            ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
            ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
            ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
            ${GSL_LIBRARY}                  boost_thread
            boost_system                    boost_serialization
            ${MKL_LIBRARY}
            ${CUDA_LIBRARIES}               ${PORTAUDIO_LIBRARY}
            ${PYTHON_LIBRARIES}             ${Boost_LIBRARIES})

        qt5_use_modules(pat_solver Xml OpenGL)
        install(TARGETS pat_solver DESTINATION ${PROJECT_BINARY_DIR}/lib/solver/)
        add_dependencies(acoustic_python pat_solver)
endif (BUILD_PYTHON_EXTENSIONS)
