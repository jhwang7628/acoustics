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
    linearalgebra/VECTOR.h        linearalgebra/VECTOR_FAST.cpp)

set(GENERIC_SRC
    generic/precision_type.cpp)

set(LOGGING_SRC
    logging/logging.cpp)

set(DEFORMABLE_SRC
    deformable/linear.cpp         deformable/ModeData.cpp
    deformable/stvk.cpp           deformable/StVKMesh.cpp)

set(IO_SRC
    io/MatrixIO.cpp               io/TglMeshReader.cpp)

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
    math/Function.h                     math/Function.cpp)

set(GEOMETRY_SRC
    geometry/RigidMesh.h                geometry/RigidMesh.cpp)

# tinyXML-based Parser
set(PARSER_SRC
    parser/Parser.cpp                   parser/Parser.h
    tinyxml/tinystr.cpp                 tinyxml/tinystr.h
    tinyxml/tinyxml.cpp                 tinyxml/tinyxmlerror.cpp
    tinyxml/tinyxml.h                   tinyxml/tinyxmlparser.cpp)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Boost_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

# =============================================================================
#   Isostuffer GUI
# =============================================================================
if ( BUILD_ISOSTUFFER_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    include(${QT_USE_FILE})

    include_directories(
        ${GLUT_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
        ${QGLViewer_INCLUDE_DIR})

    # setup the qrc files
    set(UI_QRC_FILES ui/images.qrc)
    QT4_ADD_RESOURCES(UI_QRC_SRCS ${UI_QRC_FILES})

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
    QT4_WRAP_CPP(UI_MOC_SRCS ${UI_MOC_HEADERS} 
                 OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
    QT4_WRAP_UI(UI_UI_SRCS ${UI_UI_FILES})

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
        ${ICC_LIBS}             ${OPENGL_LIBRARIES}     rt)

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

    include(${QT_USE_FILE})
    include_directories(
        ${GLUT_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
        ${QGLViewer_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR})
    
    # =================== TETVIEWER ==================
    # the qrc files
    set(TETVIEWER_QRC_FILES tools/tetviewer/images.qrc)
    QT4_ADD_RESOURCES(TETVIEWER_QRC_SRCS ${TETVIEWER_QRC_FILES})

    QT4_WRAP_UI(TETVIEWER_UI_SRC tools/tetviewer/tetviewer.ui)
    QT4_WRAP_CPP(TETVIEWER_MOC_SRC 
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
        rt)
endif ( BUILD_ELASTICITY_SOLVER )

# =============================================================================
#   Modal eigensolver
# =============================================================================
if ( BUILD_EIGENSOLVER )
    include_directories(
        #${ARPACK_INCLUDE_DIR})
        #/home/chadwick/packages/arpack++/include)
        arpack++/include)
        
    link_directories(${MKL_LIBRARY_DIR})

    add_executable(eigensolver
        tools/eigensolver/eigensolver.cpp
        ${LOGGING_SRC}        ${UTILS_SRC}
        ${IO_SRC})
    target_link_libraries(eigensolver
        ${SUPERLU_LIBRARIES}  ${ARPACK_LIBRARY}
        ${BLAS_LAPACK_LIBS}
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
endif ( BUILD_FBEM_INPUT_GEN )

# =============================================================================
#   Builds a rigid-body solver GUI
# =============================================================================
if ( BUILD_RIGIDSIM_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    include(${QT_USE_FILE})

    # Set up the QRC files
    set(RIGIDSIM_UI_QRC_FILES tools/rigidsim-gui/images.qrc)
    QT4_ADD_RESOURCES(RIGIDSIM_UI_QRC_SRCS ${RIGIDSIM_UI_QRC_FILES})

    # Headers to be moc-ified
    set(RIGIDSIM_UI_MOC_HEADERS tools/rigidsim-gui/RigidSim.h)
    # The ui files
    set(RIGIDSIM_UI_UI_FILES tools/rigidsim-gui/RigidSim.ui)

    add_definitions(-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)

    QT4_WRAP_CPP(RIGIDSIM_UI_MOC_SRCS ${RIGIDSIM_UI_MOC_HEADERS})
    QT4_WRAP_UI(RIGIDSIM_UI_UI_SRCS ${RIGIDSIM_UI_UI_FILES})

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
        ${LibConfig_LIBRARY}              ${BLAS_LAPACK_LIBS})
endif ( BUILD_RIGIDSIM_GUI )

# =============================================================================
#   Builds a GUI for the PAN Wavesolver
# =============================================================================
if ( BUILD_PAN_PRECOMPUTE_GUI )
    # turn on the features of Qt for libQGLViewer
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

    include(${QT_USE_FILE})

    # Headers to be moc-ified
    set(PAN_PRECOMPUTE_GUI_MOC_HEADERS ui/WaveViewer.h)

    QT4_WRAP_CPP(PAN_PRECOMPUTE_GUI_MOC_SRCS ${PAN_PRECOMPUTE_GUI_MOC_HEADERS})

    add_executable(precompute-pan-gui
        tools/precompute_acceleration_pulse.cpp
        ui/WaveViewer.cpp
        ${PAN_PRECOMPUTE_GUI_MOC_SRCS}
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEAR_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEAR_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC})

    target_link_libraries(precompute-pan-gui
        ${GLUT_LIBRARIES}               ${QT_LIBRARIES}
        ${QGLViewer_LIBRARY}            ${OPENGL_LIBRARIES}
        ${BLAS_LAPACK_LIBS})
endif ( BUILD_PAN_PRECOMPUTE_GUI )
