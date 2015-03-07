set(AM_SOURCES "/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tools/acoustic_transfer/precompute_acoustic_transfer.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/ui/WaveViewer.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/generic/precision_type.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/io/MatrixIO.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/io/TglMeshReader.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/io/ImpulseIO.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/logging/logging.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/eig3.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/Quaternion.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/mat4inv.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/boundingBox.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/closestPointField.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/distanceField.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/mat3d.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/objfile.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/objfileOctree.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/octree.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/sphere.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/triangle.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/tribox2.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/vec2d.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/vec3d.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/FieldBuilder.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/math/InterpolationFunction.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/math/Function.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/math/SampledFunction.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/MATRIX_FAST.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/VECTOR_FAST.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/SPARSE_MATRIX.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/RigidMesh.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/field/ScalarField.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/Laplacian.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/MAC_Grid.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/PML_WaveSolver.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/WaveSolver.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/utils/timer.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/utils/MathUtil.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/parser/Parser.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tinyxml/tinystr.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tinyxml/tinyxml.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tinyxml/tinyxmlerror.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tinyxml/tinyxmlparser.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/GTS_TriMesh.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/MeshTree.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/ClosestPointMesh.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/AccelerationNoiseModel.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/MultiTermApproximation.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/PulseApproximation.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/ProxyManager.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/RadialApproximation.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/CompressedMultiTermApproximation.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/deformable/linear.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/deformable/ModeData.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/deformable/stvk.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/deformable/StVKMesh.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/multipole/MultipoleData.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/multipole/MultipoleUtil.cpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/multipole/MultipoleMath.cpp" )
set(AM_RCC_SOURCES "" )
set(AM_SKIP_MOC "" )
set(AM_SKIP_UIC "" )
set(AM_HEADERS "/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/ui/WaveViewer.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/generic/precision_type.hpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/eig3.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/Quaternion.hpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/boundingBox.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/closestPointField.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/distanceField.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/mat3d.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/minivector.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/objfile.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/objfileOctree.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/obj.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/octree.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/sphere.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/triangle.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/tribox2.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/triple.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/vec2d.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/vec3d.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/vec_types.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/distancefield/FieldBuilder.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/math/InterpolationFunction.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/math/Function.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/math/SampledFunction.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/MATRIX.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/VECTOR.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/linearalgebra/SPARSE_MATRIX.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/RigidMesh.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/field/ScalarField.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/Laplacian.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/MAC_Grid.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/PML_WaveSolver.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/WaveSolver.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/gpusolver/cuda/cuda_PAT_wave_3d.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/utils/macros.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/utils/math.hpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/utils/nano_timer.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/utils/timer.hpp;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/parser/Parser.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tinyxml/tinystr.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/tinyxml/tinyxml.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/GTS_TriMesh.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/MeshTree.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/geometry/ClosestPointMesh.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/AccelerationNoiseModel.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/MultiTermApproximation.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/PulseApproximation.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/ProxyManager.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/RadialApproximation.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/transfer/CompressedMultiTermApproximation.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/multipole/MultipoleData.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/multipole/MultipoleUtil.h;/media/jui-hsien/General/Research/gpu-wavesolver-repo/acoustics/src/multipole/MultipoleMath.h" )
set(AM_MOC_COMPILE_DEFINITIONS "CGAL_INTERSECTION_VERSION=1;DIFF_DEFINE;QT_CORE_LIB;QT_GUI_LIB;QT_NO_DEBUG;QT_OPENGL_LIB;QT_WIDGETS_LIB;QT_XML_LIB;USE_CGAL;USE_CUDA;USE_LAPACKE;USE_MKL")
set(AM_MOC_INCLUDES "/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build2/src/tools/acoustic_transfer;/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/src/tools/acoustic_transfer;/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/src;/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build2/src;/usr/local/include;/opt/intel/composer_xe_2015.0.090/mkl/include;/usr/include/glib-2.0;/usr/lib/x86_64-linux-gnu/glib-2.0/include;/usr/local/cuda-6.5/include;/usr/include/QGLViewer;/opt/Qt5.3.1/5.3/gcc_64/include;/opt/Qt5.3.1/5.3/gcc_64/include/QtXml;/opt/Qt5.3.1/5.3/gcc_64/include/QtCore;/opt/Qt5.3.1/5.3/gcc_64/mkspecs/linux-g++;/opt/Qt5.3.1/5.3/gcc_64/include/QtOpenGL;/opt/Qt5.3.1/5.3/gcc_64/include/QtWidgets;/opt/Qt5.3.1/5.3/gcc_64/include/QtGui;/usr/include")
set(AM_MOC_OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
set(AM_CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE "")
set(AM_CMAKE_BINARY_DIR "/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build2/")
set(AM_CMAKE_SOURCE_DIR "/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/")
set(AM_QT_MOC_EXECUTABLE "/opt/Qt5.3.1/5.3/gcc_64/bin/moc")
set(AM_QT_UIC_EXECUTABLE "/opt/Qt5.3.1/5.3/gcc_64/bin/uic")
set(AM_QT_RCC_EXECUTABLE "/opt/Qt5.3.1/5.3/gcc_64/bin/rcc")
set(AM_CMAKE_CURRENT_SOURCE_DIR "/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/src/tools/acoustic_transfer/")
set(AM_CMAKE_CURRENT_BINARY_DIR "/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build2/src/tools/acoustic_transfer/")
set(AM_QT_VERSION_MAJOR "5")
set(AM_TARGET_NAME "precompute-pat-gui_automoc")
set(AM_RELAXED_MODE "FALSE")
set(AM_UIC_TARGET_OPTIONS )
set(AM_UIC_OPTIONS_FILES "")
set(AM_UIC_OPTIONS_OPTIONS "")
set(AM_RCC_OPTIONS_FILES "")
set(AM_RCC_OPTIONS_OPTIONS "")
