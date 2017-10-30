# # Need to find both Qt4 and QGLViewer if the QQL support is to be built
# FIND_PACKAGE(Qt4 COMPONENTS QtCore QtXml QtOpenGL QtGui REQUIRED)

if(APPLE)
    find_library(QGLViewer_LIBRARY QGLViewer DOC "QGLViewer lib for OSX")
    find_path(QGLViewer_INCLUDE_DIR QGLViewer/qglviewer.h DOC "Include for QGLViewer on OSX")
else()
    find_path(QGLViewer_INCLUDE_DIR qglviewer.h
            /usr/include/QGLViewer
            /opt/local/include/QGLViewer
            /usr/local/include/QGLViewer
            /sw/include/QGLViewer
        )

    find_library(QGLViewer_LIBRARY NAMES qglviewer-qt4 QGLViewer-qt5 QGLViewer
        PATHS
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /sw/lib
        )
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QGLVIEWER DEFAULT_MSG
    QGLViewer_INCLUDE_DIR QGLViewer_LIBRARY)
