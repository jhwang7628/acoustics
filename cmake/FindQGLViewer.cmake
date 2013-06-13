# - Find library libQGLViewer
IF (QGLViewer_DIR)
    FIND_PATH(QGLViewer_INCLUDE_DIR QGLViewer/qglviewer.h
            PATHS ${QGLViewer_DIR}/include)
    FIND_LIBRARY(QGLViewer_LIBRARY QGLViewer
            PATHS ${QGLViewer_DIR}
            PATHS ${QGLViewer_DIR}/lib
            PATHS ${QGLViewer_DIR}/Release)
ELSE (QGLViewer_DIR)
    FIND_PATH(QGLViewer_INCLUDE_DIR QGLViewer/qglviewer.h
        PATHS ${SYSTEM_INC_PATH}
        PATHS $ENV{INCLUDE})
    FIND_LIBRARY(QGLViewer_LIBRARY QGLViewer
        PATHS ${SYSTEM_LIB_PATH}
        PATHS $ENV{LD_LIBRARY_PATH})
ENDIF (QGLViewer_DIR)

IF (QGLViewer_INCLUDE_DIR AND QGLViewer_LIBRARY)
    IF (NOT QGLViewer_FIND_QUIETLY)
        MESSAGE(STATUS "Found libQGLViewer: ${QGLViewer_LIBRARY}")
    ENDIF (NOT QGLViewer_FIND_QUIETLY)
ELSE (QGLViewer_INCLUDE_DIR AND QGLViewer_LIBRARY)
    IF (QGLViewer_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find libQGLViewer in Path:${SYSTEM_LIB_PATH}")
    ENDIF(QGLViewer_FIND_REQUIRED)
ENDIF (QGLViewer_INCLUDE_DIR AND QGLViewer_LIBRARY)

