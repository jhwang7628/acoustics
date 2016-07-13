set(IO_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/MatrixIO.cpp               
    ${CMAKE_CURRENT_SOURCE_DIR}/TglMeshReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ImpulseIO.cpp

    CACHE FILE "IO Source")

set(RIGID_IO_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidObjImpRecorder.cpp    
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidObjRecorder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidObjDispReader.h    
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidObjDispReader.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/ImpulseSeriesReader.h    
    ${CMAKE_CURRENT_SOURCE_DIR}/ImpulseSeriesReader.cpp 

    CACHE FILE "Rigid IO Source")
