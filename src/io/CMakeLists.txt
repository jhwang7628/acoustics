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

    CACHE FILE "Rigid IO Source")

set(IMPULSE_SERIES_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ImpulseSeriesReader.h    
    ${CMAKE_CURRENT_SOURCE_DIR}/ImpulseSeriesReader.cpp 

    CACHE FILE "Impulse Series Source")

set(BEM_IO_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/FBemReader.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FBemReader.cpp
    CACHE FILE "BEM IO Source" 
    )
