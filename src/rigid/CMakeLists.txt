set(RIGID_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/CollisionRec.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Simulator.cpp           
    ${CMAKE_CURRENT_SOURCE_DIR}/Simulator_QS.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ContactGraph.cpp
    CACHE FILE "Rigid Sources")

set(RIGID_DATA_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidsimConfigData.h
    CACHE FILE "Rigidsim Config Data Sources")
