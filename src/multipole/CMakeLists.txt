set(MULTIPOLE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/MultipoleData.h           ${CMAKE_CURRENT_SOURCE_DIR}/MultipoleData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MultipoleUtil.h           ${CMAKE_CURRENT_SOURCE_DIR}/MultipoleUtil.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MultipoleMath.h           ${CMAKE_CURRENT_SOURCE_DIR}/MultipoleMath.cpp
    # ${CMAKE_CURRENT_SOURCE_DIR}/MultipolePlayer.h         ${CMAKE_CURRENT_SOURCE_DIR}/MultipolePlayer.cpp

    CACHE FILE "Multipole sources")

add_subdirectory(python)