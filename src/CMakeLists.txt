## CMake specification for src directory
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Boost_INCLUDE_DIR}
    ${Eigen_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR})

link_directories(${Boost_LIBRARY_DIRS})

if(BUILD_PYTHON_EXTENSIONS)
    add_custom_target(acoustic_python)
endif(BUILD_PYTHON_EXTENSIONS)

set(CONFIG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/config.h
    CACHE FILE "Config sources")

add_subdirectory(arpack++)
add_subdirectory(isostuffer)
add_subdirectory(utils)
add_subdirectory(linearalgebra)
add_subdirectory(generic)
add_subdirectory(logging)
add_subdirectory(deformable)
add_subdirectory(io)
add_subdirectory(rigid)
add_subdirectory(field)
add_subdirectory(graph)

add_subdirectory(wavesolver)

add_subdirectory(superlu-interface)
add_subdirectory(eigensolver)

# Sources for Jernej Barbic's flat signed distance function classes
add_subdirectory(distancefield)
add_subdirectory(math)
add_subdirectory(geometry)
add_subdirectory(transfer)

# For regular grid finite difference 
add_subdirectory(grid) 

# tinyXML-based Parser
set(PARSER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/parser/Parser.cpp                   ${CMAKE_CURRENT_SOURCE_DIR}/parser/Parser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/tinystr.cpp                 ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/tinystr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/tinyxml.cpp                 ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/tinyxmlerror.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/tinyxml.h                   ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/tinyxmlparser.cpp

    CACHE FILE "Parser sources")

# Sound generation sources
add_subdirectory(sndgen)
add_subdirectory(multipole)
add_subdirectory(scripts)

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

add_subdirectory(tools)
