if (BUILD_PYTHON_EXTENSIONS)
    find_package(PythonLibs REQUIRED)
    FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

    # Multipole
        #MultipoleMath
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(multipole_math SHARED math.cpp
                            ${MULTIPOLE_SRC})

            set_target_properties(multipole_math PROPERTIES PREFIX "_" OUTPUT_NAME multipole_math)
            set_target_properties(multipole_math PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}/multipole/math)
            target_link_libraries(multipole_math ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)

            file(WRITE "${LIBRARY_OUTPUT_PATH}/multipole/math/__init__.py" "from _multipole_math import *\n")

            add_dependencies(acoustic_python multipole_math)
        #MultipoleData
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(multipole_data SHARED data.cpp
                            ${MULTIPOLE_SRC})

            set_target_properties(multipole_data PROPERTIES PREFIX "_" OUTPUT_NAME multipole_data)
            set_target_properties(multipole_data PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}/multipole/)
            target_link_libraries(multipole_data ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)

            add_dependencies(acoustic_python multipole_data)

        #MultipoleUtil
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
            ADD_LIBRARY(multipole_util SHARED util.cpp
                            ${MULTIPOLE_SRC})

            set_target_properties(multipole_util PROPERTIES PREFIX "_" OUTPUT_NAME multipole_util)
            set_target_properties(multipole_util PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}/multipole/)
            target_link_libraries(multipole_util ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
                ${GSL_LIBRARIES} boost_serialization)
            
            add_dependencies(acoustic_python multipole_util)
            file(WRITE "${LIBRARY_OUTPUT_PATH}/multipole/__init__.py" "from _multipole_data import *\nfrom _multipole_util import *\n")
endif (BUILD_PYTHON_EXTENSIONS)
