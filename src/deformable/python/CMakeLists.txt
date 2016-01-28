if (BUILD_PYTHON_EXTENSIONS)
    find_package(PythonLibs REQUIRED)
    FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
        ADD_LIBRARY(mode_data SHARED data.cpp
                        ${DEFORMABLE_SRC})

        set_target_properties(mode_data PROPERTIES PREFIX "_" OUTPUT_NAME mode_data)
        set_target_properties(mode_data PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}/mode)
        target_link_libraries(mode_data ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}
            ${GSL_LIBRARIES} boost_serialization)

        file(WRITE "${LIBRARY_OUTPUT_PATH}/mode/__init__.py" "from _mode_data import *\n")

        add_dependencies(acoustic_python mode_data)

endif(BUILD_PYTHON_EXTENSIONS)