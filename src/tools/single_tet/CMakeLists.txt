# =============================================================================
#   Simple utility for generating a tet mesh with one tet
# =============================================================================
if ( BUILD_SINGLE_TET )
    # ===================
    # Dependencies
    # ===================


    # ===================
    # Target: single_tet
    # ===================
        add_executable(single_tet
            single_tet.cpp
            ${UTILS_SRC}          ${IO_SRC}
            ${LINEARALGEBRA_SRC}  ${GENERIC_SRC}
            )

        target_link_libraries(single_tet
            rt)

        set_target_properties(single_tet PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_SINGLE_TET )