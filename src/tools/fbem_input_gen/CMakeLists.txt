# =============================================================================
#   Builds input files for the FastBEM acoustics solver
# =============================================================================
if ( BUILD_FBEM_INPUT_GEN )
    # ===================
    # Dependencies
    # ===================


    # ===================
    # Target: fbem_input_gen
    # ===================
        add_executable(fbem_input_gen
            fbem_input_gen.cpp
            ${UTILS_SRC}          ${IO_SRC}
            ${GENERIC_SRC}
            )

        target_link_libraries(fbem_input_gen
            boost_program_options
            boost_filesystem
            boost_system
            rt)

        set_target_properties(fbem_input_gen PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_FBEM_INPUT_GEN )

# =============================================================================
#   Builds input files for the FastBEM acoustics solver
# =============================================================================
if ( BUILD_FBEM_INPUT_GEN )
    # ===================
    # Dependencies
    # ===================


    # ===================
    # Target: fieldpoints_gen
    # ===================
        add_executable(fieldpoints_gen
            fieldpoints_gen.cpp
            ${UTILS_SRC}          ${IO_SRC}
            ${GENERIC_SRC}        ${LOGGING_SRC}
            )

        target_link_libraries(fieldpoints_gen
            boost_program_options
            boost_filesystem
            boost_system
            rt)

        set_target_properties(fieldpoints_gen PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_FBEM_INPUT_GEN )