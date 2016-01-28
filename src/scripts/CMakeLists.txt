set(SCRIPTS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/acoustic_templates.py       ${CMAKE_CURRENT_SOURCE_DIR}/compute_multipole_coefficients.py
    ${CMAKE_CURRENT_SOURCE_DIR}/encode_frames.py
    ${CMAKE_CURRENT_SOURCE_DIR}/estimate_timescales.py      ${CMAKE_CURRENT_SOURCE_DIR}/evaluate_transfer.py
    ${CMAKE_CURRENT_SOURCE_DIR}/generate_sound.py           ${CMAKE_CURRENT_SOURCE_DIR}/make_frames.py
    ${CMAKE_CURRENT_SOURCE_DIR}/material_parameters.py      ${CMAKE_CURRENT_SOURCE_DIR}/mix_video.py
    ${CMAKE_CURRENT_SOURCE_DIR}/precompute_modes.py         ${CMAKE_CURRENT_SOURCE_DIR}/prep_transfer_eval.py 
    ${CMAKE_CURRENT_SOURCE_DIR}/run_sim.py                  ${CMAKE_CURRENT_SOURCE_DIR}/templates.py

    CACHE FILE "Scripts sources")

if(BUILD_PYTHON_EXTENSIONS)
	file(COPY ${SCRIPTS_SRC} DESTINATION ${LIBRARY_OUTPUT_PATH}/pipeline)
	file(WRITE "${LIBRARY_OUTPUT_PATH}/pipeline/__init__.py" "")

    configure_file("config_path.py.in" "${LIBRARY_OUTPUT_PATH}/pipeline/config_path.py")
    configure_file("set_lib_path.py.in" "${LIBRARY_OUTPUT_PATH}/helper/set_lib_path.py")
endif(BUILD_PYTHON_EXTENSIONS)