set(WAVESOLVER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/Laplacian.cpp      ${CMAKE_CURRENT_SOURCE_DIR}/Laplacian.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MAC_Grid.cpp       ${CMAKE_CURRENT_SOURCE_DIR}/MAC_Grid.h
    ${CMAKE_CURRENT_SOURCE_DIR}/PML_WaveSolver.cpp ${CMAKE_CURRENT_SOURCE_DIR}/PML_WaveSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveSolver.cpp     ${CMAKE_CURRENT_SOURCE_DIR}/WaveSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidObject.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidObject.h


    CACHE FILE "Wave solver sources")

set(WAVESOLVER_CUDA_PAN_WRAPPER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.cpp

    CACHE FILE "CUDA Wave solver wrapper sources")

set(WAVESOLVER_CUDA_PAT_WRAPPER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.cpp

    CACHE FILE "CUDA Wave solver wrapper sources")

set(WAVESOLVER_CUDA_PAN_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAN_wave_3d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAN_wave_3d.cu
    # ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAN_wave_3d_kernel.cu
    # ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAN_wave_3d_kernel_math.cu

    CACHE FILE "CUDA Wave solver sources" FORCE)

set(WAVESOLVER_CUDA_PAT_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAT_wave_3d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAT_wave_3d.cu
    # ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAT_wave_3d_kernel.cu
    # ${CMAKE_CURRENT_SOURCE_DIR}/gpusolver/cuda/cuda_PAT_wave_3d_kernel_math.cu

    CACHE FILE "CUDA Wave solver sources" FORCE)
