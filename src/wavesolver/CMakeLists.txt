set(WAVESOLVER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/Laplacian.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Laplacian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MAC_Grid.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MAC_Grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PML_WaveSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/PML_WaveSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PML_WaveSolver_Settings.h
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PressureSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussianPressureSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussianPressureSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VibrationalSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/HarmonicVibrationalSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/HarmonicVibrationalSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ModalVibrationalSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ModalVibrationalSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationNoiseVibrationalSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationNoiseVibrationalSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WaterVibrationalSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/WaterVibrationalSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_Objects.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_Objects.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_MovableObject.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_MovableObject.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidObject.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidObject.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_ShellObject.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_ShellObject.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidObject_Animator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidObject_Animator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidSoundObject.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_RigidSoundObject.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_AcousticSimulator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_AcousticSimulator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Wavesolver_ConstantsAndTypes.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SimWorld.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SimWorld.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SimWorldAuxDS.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/AudioOutput.h
    ${CMAKE_CURRENT_SOURCE_DIR}/AudioOutput.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryInterface.h
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryInterface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_PlaneConstraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FDTD_PlaneConstraint.cpp

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
