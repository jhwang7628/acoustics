set(TRANSFER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationNoiseModel.h   ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationNoiseModel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MultiTermApproximation.h   ${CMAKE_CURRENT_SOURCE_DIR}/MultiTermApproximation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PulseApproximation.h       ${CMAKE_CURRENT_SOURCE_DIR}/PulseApproximation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ProxyManager.h             ${CMAKE_CURRENT_SOURCE_DIR}/ProxyManager.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RadialApproximation.h      ${CMAKE_CURRENT_SOURCE_DIR}/RadialApproximation.cpp

    CACHE FILE "Transfer source")

set(WAVELET_TRANSFER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/CompressedMultiTermApproximation.h
    ${CMAKE_CURRENT_SOURCE_DIR}/CompressedMultiTermApproximation.cpp

    CACHE FILE "Wavelet Transfer source")