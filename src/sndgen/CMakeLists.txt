set(SNDGEN_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/FBemTransferInterp.h         ${CMAKE_CURRENT_SOURCE_DIR}/FBemTransferInterp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ImpulseInterp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidModal.h                 ${CMAKE_CURRENT_SOURCE_DIR}/RigidModal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RigidSoundObj.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SoundBuffer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SoundWriter.hpp

    CACHE FILE "SoundGen Sources")

set(SND_IO_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/WavReader.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WavReader.cpp

    CACHE FILE "Sound IO Sources")

set(SNDPLAYER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/player/SoundPlayer.h         ${CMAKE_CURRENT_SOURCE_DIR}/player/SoundPlayer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/player/SineWavePlayer.h

    CACHE FILE "SoundPlayer Sources"
    )
