#ifndef SPEAKER_PRESSURE_SOURCE_H
#define SPEAKER_PRESSURE_SOURCE_H

#include <TYPES.h>
#include <config.h>
#include <linearalgebra/Vector3.hpp>
#include <wavesolver/PressureSource.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

//##############################################################################
// Speaker point pressure source
//##############################################################################
class SpeakerPressureSource : public PressureSource
{
    private:
        std::vector<REAL> _data;
        std::string       _dataFile;
        REAL              _sampleRate;
        REAL              _startTime;
        REAL              _widthSpace; // spread out using a spatial gaussian
        Vector3d          _position;
        REAL              _c = 343.0;
        const REAL        _dataThreshold = 1E-15;


    public:
        SpeakerPressureSource() = default;
        inline void SetStartTime(const REAL &t){_startTime = t;}
        inline void SetWidthSpace(const REAL &w){_widthSpace = w;}
        inline void SetSoundSpeed(const REAL &c){_c = c;}
        void Initialize(const std::string &speakerFile);
        virtual void UpdateBoundingBox();
        virtual REAL Evaluate(const Vector3d &position,
                              const Vector3d &normal,
                              const REAL &time,
                              const int &hintTriangle);
        virtual REAL Evaluate(const int &vertexID,
                              const Vector3d &vertexNormal,
                              const REAL &time);
        virtual Vector3d Evaluate(const int &vertexID,
                                  const REAL &time);
        virtual REAL EarliestEventTime(const REAL startTime) const;
        virtual bool IsZero(const REAL t,
                            const bool checkBound = false,
                            const Vector3d &minBound = Vector3d(),
                            const Vector3d &maxBound = Vector3d()) const;
};

#endif
