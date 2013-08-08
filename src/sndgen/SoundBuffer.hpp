#ifndef SOUND_BUFFER_HPP
#   define SOUND_BUFFER_HPP

#include <fstream>
#include <assert.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "utils/print_msg.h"
#include "utils/macros.h"

#include "interp/InterpolationFunction.h"

class SoundBuffer
{
    public:
        SoundBuffer()
          : m_interp( NULL )
        {
        }

        /*!
         * Initialize with sample rate, and the timespane of the 
         * whole buffer
         */
        void init(int sr, double span)
        {
            m_sampleRate = sr;
            m_timeSpan   = span;
            m_buffer.assign(int(sr*span) + 1, 0);

            m_sampleStep = 1.0 / (Real)m_sampleRate;

            delete m_interp;
            m_interp = new InterpolationMitchellNetravali( m_sampleStep );
        }

        /*!
         * Begin to generated sound from time <ts>
         */
        void begin_sound_gen(double ts)
        {
            m_curIdx = int(ts * m_sampleRate);
        }

        /*!
         * Add the next sample of sound tick
         */
        bool add_sound_sample(double s)
        {
            if ( m_curIdx >= (int)m_buffer.size() ) return false;
            m_buffer[m_curIdx ++] += s;
            return true;
        }

        /*!
         * Adds a sample to an arbitrary point in time.  This requires
         * interpolation
         */
        bool add_sound_sample(double s, double t)
        {
          double             tMin, tMax;
          double             tSample;
          int                sampleMin, sampleMax;

          tMin = t - m_interp->supportLength();
          tMax = t + m_interp->supportLength();

          tMin = max( 0.0, tMin );
          tMax = min( m_timeSpan, tMax );

          sampleMin = (int)( tMin / m_sampleStep );
          sampleMax = (int)( tMax / m_sampleStep );

          // Add each of these samples
          for ( int sample_idx = sampleMin; sample_idx <= sampleMax;
                sample_idx++ )
          {
            tSample = m_sampleStep * (Real)sample_idx;

            // Evaluate interpolation function centered at t
            m_buffer[ sample_idx ]
              += s * m_interp->evaluate( tSample, t );
          }

          return true;
        }

        void scale(double s)
        {
            for(size_t i = 0;i < m_buffer.size();++ i)
                m_buffer[i] *= s;
        }

        /*!
         * normalize the current sound buffer
         */
        int normalize()
        {
            double vmax = 0;
            for(size_t i = 0;i < m_buffer.size();++ i)
                vmax = fmax(vmax, fabs(m_buffer[i]));

            if ( vmax < 1E-12 ) return ERROR_RETURN;

            vmax = 1. / (vmax + 1E-8);
            for(size_t i = 0;i < m_buffer.size();++ i)
                m_buffer[i] *= vmax;
            return SUCC_RETURN;
        }

        void write_to_raw(const char* file) const
        {
#if 0
            std::ofstream fout(file);
            for(size_t i = 0;i < m_buffer.size();++ i)
                fout << m_buffer[i] << std::endl;
            fout.flush();
            fout.close();
#endif
            FILE *f = fopen( file, "wb" );
            int size = m_buffer.size();
            fwrite( (void *)&size, sizeof( int ), 1, f );

            for ( size_t i = 0; i < m_buffer.size(); i++ )
            {
              double sample = m_buffer[ i ];
              fwrite( (void *)&sample, sizeof( double ), 1, f );
            }

            fclose( f );
        }

        int num_frames() const { return m_buffer.size(); }
        const double* data() const
        { 
            assert(!m_buffer.empty());
            return &m_buffer[0]; 
        }

        int sample_rate() const 
        { return m_sampleRate; }

    private:
        int                    m_sampleRate;
        double                 m_sampleStep;
        double                 m_timeSpan;
        int                    m_curIdx;
        std::vector<double>    m_buffer;

        InterpolationFunction *m_interp;
};

class StereoSndBuffer
{
    public:
        void init(int sr, double span, double delay=0)
        {
            m_sampleRate = sr;
            m_timeSpan   = span + delay;
            m_numFrames  = int(sr * m_timeSpan) + 1;
            m_delay      = delay;
            m_buffer.assign(m_numFrames*2, 0);
        }

        /*!
         * Begin to generated sound from time <ts>
         */
        void begin_sound_gen(double ts)
        {
            m_curIdx = int((ts + m_delay) * m_sampleRate)*2;
        }

        /*!
         * Add the next sample of sound tick
         */
        bool add_sound_sample(double ls, double rs)
        {
            if ( m_curIdx >= (int)m_buffer.size() ) return false;

            m_buffer[m_curIdx ++] += ls;
            m_buffer[m_curIdx ++] += rs;
            return true;
        }

        /*!
         * normalize the current sound buffer
         */
        int normalize()
        {
            double vmax = 0;
            for(size_t i = 0;i < m_buffer.size();++ i)
                vmax = fmax(vmax, fabs(m_buffer[i]));

            if ( vmax < 1E-12 ) return ERROR_RETURN;

            vmax = 1. / (vmax + 1E-8);
            for(size_t i = 0;i < m_buffer.size();++ i)
                m_buffer[i] *= vmax;
            return SUCC_RETURN;
        }

        int num_frames() const 
        { return m_numFrames; }

        int sample_rate() const 
        { return m_sampleRate; }

        const double* data() const
        { 
            assert(!m_buffer.empty());
            return &m_buffer[0]; 
        }

        void write_to_raw(const char* file) const
        {
            using namespace std;

            ofstream fout(file, ios::binary);
            fout.write((char *)&m_sampleRate, sizeof(int));
            fout.write((char *)&m_timeSpan,   sizeof(double));
            fout.write((char *)&m_delay,      sizeof(double));
            fout.write((char *)&m_numFrames,  sizeof(int));
            for(int i = 0;i < (int)m_buffer.size();++ i)
                fout.write((char *)&(m_buffer[i]), sizeof(double));
            fout.close();
        }

        int load_from_raw(const char *file)
        {
            using namespace std;

            ifstream fin(file, ios::binary);
            fin.read((char *)&m_sampleRate, sizeof(int));
            fin.read((char *)&m_timeSpan, sizeof(double));
            fin.read((char *)&m_delay, sizeof(double));
            fin.read((char *)&m_numFrames, sizeof(int));
            m_curIdx = 0;
            m_buffer.resize(m_numFrames*2);
            for(int i = 0;i < (int)m_buffer.size();++ i)
                fin.read((char *)&(m_buffer[i]), sizeof(double));
            if ( fin.fail() )
            {
                PRINT_WARNING("Fail to read file: %s\n", file);
                return ERROR_RETURN;
            }
            fin.close();
            return SUCC_RETURN;
        }

    private:
        int                 m_sampleRate;
        double              m_timeSpan;
        double              m_delay;
        int                 m_curIdx;
        int                 m_numFrames;
        std::vector<double> m_buffer;
};

#endif
