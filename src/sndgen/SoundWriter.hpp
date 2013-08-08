/*
 * SoundWriter.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef SOUND_WRITER_HPP
#   define SOUND_WRITER_HPP

#include <sndfile.h>
#include <assert.h>

/*!
 * Write the data into a single channel WAV file
 */
class SingleChannelWavWriter
{
    public:
        /*!
         * \param file  name of the output file
         * \param sr    sample rate
         * \param BIT   bit per frame for output file
         */
        SingleChannelWavWriter(const char* file, int sr, int BIT)
        {
            memset(&m_sfInfo, 0, sizeof(m_sfInfo));
            m_sfInfo.samplerate = sr;
            m_sfInfo.channels = 1;
            m_sfInfo.format = SF_FORMAT_WAV; 

            switch ( BIT )
            {
                case 8:
                    m_sfInfo.format |= SF_FORMAT_PCM_S8;
                    break;
                case 16:
                    m_sfInfo.format |= SF_FORMAT_PCM_16;
                    break;
                case 24:
                    m_sfInfo.format |= SF_FORMAT_PCM_24;
                    break;
                case 32:
                    m_sfInfo.format |= SF_FORMAT_PCM_32;
                    break;
                case 64:
                    m_sfInfo.format |= SF_FORMAT_DOUBLE;
                    break;
                default:
                    fprintf(stderr, "Unsupported bit rate! %d\n", BIT);
                    exit(1);
            }

            m_sndFile = sf_open(file, SFM_WRITE, &m_sfInfo);
            assert(m_sndFile);
        }

        ~SingleChannelWavWriter()
        {
            if ( m_sndFile )
                sf_close(m_sndFile);
        }

        template<typename T>
        int write(const T* data, int numFrame);

        template<typename T>
        void set_normalized(bool normalizeIt);

    private:
        SF_INFO     m_sfInfo;
        SNDFILE*    m_sndFile;
};

template<>
int SingleChannelWavWriter::write(const double* data, int numFrame)
{
    return sf_write_double(m_sndFile, data, numFrame);
}

template<>
void SingleChannelWavWriter::set_normalized<double>(bool normalizeIt)
{
    sf_command(m_sndFile, SFC_SET_NORM_DOUBLE, NULL, normalizeIt);
}


/********************************************************************/
/**
 * Writer two sound buffer (or stereo sound buffer) into a Wav file
 * as a stereo sound
 */
class StereoWavWriter
{
    public:
        /*!
         * \param file  name of the output file
         * \param sr    sample rate
         * \param BIT   bit per frame for output file
         */
        StereoWavWriter(const char* filename, int sr, int bit)
        {
            memset(&m_sfInfo, 0, sizeof(m_sfInfo));
            m_sfInfo.samplerate = sr;
            m_sfInfo.channels   = 2;
            m_sfInfo.format     = SF_FORMAT_WAV; 

            switch ( bit )
            {
                case 8:
                    m_sfInfo.format |= SF_FORMAT_PCM_S8;
                    break;
                case 16:
                    m_sfInfo.format |= SF_FORMAT_PCM_16;
                    break;
                case 24:
                    m_sfInfo.format |= SF_FORMAT_PCM_24;
                    break;
                case 32:
                    m_sfInfo.format |= SF_FORMAT_PCM_32;
                    break;
                default:
                    fprintf(stderr, "Unsupported bit rate! %d\n", bit);
                    exit(1);
            }

            m_sndFile = sf_open(filename, SFM_WRITE, &m_sfInfo);
            assert(m_sndFile);
        }

        ~StereoWavWriter()
        {
            if ( m_sndFile ) sf_close(m_sndFile);
        }

        template<typename T>
        int write(const T* data1, const T* data2, int numFrame)
        {
            std::vector<T> buf(numFrame*2);
            for(int i = 0;i < numFrame;++ i)
            {
                buf[i*2]   = data1[i];
                buf[i*2+1] = data2[i];
            }
            return sf_writef_double(m_sndFile, &buf[0], numFrame);
        }

        template<typename T>
        int write(const T* data, int numFrame)
        {
            return sf_writef_double(m_sndFile, data, numFrame);
        }

        template<typename T>
        void set_normalized(bool normalizeIt)
        {
            sf_command(m_sndFile, SFC_SET_NORM_DOUBLE, NULL, normalizeIt);
        }

    private:
        SF_INFO     m_sfInfo;
        SNDFILE*    m_sndFile;
};
#endif
