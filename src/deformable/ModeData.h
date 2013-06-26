///////////////////////////////////////////////////////////////////////////////
// ModeData.h: Simple data structure for storing modal displacements and
//             frequencies
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __MODE_DATA_H__
#define __MODE_DATA_H__

#include <config.h>

#include <vector>

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
struct ModeData {
    public:
        // Eigen-values produced by modal analysis
        std::vector<REAL>                    _omegaSquared;

        std::vector<std::vector<REAL> >      _modes;

    public:
        std::vector<REAL>       &mode( int modeIndex )
                                 {
                                    return _modes[ modeIndex ];
                                 }

        const std::vector<REAL> &mode( int modeIndex ) const
                                 {
                                    return _modes[ modeIndex ];
                                 }

        REAL                     omegaSquared( int modeIndex ) const
                                 {
                                    return _omegaSquared[ modeIndex ];
                                 }

        int                      numModes() const
                                 {
                                     return _omegaSquared.size();
                                 }

        int                      numDOF() const
                                 {
                                     return ( numModes() > 0 ) ? _modes[ 0 ].size() : 0;
                                 }

        void                     read( const char *filename );
        void                     write( const char *filename ) const;

};

#endif // __MODE_DATA_H__
