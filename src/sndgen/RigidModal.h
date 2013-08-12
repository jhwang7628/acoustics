#ifndef RIGID_MODAL_H
#   define RIGID_MODAL_H
/*
 * Load rigid model analysis from file
 */
#include <list>
#include <algorithm>
#include <vector>
#include <libconfig.h++>

#include "linearalgebra/Vector3.hpp"

class RigidModal
{
    public:
        RigidModal(libconfig::Setting& s);


        void modal_impulses(const int* impvtx, const double* imp, 
                            int impLen, double* modalImp) const;

        /* ------- getter methods ------- */
        int num_modes() const
        {   return m_numModes; }
        int len_eigvec() const  // length of eigen vector
        {   return m_n3; }
        int id() const
        {   return m_id; }
        double inv_density() const
        {   return m_invDensity; }
        const std::vector<double>& omega() const
        {   return m_omega; }
        const std::vector<double>& omegad() const
        {   return m_omegad; }
        const std::vector<double>& xi() const
        {   return m_xi; }
        const std::vector<double>& gaussian_wide() const
        {   return m_wides; }
        const std::vector<double>& eigenvec() const
        {   return m_eigenvec; }    // eigenvector in column order
        double alpha() const
        {   return m_alpha; }
        double beta() const
        {   return m_beta; }
        double gamma() const
        {   return m_gamma; }

    private:
        void load_eigenmodes(const char* file);

    protected:
        double      m_density;
        double      m_invDensity;
        double      m_alpha;
        double      m_beta;
        double      m_gamma;

        int         m_n3;
        int         m_numModes;

        int         m_id;

        /* v1_1 v1_2 ... v1_n v2_1 v2_2 ... */
        std::vector<double>     m_eigenvec;     // eigenvector in column order
        std::vector<double>     m_eigenmodes;

        std::vector<double>     m_omega;
        std::vector<double>     m_omegad;
        std::vector<double>     m_xi;
        std::vector<double>     m_wides;        // wide for gaussian impulse lobe

};

struct ImpulseRecord {
  ImpulseRecord()
  {
  }

  ImpulseRecord( int vertexID, double impulseStart, double impulseLength,
                 double impulseScale, const Vector3d &impulse )
    : m_impulseVertex( vertexID ),
      m_impulseStartTime( impulseStart ),
      m_impulseLength( impulseLength ),
      m_impulseScale( impulseScale )
  {
    double               impulseNorm = impulse.norm();

    m_impulseDirection[ 0 ] = impulse[ 0 ] / impulseNorm;
    m_impulseDirection[ 1 ] = impulse[ 1 ] / impulseNorm;
    m_impulseDirection[ 2 ] = impulse[ 2 ] / impulseNorm;
  }

  ImpulseRecord( const ImpulseRecord &rec )
    : m_impulseVertex( rec.m_impulseVertex ),
      m_impulseStartTime( rec.m_impulseStartTime ),
      m_impulseLength( rec.m_impulseLength ),
      m_impulseScale( rec.m_impulseScale )
  {
    m_impulseDirection[ 0 ] = rec.m_impulseDirection[ 0 ];
    m_impulseDirection[ 1 ] = rec.m_impulseDirection[ 1 ];
    m_impulseDirection[ 2 ] = rec.m_impulseDirection[ 2 ];
  }

  ImpulseRecord &operator=( const ImpulseRecord &rec )
  {
    m_impulseVertex = rec.m_impulseVertex;
    m_impulseStartTime = rec.m_impulseStartTime;
    m_impulseLength = rec.m_impulseLength;
    m_impulseScale = rec.m_impulseScale;

    m_impulseDirection[ 0 ] = rec.m_impulseDirection[ 0 ];
    m_impulseDirection[ 1 ] = rec.m_impulseDirection[ 1 ];
    m_impulseDirection[ 2 ] = rec.m_impulseDirection[ 2 ];

    return *this;
  }

  int                    m_impulseVertex;
  double                 m_impulseStartTime;
  double                 m_impulseLength;
  double                 m_impulseScale;

  double                 m_impulseDirection[ 3 ];

  bool operator<( const ImpulseRecord &otherRecord )
  {
    return m_impulseStartTime < otherRecord.m_impulseStartTime;
  }
};

typedef std::vector<ImpulseRecord>   ImpulseTable;
typedef std::list<ImpulseRecord>     ImpulseList;

template <class TModalAnalysis>
class HertzImpulseSeries
{
  public:
    HertzImpulseSeries( TModalAnalysis *rigidModal );

    ~HertzImpulseSeries();

    void add_impulse( int vertexID, double impulseStart, double impulseLength,
                      double impulseScale, const Vector3d &impulse );

    void init();

    void initInterp( int modeID );

    double getForce( double t );

    void clearImpulses()
    {
      m_impulses.clear();
    }

    void save( const char *filename );
    bool load( const char *filename );

    double startTime() const
    {
      return m_startTime;
    }

    double endTime() const
    {
      return m_endTime;
    }

  private:
    // Updates the set of active impulses to reflect the current
    // position in time
    void updateActiveImpulses( double t );

    void addNewImpulses( double t );
    void removeInactiveImpulses( double t );

  private:
    TModalAnalysis              *m_rigidModal; 

    ImpulseTable             m_impulses;

    // We will walk through the impulses sequentially in time
    // and keep track of the current position
    ImpulseList              m_activeImpulses;
    int                      m_activePosition;

    int                      m_currentModeID;

    double                   m_startTime;
    double                   m_endTime;

};

static bool operator<( const ImpulseRecord &r1, const ImpulseRecord &r2 )
{
  return r1.m_impulseStartTime < r2.m_impulseStartTime;
}

/////////////////////////////////////////////////////////////////////
// HertzImpulseSeries implementation
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
HertzImpulseSeries<TModalAnalysis>::HertzImpulseSeries(
                                              TModalAnalysis *rigidModal )
  : m_rigidModal( rigidModal ),
    m_startTime( 0.0 ),
    m_endTime( 0.0 )
{
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
HertzImpulseSeries<TModalAnalysis>::~HertzImpulseSeries()
{
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::add_impulse( int vertexID,
                                                      double impulseStart,
                                                      double impulseLength,
                                                      double impulseScale,
                                                      const Vector3d &impulse )
{
#if 0
  printf( "Adding impulse for object with start time %f, "
          "length %f and scale %f\n",
          impulseStart, impulseLength, impulseScale );
#endif
  m_impulses.push_back( ImpulseRecord( vertexID, impulseStart,
                                       impulseLength, impulseScale,
                                       impulse ) );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::init()
{
#if 0
  for ( int i = 0; i < m_impulses.size(); i++ )
  {
    printf( "Old start time = %f, length = %f, scale = %f\n",
            m_impulses[ i ].m_impulseStartTime,
            m_impulses[ i ].m_impulseLength,
            m_impulses[ i ].m_impulseScale );
  }
#endif
  // Makes sure the impulses are sorted temporally
  std::sort( m_impulses.begin(), m_impulses.end() );

  for ( int i = 0; i < m_impulses.size(); i++ )
  {
    m_startTime = std::min( m_startTime, m_impulses[ i ].m_impulseStartTime );
    m_endTime = std::max( m_endTime,
                          m_impulses[ i ].m_impulseStartTime
                            + m_impulses[ i ].m_impulseLength );
#if 0
    printf( "Impulse start time = %f, length = %f, scale = %f\n",
            m_impulses[ i ].m_impulseStartTime,
            m_impulses[ i ].m_impulseLength,
            m_impulses[ i ].m_impulseScale );
#endif
  }

  printf( "Hertz Impulse series has %d impulses\n", (int)m_impulses.size() );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::initInterp( int modeID )
{
  m_currentModeID = modeID;
  m_activeImpulses.clear();
  m_activePosition = 0;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
double HertzImpulseSeries<TModalAnalysis>::getForce( double t )
{
  updateActiveImpulses( t );

  // Get the eigenvectors
  const std::vector<double> &eigenvectors = m_rigidModal->eigenvec();
  int                        length = m_rigidModal->len_eigvec();
  const double              *modeShape;
  double                     modeForce = 0.0;
  double                     impulseForce;
  double                     impulseScale;

  ImpulseList::iterator      itr = m_activeImpulses.begin();
  
  modeShape = &eigenvectors[ m_currentModeID * length ];

  for ( ; itr != m_activeImpulses.end(); itr++ )
  {
    const ImpulseRecord     &rec = *itr;
    int                      vid = rec.m_impulseVertex;

    impulseForce = 0.0;

    impulseForce += rec.m_impulseDirection[ 0 ] * modeShape[ 3 * vid + 0 ];
    impulseForce += rec.m_impulseDirection[ 1 ] * modeShape[ 3 * vid + 1 ];
    impulseForce += rec.m_impulseDirection[ 2 ] * modeShape[ 3 * vid + 2 ];

#if 0
    printf( "impulseForce = %f\n", impulseForce );
#endif

    // Scale by a half-sine pulse
    impulseScale = t - rec.m_impulseStartTime;
    impulseScale *= M_PI / rec.m_impulseLength;
    impulseScale = sin( impulseScale );
    impulseScale *= rec.m_impulseScale;
    impulseScale *= m_rigidModal->inv_density();

#if 0
    printf( "rec.m_impulseScale = %f\n", rec.m_impulseScale );
    printf( "impulseScale = %f\n", impulseScale );
#endif

    modeForce += impulseForce * impulseScale;
  }

  return modeForce;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::save( const char *filename )
{
  int                        size = m_impulses.size();

  FILE                      *f = fopen( filename, "wb" );

  if ( !f ) {
    std::cerr << "Failed to open " << filename << " for writing" << std::endl;
    abort();
  }

  fwrite( &size, sizeof( int ), 1, f );
  fwrite( m_impulses.data(), sizeof( ImpulseRecord ), size, f );

  fclose( f );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
bool HertzImpulseSeries<TModalAnalysis>::load( const char *filename )
{
  int                        size;

  FILE                      *f = fopen( filename, "rb" );

  if ( !f ) {
    std::cerr << "Couldn't read from " << filename << std::endl;
    return false;
  }
  std::cout << "Loading impulses from " << filename << std::endl;

  fread( &size, sizeof( int ), 1, f );
  m_impulses.clear();
  m_impulses.resize( size );
  fread( m_impulses.data(), sizeof( ImpulseRecord ), size, f );

  fclose( f );

  return true;
}

/////////////////////////////////////////////////////////////////////
// Updates the set of active impulses to reflect the current
// position in time
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::updateActiveImpulses( double t )
{
  addNewImpulses( t );
  removeInactiveImpulses( t );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::addNewImpulses( double t )
{
  // Add all impulses starting before time t
  while ( m_activePosition < m_impulses.size()
       && m_impulses[ m_activePosition ].m_impulseStartTime < t )
  {
#if 0
    printf( "adding impulse with start time%f\n",
            m_impulses[ m_activePosition ].m_impulseStartTime );
#endif
    m_activeImpulses.push_back( m_impulses[ m_activePosition ] );
    m_activePosition++;
  }
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
template <class TModalAnalysis>
void HertzImpulseSeries<TModalAnalysis>::removeInactiveImpulses( double t )
{
  // Remove any impulses that have ended by time t
  ImpulseList::iterator      itr = m_activeImpulses.begin();

  while ( itr != m_activeImpulses.end() )
  {
    if ( itr->m_impulseStartTime + itr->m_impulseLength <= t )
    {
      // We're done with this impulse
      itr = m_activeImpulses.erase( itr );
    }
    else
    {
      itr++;
    }
  }
}

#endif
