#include "FBemTransferInterp.h"
#include <assert.h>
#include <string>
#include <fstream>
#include <sstream>

#include "utils/print_msg.h"
#include "utils/math.hpp"

using namespace std;

FBemTransferInterp::FBemTransferInterp(
        libconfig::Setting& s, int nmodes, 
        const vector<double>& omega,
        int objectID,
        const Vector3d *centerOfMass)
  : m_nModes(nmodes),
    m_curMId(0)
{
    using namespace libconfig;

    string txt;
    if ( !s.lookupValue("tick_file", txt) )
    {
        PRINT_ERROR("Cannot find 'tick_file' in configuration\n");
        exit(3);
    }
    load_time_ticks(txt.c_str(), centerOfMass);
    printf( "Loaded ticks\n" );

    if ( !s.lookupValue("transfer_file", txt) )
    {
        PRINT_ERROR("Cannot find 'transfer_file' in configuration\n");
        exit(3);
    }
    load_transfer(txt.c_str(), omega, objectID);
    printf( "Loaded transfer\n" );
}

void FBemTransferInterp::load_time_ticks(const char* file,
                                         const Vector3d *centerOfMass)
{
    ifstream fin(file);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot open file[%s] to read", file);
        exit(3);
    }

#if 0
    int n;
    fin >> n;
    m_transferTs.resize(n);
#endif
    m_transferTs.clear();

    double tx, ty, tz, time;
#if 0
    for(int i = 0;i < n;++ i)
        fin >> m_transferTs[i] >> tx >> ty >> tz; 
#endif
    while ( fin >> time >> tx >> ty >> tz ) {
        m_transferTs.push_back( time );

        if ( centerOfMass ) {
          Vector3d pos( tx, ty, tz );
          pos = pos - *centerOfMass;

          m_distances.push_back( pos.length() );
        }
    }

    if ( m_distances.size() > 0 ) {
      m_distanceInterp.init(m_distances.size(), &m_transferTs[0], &m_distances[0]);
    }

    fin.close();
}

void FBemTransferInterp::load_transfer(const char* fileptn,
                                       const vector<double>& omega,
                                       int objectID)
{
    double realV, imagV;
    char c;
    char file[128];

    if ( load_transfer_binary( fileptn, objectID ) ) {
      cout << "Loaded from binary" << endl;
      return;
    }

    cout << "Building transfer explicitly for object " << objectID << endl;
    cout << "Object has " << m_nModes << " modes" << endl;

    string inputData;
    string nanString = "NaN";

    double lastNorm = 0.0;

    m_transferNrm.resize(m_nModes);
    for(int i = 0, j = 0;i < m_nModes;++ i)
    {
        lastNorm = 0.0;
        m_transferNrm[i].clear();
        // FIXME
        if ( omega[i] < 125.6637 ) // || omega[i] > 121000 )
        //if ( omega[i] < 125.6637 || omega[i] > 121000 )
        {
          printf("IGNORE omega=%lf\n", omega[i]);
          m_transferNrm[i].resize( m_transferTs.size(), 0.0 );
          continue;
        }
        sprintf(file, fileptn, j++);

        ifstream fin(file);
        if ( fin.fail() )
        {
            printf( "omega = %f\n", omega[i] );
            PRINT_ERROR("Cannot open transfer function file: %s\n", file);
            m_transferNrm[i].clear();
            m_transferNrm[i].resize( m_transferTs.size() );
            //exit(3);
        }
#if 0
        while ( fin >> realV >> c >> imagV ) 
            m_transferNrm[i].push_back(sqrt(M_SQR(realV) + M_SQR(imagV)));
#endif
        while ( !fin.eof() )
        {
          getline( fin, inputData );

          if ( inputData.length() == 0 )
          {
            break;
          }

          if ( inputData.find( nanString ) != string::npos )
          {
            printf( "Warning: NaN encountered\n" );
            // Found a NaN, so cheat and copy the last value
            m_transferNrm[i].push_back( lastNorm );
          }
          else
          {
            istringstream os( inputData );

            os >> realV >> c >> imagV;
            m_transferNrm[i].push_back(sqrt(M_SQR(realV) + M_SQR(imagV)));
            lastNorm = m_transferNrm[i].back();

            // FIXME
            if ( lastNorm > 1e20 ) {
              cout << "Huge norm!!! object " << objectID << " mode " << i << endl;
            }
          }
        }
        fin.close();

        // FIXME
        if ( objectID == 140 ) {
          cout << "Fucking object 140, mode " << i << endl;
          for ( int k = 0; k < m_transferNrm[i].size(); k++ ) {
            printf( "   %i: %f\n", k, m_transferNrm[i][k] );
          }
        }

        if ( m_transferNrm[i].size() != m_transferTs.size() )
        {
            PRINT_ERROR("inconsistent number of transfer function samples m#%d [%d,%d]\n",
                    i, (int)m_transferNrm[i].size(), (int)m_transferTs.size());
            exit(3);
        }
    }

    save_transfer( fileptn, objectID );
}

/*
 * initialize the transfer interpolation
 */
void FBemTransferInterp::init_interp(int mid)
{
    m_curMId = mid;
    if ( m_transferTs.size() > 1 )
        m_transferInterp.init(m_transferTs.size(), &m_transferTs[0], &m_transferNrm[mid][0]);
}

double FBemTransferInterp::transfer_norm(double ts)
{
    assert(!m_transferTs.empty());
    if ( m_transferNrm[m_curMId].size() == 1 )
        return m_transferNrm[m_curMId][0];
    else
        return ts >= m_transferTs.back() ? m_transferNrm[m_curMId].back() : (
               ts <= m_transferTs.front() ? m_transferNrm[m_curMId].front() :
               m_transferInterp.eval(ts));
}

double FBemTransferInterp::transfer_distance(double ts)
{
  if ( m_distances.empty() ) {
    // No delays
    return 0.0;
  }

  if ( m_distances.size() == 1 ) {
    return m_distances[ 0 ];
  } else {
    return ts >= m_transferTs.back() ? m_distances.back() : (
           ts <= m_transferTs.front() ? m_distances.front() :
           m_distanceInterp.eval(ts));
  }
}

void FBemTransferInterp::save_transfer(const char* fileptn,
                                       int objectID)
{
  int size = m_transferNrm.size();

  char buf1[1024];
  char filename[1024];

  sprintf( buf1, fileptn, 0 );
  sprintf( filename, "%s_alltransfer", buf1 );

  FILE *f = fopen( filename, "wb" );
  if ( !f ) {
    cerr << "Error opening " << filename << endl;
    abort();
  }
  cout << "Saving transfer to " << filename << endl;

  fwrite( &size, sizeof(int), 1, f );
  for ( int i = 0; i < m_transferNrm.size(); i++ ) {
    size = m_transferNrm[i].size();
    fwrite( &size, sizeof(int), 1, f );

    fwrite( m_transferNrm[i].data(), sizeof(double), size, f );
  }

  fclose( f );
}

bool FBemTransferInterp::load_transfer_binary(const char* fileptn,
                                              int objectID)
{
  char buf1[1024];
  char filename[1024];

  sprintf( buf1, fileptn, 0 );
  sprintf( filename, "%s_alltransfer", buf1 );

  printf( "Trying to load binary transfer: %s\n", filename );

  FILE *f = fopen( filename, "rb" );
  if ( !f ) {
    cerr << "Could not open " << filename << endl;
    return false;
  }

  int size; 
  fread( &size, sizeof(int), 1, f );
  m_transferNrm.clear();
  m_transferNrm.resize( size );
  for ( int i = 0; i < m_transferNrm.size(); i++ ) {
    fread( &size, sizeof(int), 1, f );

    m_transferNrm[i].resize( size );

    fread( m_transferNrm[i].data(), sizeof(double), size, f );

    if ( m_transferNrm[i].size() != m_transferTs.size() )
    {
      PRINT_ERROR("inconsistent number of transfer function samples m#%d [%d,%d]\n",
                  i, (int)m_transferNrm[i].size(), (int)m_transferTs.size());
      exit(3);
    }
  }

  fclose( f );

  return true;
}
