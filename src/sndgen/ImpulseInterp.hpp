#ifndef IMPULSE_INTERP_H
#   define IMPULSE_INTERP_H

#include <assert.h>
#include <fstream>
#include <string>
#include <vector>
#include <libconfig.h++>

#include "linearalgebra/Vector3.hpp"
#include "interp/CSpline.hpp"

#include <utils/MathUtil.h>

template<class TModalAnalysis>
class ImpulseInterp
{
    public:
        typedef CSpline<double, false>  TInterp;

        ImpulseInterp(libconfig::Setting& s, double stepsize);
        ~ImpulseInterp()
        {   clean_impulse_data(); }

        /* -------- Load impulses -------- */
        void begin_read_impulse(TModalAnalysis* pmodal)
        {
            assert(pmodal);
            mp_modal = pmodal;

            m_lastImpTs = -1E+99;
            m_impTs.clear();
            clean_impulse_data();

            m_curImp.resize(pmodal->len_eigvec() * 3);
            m_impVtx.resize(pmodal->len_eigvec());
        }
        void add_impulse(double ts, int vtxid, const Vector3d& imp, char onsurf);
        void add_modal_impulse();
        void end_read_impulse()
        {   if ( m_lastImpTs > -0.99E+99 ) add_modal_impulse(); }

        // ---------------------------------- 
        bool has_impulses() const
        {   return m_impTs.size() > 2; }
        double begin_time() const
        {   return m_impTs[0]; }
        double end_time() const
        {   return m_impTs.back(); }
        double sim_step_size() const
        {   return m_simStep; }
        int surfid( int vtxid, char onsurf ) const
        {
          int vid = onsurf == 'T' ?  (vtxid - m_numFixed):
                              (m_surfid[vtxid] - m_numFixed);
          
          return vid;
        }

        // -----------------------------------
        void init_interp(int modeid, double gwide)
        {
            m_curMId    = modeid;
            m_gaussWide = gwide*0.01;
            //m_gaussWide = gwide;

            m_yyImp.clear();
            m_lastImp = 0;
            m_vhi = m_tsid = -1;
        }

        double get_impulse(double ts);

        /* the following two functions are used in contact-damping sound synthesis */
        double impulse_ts(int id) const
        {   return m_impTs[id]; }
        void current_interval_info(int& st, int& ed) const
        {   st = m_vlo;  ed = m_vhi; }

    private:
        void update_impulse_interval();
            
        void load_geometry(const char* file);
        void clean_impulse_data()
        {
            for(size_t i = 0;i < m_modeImp.size();++ i)
                delete [](m_modeImp[i]);
            m_modeImp.clear();
        }

    private:
        int                     m_curMId;
        double                  m_gaussWide;
        int                     m_numFixed;     // # of fixed vertices
        double                  m_lastImpTs;
        // all the impulses on the modes (V^T*f)
        std::vector<double>     m_impTs;        // impulse timestamps
        std::vector< double* >  m_modeImp;      // impulse on each mode
        std::vector<double>     m_curImp;       // m_n3 x 1 vector
        std::vector<int>        m_impVtx;       // the vertices applied impulses
        int                     m_impPtr;       // pointer to the first empty slot in m_impVtx

        /* ------------ geometry ----------
         * map the vertex from surface triangle mesh to tet mesh
         * m_surfid[i] is the id in tet mesh of the i-th vertex in surface
         * triangle mesh
         */
        std::vector<int>        m_surfid;

        /* ------------ For impulse interpolation ----------- */
        double                  m_simStep;      // step size for rigid simulation
        int                     m_vlo, m_vhi, m_tsid;
        double                  m_lastImp;
        std::vector<double>     m_yyImp;        // impulse used in cspline interpolation
        TInterp                 m_impInterp;

        TModalAnalysis* mp_modal;
};

/////////////////////////////////////////////////////////////////////

template<class TModalAnalysis>
ImpulseInterp<TModalAnalysis>::ImpulseInterp(libconfig::Setting& s,
                                             double stepsize):
        m_curMId(0), m_simStep(stepsize)
{
    using namespace libconfig;

    if ( !s.lookupValue("numfixed", m_numFixed) )
    {
        PRINT_ERROR("Cannot find the property 'numfixed'\n");
        exit(3);
    }

    std::string txt;
    if ( !s.lookupValue("geometry", txt) )
    {
        PRINT_ERROR("Cannot find the property 'geometry'\n");
        exit(3);
    }
    printf("Load geometry [%s] for obj\n", txt.c_str());
    load_geometry(txt.c_str());
}

template<class TModalAnalysis>
void ImpulseInterp<TModalAnalysis>::load_geometry(const char* file)
{
    int id1, id2;
    double nx, ny, nz, a;

    std::ifstream fin(file);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file: %s\n", file);
        exit(4);
    }

    int n;
    fin >> n;   // how many surface vertices
    m_surfid.resize(n);

    for(int i = 0;i < n;++ i)
    {
        // id1: id in tetrahedron mesh
        // id2: id in surface triangle mesh
        fin >> id1 >> id2 >> nx >> ny >> nz >> a;   // both id1 and id2 are 0-based
        if ( id2 >= n )
        {
            PRINT_ERROR("id2 is out of range in geometry file\n");
            exit(4);
        }
        m_surfid[id2] = id1;
    }

    printf("%d vertices were read\n", n);
    fin.close();
}

template<class TModalAnalysis>
void ImpulseInterp<TModalAnalysis>::add_impulse(
        double ts, int vtxid, const Vector3d& imp, char onsurf)
{
    //const double IMP_SCALE = 1E+3;
    const double IMP_SCALE = 1.0;

    if ( ts > m_lastImpTs ) 
    {
        if ( m_lastImpTs > -0.99E+99 ) add_modal_impulse();
        // reset the accumulator of the current impulses
        m_impPtr = 0;
        m_lastImpTs = ts;
    }

    const int vid = onsurf == 'T' ?  (vtxid - m_numFixed):
                            (m_surfid[vtxid] - m_numFixed);    // idx in m_curImp for the vtx
    if ( vid < 0 ) return;

    if ( m_impPtr >= m_impVtx.size() ) 
    {
        ++ m_impPtr;
        m_impVtx.push_back(vid);
        m_curImp.push_back(imp.x * IMP_SCALE);
        m_curImp.push_back(imp.y * IMP_SCALE);
        m_curImp.push_back(imp.z * IMP_SCALE);
    }
    else
    {
        m_curImp[m_impPtr*3]   = imp.x * IMP_SCALE;
        m_curImp[m_impPtr*3+1] = imp.y * IMP_SCALE;
        m_curImp[m_impPtr*3+2] = imp.z * IMP_SCALE;
        m_impVtx[m_impPtr ++]  = vid;
    }
}

template<class TModalAnalysis>
void ImpulseInterp<TModalAnalysis>::add_modal_impulse()
{
    m_impTs.push_back(m_lastImpTs);
    const int NM = mp_modal->num_modes();
    double * ptr = new double[NM];
    memset(ptr, 0, sizeof(double)*NM);

    mp_modal->modal_impulses(&m_impVtx[0], &m_curImp[0], m_impPtr, ptr);

    m_modeImp.push_back(ptr);
}

/*
 * get the impulse value for the current mode at the given time
 */
template<class TModalAnalysis>
double ImpulseInterp<TModalAnalysis>::get_impulse(double ts)
{
    //// update the impulse interval
    if ( m_vhi < 0 || m_impTs[m_vhi] <= ts ) update_impulse_interval();  // update the current impulse interval

    //// update m_tsid
    if ( m_tsid < (int)m_impTs.size() ) m_tsid += (m_impTs[m_tsid+1] <= ts);

    /* Now ts should always be in [m_impTs[m_tsid], m_impTs[vhi]] */
    //// compute impuse for next tick
    if ( m_vlo >= (int)m_impTs.size() ) // ts is after the last impulse record
    {
#if 0
        printf( "gaussian = %f\n",
                MathUtil::gaussian(ts, m_impTs.back(), m_gaussWide) );
#endif
        //return MathUtil::gaussian(ts, m_impTs.back(), m_gaussWide) * m_lastImp;
        return MathUtil::gaussian_normalized(ts, m_impTs.back(), m_gaussWide) * m_lastImp;
    }
    else if ( m_tsid < 0 )              // ts is before the first impuse record
    {
#if 0
        printf( "gaussian = %f\n",
                MathUtil::gaussian(ts, m_impTs[0], m_gaussWide) );
#endif
        //return MathUtil::gaussian(ts, m_impTs[0], m_gaussWide) * m_yyImp[0];
        return MathUtil::gaussian_normalized(ts, m_impTs[0], m_gaussWide) * m_yyImp[0];
    }
    else if ( ts <= m_impTs[m_vlo] )
    {
#if 0
        printf( "gaussian1 = %f\n",
                MathUtil::gaussian(ts, m_impTs[m_vlo], m_gaussWide) );
        printf( "gaussian2 = %f\n",
                MathUtil::gaussian(ts, m_impTs[m_tsid], m_gaussWide) );
#endif
#if 0
        return MathUtil::gaussian(ts, m_impTs[m_vlo], m_gaussWide) * m_yyImp[0]  +
               MathUtil::gaussian(ts, m_impTs[m_tsid], m_gaussWide) * m_lastImp;
#endif
        return MathUtil::gaussian_normalized(ts, m_impTs[m_vlo], m_gaussWide) * m_yyImp[0]
          + MathUtil::gaussian_normalized(ts, m_impTs[m_tsid], m_gaussWide) * m_lastImp;
    }
    else
    {
        return m_impInterp.eval(ts);
    }
}

template<class TModalAnalysis>
void ImpulseInterp<TModalAnalysis>::update_impulse_interval()
{
    m_vlo = m_vhi + 1;
    for(m_vhi = m_vlo;m_vhi+1 < (int)m_impTs.size() && 
        m_impTs[m_vhi+1] - m_impTs[m_vhi] < 2.0*m_simStep;
        ++ m_vhi) ;

    if ( !m_yyImp.empty() ) m_lastImp = m_yyImp.back();
    if ( m_vlo >= (int)m_impTs.size() ) return;

    m_yyImp.resize(m_vhi - m_vlo + 1);

    for(int i = m_vlo;i <= m_vhi;++ i)
        m_yyImp[i - m_vlo] = m_modeImp[i][m_curMId];

    //// prepare the cspline 
    if ( m_vhi > m_vlo ) 
        m_impInterp.init(m_yyImp.size(), &m_impTs[m_vlo], &m_yyImp[0], 0, 0);
}

#endif

