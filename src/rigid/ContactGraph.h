#ifndef RIGID_CONTACT_GRAPH_H
#   define RIGID_CONTACT_GRAPH_H

#include <vector>
#include "geometry/FixVtxTetMesh.hpp"
#include "rigid/LSCollisionRigidBody.hpp"

class ContactGraph
{
    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;
        typedef TRigidBody::TCollProc                   TCollProc;
        typedef CollisionConstraint<REAL, TCollProc>    TConstraint;

    public:
        ContactGraph():m_cmp(this) { }

        /*!
         * initialize the contact graph. 
         * NOTE: This method should be called each time some old rigid objects
         *       are removed or some new ones are added.
         */
        void init(const std::vector<TRigidBody*>& rigids,
                  const std::vector<TConstraint*>& cons);

        void update();

        void print_graph() const; 

        const std::vector<int>& levels() const
        { return m_level; }
        const std::vector<int>& sorted_objs() const
        { return m_sortedObjs; }
        const std::vector<int>& root_objs() const
        { return m_root; }

    private:
        void update_level();

    private:
        struct LevelCmp__
        {
            const ContactGraph* graph;

            LevelCmp__(ContactGraph* g):graph(g) { }

            bool operator () (int a, int b) const
            {
                return graph->m_level[a] == graph->m_level[b] ? 
                        graph->m_priori[a] < graph->m_priori[b] : 
                        graph->m_level[a]  < graph->m_level[b]; 
            }
        };

        std::vector<TRigidBody*>           m_rigids;
        /* currently we only consider the fixed constraint */
        std::vector<TConstraint*>          m_fixedConstraints;

        std::vector<int>                   m_priori;
        /* the set of objs contacting with some fixed constraint */
        std::vector<int>                   m_level;     // level of the i-th obj
        std::vector< std::vector<int> >    m_graph;
        std::vector<int>                   m_root;
        /* the IDs of the rigid objs sorted based on their level */
        std::vector<int>                   m_sortedObjs;
        LevelCmp__                         m_cmp; 
};

#endif
