#include "ContactGraph.h"
#include "rigid/LSCollisionDetect.hpp"
#include <queue>

void ContactGraph::init(const std::vector<TRigidBody*>& rigids,
                        const std::vector<TConstraint*>& cons)
{
    m_rigids.resize(rigids.size());
    m_priori.resize(rigids.size());
    m_sortedObjs.resize(rigids.size());
    m_level.resize(rigids.size());
    m_graph.resize(rigids.size());

    for(size_t i = 0;i < rigids.size();++ i)
    {
        m_rigids[i] = rigids[i];
        m_priori[i] = i;
        m_sortedObjs[i] = i;
    }

    m_fixedConstraints.clear();
    for(size_t i = 0;i < cons.size();++ i)
        if ( cons[i]->is_fixed() ) m_fixedConstraints.push_back(cons[i]);
}

void ContactGraph::update()
{
    //// shuffle the priority
    int tida = rand() % m_priori.size();
    int tidb = rand() % m_priori.size();
    if ( tida != tidb ) std::swap(m_priori[tida], m_priori[tidb]);

    m_root.clear();

    //// detect the objs contacted with fixed constraints
    //   m_root stores id of the obj that is colliding with and fixed
    //   constraints
    for(size_t objId = 0;objId < m_rigids.size();++ objId)
    {
        if ( m_rigids[objId]->is_fixed() ) continue;

        for(size_t cId = 0;cId < m_fixedConstraints.size();++ cId)
            if ( m_fixedConstraints[cId]->is_colliding(
                        m_rigids[objId]->collision_processor(),
                        m_rigids[objId]->collision_processor()->bounding_tree_root()) )
            {
                m_root.push_back(objId);
                break;
            }
    }

    for(size_t i = 0;i < m_rigids.size();++ i) m_graph[i].clear();

    TCollProc *cp1, *cp2;
    //// detect each pair of objs that is colliding
    for(size_t i = 0;i < m_rigids.size();++ i)
    for(size_t j = i+1;j < m_rigids.size();++ j)
    {
        cp1 = m_rigids[i]->collision_processor();
        cp2 = m_rigids[j]->collision_processor();

        cp1->m_collidedVtx.clear();
        cp2->m_collidedVtx.clear();

        //// get the candidate vertices for collision detection
        detect_tree_node_collisions(
                cp1->bounding_tree_root(),
                cp2->bounding_tree_root(),
                cp1, cp2);

        if ( cp1->is_colliding(m_rigids[j]) || cp2->is_colliding(m_rigids[i]) )
        {
            m_graph[i].push_back(j);
            m_graph[j].push_back(i);
        }
    }

    //// update level number for each object
    update_level();
    std::sort(m_sortedObjs.begin(), m_sortedObjs.end(), m_cmp);
}

/*!
 * update the level value for each object using BFS
 */
void ContactGraph::update_level()
{
    std::queue<int>  objQueue;

    for(size_t i = 0;i < m_rigids.size();++ i) 
    {
        if ( m_rigids[i]->is_fixed() )
        {
            m_level[i] = -1;
            objQueue.push(i);
        }
        else
            m_level[i] = 10000000;
    }

    for(size_t i = 0;i < m_root.size();++ i) 
        if ( m_level[m_root[i]] > 0 )
        {
            m_level[m_root[i]] = 0;
            objQueue.push(i);
        }

    while ( !objQueue.empty() )
    {
        int objId = objQueue.front();
        objQueue.pop();

        for(size_t i = 0;i < m_graph[objId].size();++ i)    // iterate all the neighbors
        {
            int clevel = m_level[objId] + 1;
            if ( clevel < m_level[m_graph[objId][i]] )
            {
                m_level[m_graph[objId][i]] = clevel;
                objQueue.push(m_graph[objId][i]);
            }
        }
    } // end while
}

void ContactGraph::print_graph() const
{
    int last = -100;
    for(size_t i = 0;i < m_level.size();++ i)
    {
        int lvl = m_level[m_sortedObjs[i]];
        if ( lvl != last )
        {
            if ( last >= 0 ) printf(")  ");
            last = lvl;
            printf("([%d]%d", lvl, m_sortedObjs[i]);
        }
        else
            printf(", %d", m_sortedObjs[i]);
    }
    printf(")\n");
}

