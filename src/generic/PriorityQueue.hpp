/*
 * PriorityQueue.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef PRIORITY_QUEUE
#   define PRIORITY_QUEUE

#include <functional>
#include <memory>
#include <assert.h>
#include "utils/macros.h"

#ifdef USE_OPENMP_NA
#include <omp.h>
#endif

/*!
 * Data should have a field Data.qIdx
 * Compare returns true if the first argument is less than the second one
 */
template<
        typename Data, 
        typename Compare = std::less<Data>,
        typename Alloc = std::allocator<Data*>
        >
class PriorityQueue
{
    public:
        // =========== Constructors ==========
        PriorityQueue():mptr_queue(NULL), m_used(0) 
        { 
#ifdef USE_OPENMP_NA
            omp_init_lock(&m_lock);
#endif
        }
        PriorityQueue(size_t s):m_size(s), m_used(0)
        {
            mptr_queue = _alloc.allocate(s);
#ifdef USE_OPENMP_NA
            omp_init_lock(&m_lock);
#endif
        }
        ~PriorityQueue()
        {
            if ( m_size > 0 ) 
                _alloc.deallocate(mptr_queue, m_size);
#ifdef USE_OPENMP_NA
            omp_destroy_lock(&m_lock);
#endif
        }

        void  resize(size_t size);
        void  push(Data* ptr);
        //! Move up or down 
        void  update_node(Data* val);
        //! Move down along the tree
        void  move_up_node(Data* ptr);
        //! Move up along the tree
        void  move_down_node(Data* ptr);

        Data* pop();
        Data* peek()
        { return m_used ? mptr_queue[0] : NULL; }

        bool  empty() const
        {
            return m_used == 0; 
        }
        void clear() 
        {
            m_used = 0;
        }
        size_t size() const { return m_used; }

    private:
        /* 
         * These private methods are not thread-safe.
         *
         * The returned boolean indicates whether or not the queue position of 
         * the data has been changed due to the update/move operation.
         */
        bool private_move_up_node(Data* ptr);
        bool private_move_down_node(Data* ptr);

        Data**          mptr_queue;
        size_t          m_size;
        size_t          m_used;

#ifdef USE_OPENMP_NA
        omp_lock_t      m_lock;
#endif
        Compare         _comp;
        Alloc           _alloc;
};

// ==================== implementation ===================
template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::resize(size_t size)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&m_lock);
#endif
    if ( m_size > 0 ) _alloc.deallocate(mptr_queue, m_size);
    m_size = size;
    m_used = 0;
    mptr_queue = _alloc.allocate(m_size);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::push(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&m_lock);
#endif
    assert(m_used < m_size && ptr != NULL);
    mptr_queue[m_used] = ptr;
    ptr->qIdx = m_used ++;
    private_move_up_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
bool PriorityQueue<Data, Compare, Alloc>::private_move_up_node(Data* ptr)
{
    bool ret = false;
    int parent = (ptr->qIdx - 1) / 2;
    while ( parent >= 0 && _comp(*ptr, *mptr_queue[parent]) )
    {
        mptr_queue[parent]->qIdx = ptr->qIdx;
        std::swap(mptr_queue[parent], mptr_queue[ptr->qIdx]);
        ptr->qIdx = parent;
        parent = (ptr->qIdx - 1) / 2;
        ret = true;
    }
    return ret;
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::update_node(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&m_lock);
#endif
    if ( !private_move_down_node(ptr) )
        private_move_up_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::move_up_node(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&m_lock);
#endif
    private_move_up_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::move_down_node(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&m_lock);
#endif
    private_move_down_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
bool PriorityQueue<Data, Compare, Alloc>::private_move_down_node(Data* ptr)
{
    bool ret = false;
    int child = ptr->qIdx * 2 + 1;
    child += (child+1 < m_used && 
            _comp(*mptr_queue[child+1], *mptr_queue[child]));
    while ( child < m_used && _comp(*mptr_queue[child], *ptr) )
    {
        mptr_queue[child]->qIdx = ptr->qIdx;
        std::swap(mptr_queue[child], mptr_queue[ptr->qIdx]);
        ptr->qIdx = child;
        child = child * 2 + 1;
        child += (child+1 < m_used &&
                _comp(*mptr_queue[child+1], *mptr_queue[child]));
        ret = true;
    }
    return ret;
}

template<typename Data, typename Compare, typename Alloc>
Data* PriorityQueue<Data, Compare, Alloc>::pop()
{
    if ( !m_used ) return NULL;

#ifdef USE_OPENMP_NA
    omp_set_lock(&m_lock);
#endif
    Data* ret = mptr_queue[0];
    mptr_queue[0] = mptr_queue[-- m_used];
    mptr_queue[0]->qIdx = 0;

    Data* cur = mptr_queue[0];
    private_move_down_node(cur);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&m_lock);
#endif
    return ret;
}

#endif
