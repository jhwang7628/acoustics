/*
 * =====================================================================================
 *
 *       Filename:  PQTests.cpp
 *
 *        Version:  1.0
 *        Created:  12/06/10 23:56:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include <gtest/gtest.h>
#include "generic/PriorityQueue.hpp"

using namespace std;

struct PQEle
{
    int qIdx;
    double u;

    PQEle(double u):u(u) { }

    bool operator < (const PQEle& b) const
    {   return u < b.u; }
};

TEST(TestPQ, push_and_pop)
{
    PriorityQueue<PQEle>    pq;
    pq.resize(100);
    pq.push(new PQEle(10.));
    pq.push(new PQEle(20.));
    pq.push(new PQEle(15.));
    pq.push(new PQEle(5.));
    ASSERT_EQ(5., pq.peek()->u);
    pq.pop();
    ASSERT_EQ(10., pq.peek()->u);
    pq.pop();
    ASSERT_EQ(15., pq.peek()->u);
    pq.push(new PQEle(7.));
    ASSERT_EQ(7., pq.peek()->u);
    pq.push(new PQEle(9.));
    ASSERT_EQ(7., pq.peek()->u);
    pq.pop();
    ASSERT_EQ(9., pq.peek()->u);
}

TEST(TestPQ, move_up_down)
{
    PriorityQueue<PQEle>    pq;
    pq.resize(100);
    PQEle* eles[] = {
        new PQEle(10.),
        new PQEle(20.),
        new PQEle(15.),
        new PQEle(5.)};
    for(int i = 0;i < 4;++ i) pq.push(eles[i]);
    ASSERT_EQ(5., pq.peek()->u);
    eles[3]->u = 12;
    pq.move_down_node(eles[3]);
    ASSERT_EQ(10., pq.peek()->u);
    eles[2]->u = 9;
    pq.move_up_node(eles[2]);
    ASSERT_EQ(9., pq.peek()->u);
}

