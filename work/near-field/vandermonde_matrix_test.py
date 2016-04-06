#!/usr/bin/env python 

import numpy as np 


def ComputeRow(x): 
    return np.array([x[0]*x[1]*x[2], x[0]*x[1], x[0]*x[2], x[1]*x[2], x[0], x[1], x[2], 1])


def Transform(a,b,c,x0,y0,z0): 
    A = np.zeros((8,8))

    A[0,0] = a*b*c

    A[1,0] = a*b*z0
    A[1,1] = a*b

    A[2,0] = a*c*z0
    A[2,2] = a*c

    A[3,0] = c*b*x0
    A[3,3] = c*b

    A[4,0] = a*y0*z0
    A[4,1] = a*y0
    A[4,2] = a*z0
    A[4,4] = a

    A[5,0] = b*x0*z0
    A[5,1] = b*x0
    A[5,3] = b*z0
    A[5,5] = b
    
    A[6,0] = c*x0*y0
    A[6,2] = c*x0
    A[6,3] = c*y0
    A[6,6] = c

    A[7,0] = x0*y0*z0
    A[7,1] = x0*y0
    A[7,2] = x0*z0
    A[7,3] = y0*z0
    A[7,4] = x0
    A[7,5] = y0
    A[7,6] = z0
    A[7,7] = 1


    print 'transform'
    print A
    print np.linalg.cond(A)




low = np.array([0.085,0.085,0.085]) 
h = 0.01

A = np.zeros((8,8))

A[0,:] = ComputeRow(np.array([low[0]+0,low[1]+0,low[2]+0]))
A[1,:] = ComputeRow(np.array([low[0]+0,low[1]+0,low[2]+h]))
A[2,:] = ComputeRow(np.array([low[0]+0,low[1]+h,low[2]+0]))
A[3,:] = ComputeRow(np.array([low[0]+0,low[1]+h,low[2]+h]))
A[4,:] = ComputeRow(np.array([low[0]+h,low[1]+0,low[2]+0]))
A[5,:] = ComputeRow(np.array([low[0]+h,low[1]+0,low[2]+h]))
A[6,:] = ComputeRow(np.array([low[0]+h,low[1]+h,low[2]+0]))
A[7,:] = ComputeRow(np.array([low[0]+h,low[1]+h,low[2]+h]))

print A

print np.linalg.cond(A)


Transform(h,h,h,0,0,0)


