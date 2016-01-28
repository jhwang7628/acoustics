#!/usr/bin/env python

import sys
import numpy as np 
import matplotlib.pyplot as plt 

def feval(x, y, t, sigma_x, sigma_t, c): 
    norm_y   = np.linalg.norm(y) 
    norm_x_y = np.linalg.norm(x-y)
    if norm_x_y < 1E-12:
        return 0.
    return 1./norm_x_y * np.exp( -(norm_y/sigma_x)**2/2. -((t-norm_y/c)/sigma_t)**2/2. ) 

def main(): 

    # listening point
    X = np.loadtxt( 'head_listening.txt' ) 

    minBound = -1.701
    maxBound =  1.701
    N = 100 
    dx = (maxBound-minBound)/N
    dx3 = dx**3

    minBoundPlusdxOver2 = minBound+dx/2.

    t = np.linspace(0,0.0365,X.shape[1])

    value_x = np.zeros( (X.shape[0],t.size) ) 

    for pp in range(0,X.shape[0],5): 
        print pp
        x = X[pp,:]
        for ii in range( N ): 
            for jj in range( N ): 
                for kk in range( N ): 
                    y = np.ndarray((3,), buffer=np.array([ii*dx,jj*dx,kk*dx])) 
                    y+= minBoundPlusdxOver2
                    value_x[pp,:] += feval(x,y,t,0.5,1./10000.,343.) * dx3

    value_x /= 4./np.pi

    np.savetxt('value_x.txt', value_x)


if __name__ == "__main__": 
    main()
