#!/usr/bin/env python 

import numpy as np 
import matplotlib.pyplot as plt 


test1 = np.loadtxt('timing_fixeddim.txt', skiprows=1) 
test2 = np.loadtxt('timing_fixednnz.txt', skiprows=1) 

test3 = np.loadtxt('timing_fixednnz_noDense.txt', skiprows=1) 

test5 = np.loadtxt('timing_fixeddim_biCGStab.txt', skiprows=1) 
test4 = np.loadtxt('timing_fixednnz_biCGStab.txt', skiprows=1) 



plt.figure() 
plt.title('Fixed dimension = %u' %(test1[0,0])) 
plt.loglog(test1[:,1], test1[:,2],label='denseQR solve') 
plt.loglog(test1[:,1], test1[:,3],label='sparseQR solve') 
plt.loglog(test1[:,1], test1[:,4],label='decoupled solve') 
plt.loglog(test5[:,1], test5[:,2],label='sparseBiCGStab solve') 
plt.xlabel('number of non-zero entries')
plt.ylabel('time (sec)')
plt.legend(loc=2)


plt.figure() 
plt.title('Fixed non-zero elements = %u' %(test2[0,1])) 
plt.loglog(test2[:,0], test2[:,2],label='denseQR solve') 
plt.loglog(test2[:,0], test2[:,3],label='sparseQR solve') 
plt.loglog(test2[:,0], test2[:,4],label='decoupled solve') 
plt.loglog(test4[:,0], test4[:,2],label='sparseBiCGStab solve') 
plt.xlabel('dimension N')
plt.ylabel('time (sec)')
plt.legend(loc=2)

plt.figure() 
plt.title('Fixed non-zero elements = %u' %(test3[0,1])) 
plt.loglog(test3[:,0], test3[:,2],label='sparseQR solve') 
plt.loglog(test3[:,0], test3[:,3],label='decoupled solve') 
plt.loglog(test4[:,0], test4[:,2],label='sparseBiCGStab solve') 
plt.xlabel('dimension N')
plt.ylabel('time (sec)')
plt.legend(loc=2)

plt.show()
