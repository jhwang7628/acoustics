#!/usr/bin/env python 
import my_io as io 
import scipy
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import sys

################################################################################
################################################################################
def FillQ(Q, q, ind): 
    Q[:, ind] = q

################################################################################
################################################################################
def LSQ(A, b): 
    Q, R = np.linalg.qr(A, 'full')
    c = scipy.linalg.solve_triangular(R, Q.transpose().dot(b)) 
    return c

################################################################################
################################################################################
def Push_back(A, q): 
    A_copy = A.copy()
    N_cols = A.shape[1]
    for c in range(N_cols-1): 
        A[:, c] = A_copy[:, c+1]
    A[:, N_cols-1] = q

################################################################################
################################################################################
def ReadTimeSignal(sampleID, filename): 
    data = io.readMatrixXdBinary(filename)
    # sampleID = 49757
    return data[:, sampleID]

################################################################################
################################################################################
if __name__ == '__main__':

    q_file = 'test_q.txt'
    t_file = 'test_displacement.txt'
    sampleID = 40000
    qSeries = io.readMatrixXdBinary(q_file)
    tSeries = ReadTimeSignal(sampleID, t_file)
    print qSeries
    
    N_steps = qSeries.shape[0]
    N_modes = qSeries.shape[1]
    print N_steps, N_modes
    relTol = 1E-3
    rank = 3
    
    Q = np.zeros((N_modes, rank))
    
    # initialize
    for c_idx in range(rank): 
        FillQ(Q, qSeries[c_idx, :], c_idx)

    # solve
    tsStart = rank
    tsSteps = N_steps - rank
    tsStop = tsStart + tsSteps
    RESIDUAL = np.zeros((N_modes, tsSteps))
    ITERATION = np.zeros(tsSteps)
    ERROR = np.zeros(tsSteps)
    count = 0
    for ts_idx in range(rank + tsStart, tsStop): 

        q = qSeries[ts_idx, :]
        eps = relTol * np.linalg.norm(q)
        # get LSQ estimate 
        # print 'rank = ', np.linalg.matrix_rank(Q) 
        if np.linalg.matrix_rank(Q) == 0:
            print 'matrix is of 0-rank, skip this iterate'
            c = np.zeros(rank)
        else: 
            c, r, rank, s = np.linalg.lstsq(Q, q)
            # c = LSQ(Q, qSeries[ts_idx, :])
        residual = q - Q.dot(c) 

        dq = np.zeros(N_modes)
        N_iteration = 0
        while np.linalg.norm(residual) > eps: 
            k = np.argmax(np.absolute(residual))
            dq[k] += residual[k]
            residual[k] = 0
            N_iteration += 1

        # print c, residual
        ITERATION[count] = N_iteration
        ERROR[count] = np.linalg.norm(residual)
        RESIDUAL[:, count] = residual
        print 'step = %u, N_iterations = %u, error = %f' %(ts_idx, N_iteration, np.linalg.norm(residual))
        print ' c = ', c

        Push_back(Q, Q.dot(c) + dq)
        count += 1

    # fig = plt.figure(figsize=(16, 12))
    plt.figure()
    gs = gridspec.GridSpec(5, 1, width_ratios=[3, 1])
    plt.subplot(gs[0])
    # plt.pcolor(RESIDUAL, cmap='RdBu', vmin=-1E-9, vmax=1E-9)
    plt.pcolor(RESIDUAL, cmap='RdBu', vmin=-eps/3.0, vmax=eps/3.0)
    plt.ylabel('Modes')
    plt.title('Signed residual q - Qc')
    # plt.colorbar()
    plt.subplot(gs[1])
    plt.plot(range(tsSteps), ITERATION, 'b', label='iteration count (NNZ)')
    plt.plot(range(tsSteps), np.ones(tsSteps)*N_modes, 'k--', label='Total mode count')
    plt.ylabel('Modes')
    plt.legend(loc=3)
    plt.subplot(gs[2])
    plt.plot(range(tsSteps), ERROR, 'r', label='error')
    plt.ylabel('||q - Qc - dq||')
    plt.legend(loc=3)
    # plt.ylim([-20, 20])
    plt.subplot(gs[3])
    plt.plot(range(tsSteps),tSeries[tsStart:tsStop], label='Vertex %u' %(sampleID))
    plt.legend(loc=3)
    plt.xlabel('Time')
    plt.ylabel('Sample vertex displacement')
    plt.show()
