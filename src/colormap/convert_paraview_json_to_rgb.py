#!/usr/bin/env python
import numpy as np
import sys

if len(sys.argv) != 2:
    print '**Usage: %s <paraview_colormap_json>' %(sys.argv[0])
    sys.exit(1)

print 'Reading json file exported from paraview, assuming ranging from [-1,1] ...'
with open(sys.argv[1], 'r') as stream:
    lines = stream.readlines()
    idxs = []
    for ii in range(len(lines)):
        if lines[ii].find('RGBPoints') != -1:
            idxs.append(ii + 1)
        if lines[ii].find(']') != -1:
            idxs.append(ii - 1)
    N = (idxs[1]-idxs[0]+1) / 4
    colormap = np.zeros((N,4))
    for ii in range(N):
        for jj in range(4):
            line = lines[idxs[0] + ii*4 + jj]
            colormap[ii,jj] = float(line.split(',')[0])

print 'Interpolating map'
N_out = 255
dt = 2./float(N_out+1)
t = -1.0
colormap_out = np.zeros((N_out+1, 3))
for ii in range(N_out+1):
    idx0 = 0
    for jj in range(colormap.shape[0]):
        if colormap[jj,0] > t:
            idx0 = max(jj - 1, 0)
            break
    idx1 = min(idx0+1, N-1)
    t0 = colormap[idx0,0]
    t1 = colormap[idx1,0]
    alpha = (t-t0)/(t1-t0)

    for jj in range(3):
        colormap_out[ii,jj] = colormap[idx0,jj+1]*(1.-alpha) \
                            + colormap[idx1,jj+1]*(   alpha)
    t += dt

print 'N = ', N_out
for ii in range(colormap_out.shape[0]):
    print '  {%f,%f,%f},' %(colormap_out[ii,0], colormap_out[ii,1], colormap_out[ii,2])
