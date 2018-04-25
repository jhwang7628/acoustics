#!/usr/bin/env python
import sys,math

def Compute_tmax(a, h, c=343.):
    # tmax = h / (math.sqrt(float(3.))*343.)
    # from Webb 2011, COMPUTING ROOM ACOUSTICS WITH CUDA - 3D FDTD SCHEMES WITH BOUNDARY LOSSES AND VISCOSITY
    return (-6.*a*c + math.sqrt((6.*a*c)**2 + 12.*(c*h)**2)) / (6.*c**2)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print '**Usage: %s <cell_size> [damping_alpha=0]' %(sys.argv[0])
        sys.exit(1)

    h = float(sys.argv[1])
    if len(sys.argv) == 3: a = float(sys.argv[2])
    else:                  a = 0.0
    c = 343.

    tmax = Compute_tmax(a, h, c)
    print 'Min sampling frequency = %u (max time step size = %.12f), if h = %f for damping parameter %f' %(int(1./tmax), tmax, h, a)

