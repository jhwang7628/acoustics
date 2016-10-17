#!/usr/bin/env python 
import sys

def UpdateVelocity(p_0, p_1, dt=1./352800., dx=0.005, rho=1.184): 
    return -(p_1 -p_0)*dt/dx/rho


if len(sys.argv) != 5: 
    print '**Usage %s <mode v=0 p=1> <p_0/v_0> <p_1/v_1> <old_v/old_p>' %(sys.argv[0])
    sys.exit()

mode = int(sys.argv[1])
value_0 = float(sys.argv[2])
value_1 = float(sys.argv[3])
old_value = float(sys.argv[4])
if mode == 0: 
    dv = UpdateVelocity(value_0, value_1)
    print 'v_o = ', old_value
    print 'dv  = ', dv
    print 'v_n = ', old_value + dv
