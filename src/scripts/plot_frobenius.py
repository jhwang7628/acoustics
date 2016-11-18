#!/usr/bin/env python
import numpy,sys
import matplotlib.pyplot as plt
if len(sys.argv) < 2:
    print '**Usage: %s <log_file> <log_file2> ...' %(sys.argv[0])
    sys.exit()
N_files = len(sys.argv) -1
plt.figure(figsize=[14,8])
for filename in sys.argv[1:]:
    lines = open(filename).readlines()
    frobenius = []
    for l in lines: 
        tokens = l.split()
        if len(tokens) == 0:
            continue
        if tokens[0] == 'total' and tokens[1] == 'energy': 
            frobenius.append(float(tokens[-1]))
    label=filename.split('/')[0].split('_')[-1]
    plt.plot(frobenius[0:], label=label)
plt.legend()
plt.ylim([-1E9, 1E9])
plt.title('Mug impulse response test: Total Energy')
plt.xlabel('Frame #')
plt.ylabel('Unnormalized total energy')
plt.savefig('mug_ir_energy.png')
# plt.show()
