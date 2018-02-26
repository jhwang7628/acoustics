#!/usr/bin/env python
import numpy as np
import os,sys

################################################################################
# Function FilterImpulses
#   Filter out impulses whose impact velocity is above given threshold
################################################################################
def FilterImpulses(innfile, outfile, thres):
    # read and parse input
    records = []
    colpairstack = []
    with open(innfile) as stream:
        lines = stream.readlines()
        for i, line in enumerate(lines):
            if i%1000 == 0:
                print i
            tokens = line.split()

            # If its collision pair, first push it in stack and check when
            # there's a pair. If its collision constraint, then directly append
            # the qualiified ones.
            if tokens[8] == 'C':
                if abs(float(tokens[3])) > thres:
                    records.append(line)
            elif tokens[8] == 'P':
                colpairstack.append(line)
                if len(colpairstack) == 2:
                    qualified = False
                    if abs(float(colpairstack[0].split()[3])) > thres or \
                       abs(float(colpairstack[1].split()[3])) > thres:
                           qualified = True
                    if qualified:
                        records.append(colpairstack.pop())
                        records.append(colpairstack.pop())
                    else:
                        colpairstack = []

    # write results
    with open(outfile, 'w') as stream:
        for r in records:
            stream.write(r)

    print records, len(records)

################################################################################
################################################################################
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print '**Usage: %s <in_file> <out_file> <vel_threshold>' %(sys.argv[0])
        sys.exit(1)
    FilterImpulses(sys.argv[1], sys.argv[2], float(sys.argv[3]))
