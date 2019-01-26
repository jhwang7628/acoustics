#!/usr/bin/env python
import struct,sys,math
import numpy as np
from mode_data_structures import Modes
################################################################################
################################################################################
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print '**Usage: %s <modes_file> <density> <out_file>' %(sys.argv[0])
        sys.exit(1)

    modes_file = sys.argv[1]
    density = float(sys.argv[2])
    modes = Modes.Read_Modes(modes_file)
    print modes.num_modes
    freqs = modes.Get_All_Frequencies(density)

    np.savetxt(sys.argv[3], freqs)
