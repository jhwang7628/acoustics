#!/usr/bin/env python
import struct,sys,math
from mode_data_structures import Modes
################################################################################
################################################################################
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print '**Usage: %s <modes_file> <density>' %(sys.argv[0])
        sys.exit(1)

    modes_file = sys.argv[1]
    density = float(sys.argv[2])
    Modes.Print_Eigen_Values_From_File(modes_file, density)

