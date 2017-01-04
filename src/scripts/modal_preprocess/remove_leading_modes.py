#!/usr/bin/env python
import sys,os
from mode_data_structures import Modes
################################################################################
################################################################################
if __name__ == "__main__":
    if len(sys.argv) != 3: 
        print '**Usage: %s <in_mode_file> <out_mode_file>' %(sys.argv[0])
        sys.exit(1)

    mode_file_i = sys.argv[1]
    mode_file_o = sys.argv[2]

    modes = Modes.Read_Modes(mode_file_i)
    modes.Remove_Leading_Modes(1)
    modes.Write_Modes(mode_file_o)
