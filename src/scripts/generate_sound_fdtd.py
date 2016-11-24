#!/usr/bin/env python 
import binary_io as IO
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm 
import argparse,sys,glob,scipy,os
################################################################################
## Global variables
################################################################################
cases = []

################################################################################
################################################################################
def ParseArguments(argv):
    parser = argparse.ArgumentParser(description='Generate Sound FDTD')
    parser.add_argument('--plot', action='store_true', help='Plot the results')
    parser.add_argument('--case', type=ParseCase, nargs='*', help='Single run case. To use parser, input something like DIR=dir, START_TIME=start_time')
    parser.add_argument('--output_dir', default='output', help='Output directory')
    parser.add_argument('--file_prefix', default='test', help='File prefix')
    parser.add_argument('--sampling_rate', required=True, help='Sampling frequency of the fdtd simulation')
    parser.add_argument('--output_wav_rate', type=int, default=44100, help='Sampling frequency of the output wav files')
    parser.add_argument('--point', help='Listening point')
    args = parser.parse_args()
    return args

################################################################################
################################################################################
def ParseCase(s):
    try: 
        case = map(str, s.split(','))
        directory = case[0].split('=')[-1]
        start_time = float(case[1].split('=')[-1])
        cases.append(Case(directory, start_time))
    except:
        raise argparse.ArgumentTypeError('**ERROR** Parser error')

################################################################################
################################################################################
def PadZero(data, N_front_pad, N_end_pad):
    N_rows = len(data)
    newData = np.zeros(N_front_pad+N_rows+N_end_pad)
    newData[N_front_pad:N_front_pad+N_rows] = data

    return newData

################################################################################
################################################################################
class Case: 
    def __init__(self, d, t):
        self.directory = d
        self.start_time = t

    def __str__(self):
        return 'Case: directory=%s; start_time=%f' %(self.directory, self.start_time)

    def ReadListenedData(self, files, args):
        self.listeningPositions = IO.readMatrixXdBinary(files.listening_position_dat) 
        self.N_points = self.listeningPositions.shape[0]
        self.N_steps  = files.N_steps
        self.stepRate = args.sampling_rate
        self.wavRate  = args.output_wav_rate
        self.rateRatio = float(self.stepRate) / float(self.wavRate)

        # get listened data for each step
        self.listenedData = np.zeros((self.N_points, self.N_steps))
        step = 0
        for f in files.filenames: 
            self.listenedData[:, step] = IO.readMatrixXdBinary(f, verbose=0).flatten()
            step += 1

################################################################################
################################################################################
class SimFiles:
    def __init__(self, case, prefix): 
        pathprefix = os.path.join(case.directory, prefix)
        f_listening_position = '%s_listening_position.dat' %(pathprefix)
        f_data_prefix = '%s_data_listening_' %(pathprefix)
        self.listening_position_dat = f_listening_position
        self.data_prefix = f_data_prefix
        self.filenames = sorted(glob.glob(f_data_prefix+'[0-9]*'))
        self.N_steps = len(self.filenames)


################################################################################
################################################################################
def Plot():
    print 'Plot'

if __name__ == '__main__': 
    ## initialize pyplot
    plt.rcParams['agg.path.chunksize'] = 20000
    ##
    args = ParseArguments(sys.argv)
    for case in cases: 
        print '=== ', case, ' ==='
        # parse, initialization
        files = SimFiles(case, args.file_prefix)
        case.ReadListenedData(files, args)

    if (args.plot):
        Plot()
