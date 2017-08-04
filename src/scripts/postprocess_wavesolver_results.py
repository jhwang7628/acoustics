#!/usr/bin/env python
import sys,os
import matplotlib.pyplot as plt
from wavesolver_results import *
from ConfigParser import SafeConfigParser
import scipy
from scipy import signal

if len(sys.argv) != 2: 
    print '**Usage: %s <config_ini_file>' %sys.argv[0]
    sys.exit()
out_config = sys.argv[1]
if not os.path.isfile(out_config): 
    print '**ERROR** Config file %s does not exist' %(out_config)
    sys.exit()

##
def Parse(parser, section, tag, val_type = 's'): 
    if val_type == 'i': 
        tmp = parser.getint(section, tag)
    elif val_type == 'f':
        tmp = parser.getfloat(section, tag)
    elif val_type == 'b':
        tmp = parser.getboolean(section, tag)
    elif val_type == 's':
        tmp = parser.get(section, tag)
    else: 
        raise TypeError
    print 'Parsing \'%s\' in \'%s\': ' %(tag, section), tmp
    return tmp
##
parser = SafeConfigParser()
parser.read(out_config)

data_dir = Parse(parser, 'general', 'data_dir', 's')
out_dir = Parse(parser, 'general', 'out_dir', 's')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
##
results = Wavesolver_Results()
results.Set_Folder(data_dir)
all_data = results.Read_All_Audio()
N_points = all_data.shape[1]
N_steps  = all_data.shape[0] 

if Parse(parser, 'general', 'plot', 'b'):
    print '\n------ PLOTTING ------'
    plt.figure()
    for ii in range(N_points): 
        plt.plot(all_data[:,ii])
    plt.show()

if Parse(parser, 'general', 'write_wav', 'b'):
    print '\n------ WRITING WAV FILES ------'
    sampfreq = Parse(parser, 'general', 'sampfreq', 'i')
    wavfreq  = Parse(parser, 'wav'    , 'wavfreq' , 'i')
    prefix   = Parse(parser, 'wav'    , 'prefix'  , 's')
    rateRatio = float(sampfreq) / float(wavfreq)
    for ii in range(N_points): 
        print 'point %u' %(ii)
        outputdata = signal.resample(all_data[:,ii], int(float(N_steps)/rateRatio))
        normalization = np.absolute(outputdata).max()
        scipy.io.wavfile.write('%s/%s_%u.wav' %(out_dir, prefix, ii), wavfreq, outputdata/normalization)
