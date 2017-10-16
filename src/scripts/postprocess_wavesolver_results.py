#!/usr/bin/env python
import sys,os
import matplotlib.pyplot as plt
from wavesolver_results import *
from ConfigParser import SafeConfigParser,NoOptionError
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
def ReqParse(parser, section, tag, val_type = 's'):
    if val_type == 'i':   # int
        tmp = parser.getint(section, tag)
    elif val_type == 'f': # float
        tmp = parser.getfloat(section, tag)
    elif val_type == 'b': # boolean
        tmp = parser.getboolean(section, tag)
    elif val_type == 's': # string
        tmp = parser.get(section, tag)
    elif val_type == 'li': # list of int (comma separated)
        tmp = map(int, parser.get(section, tag).split(','))
    else:
        raise TypeError
    print 'Parsing \'%s\' in \'%s\': ' %(tag, section), tmp
    return tmp

def OptParse(parser, section, tag, val_type = 's', val_default = None):
    try: val = ReqParse(parser, section, tag, val_type)
    except NoOptionError: val = val_default
    return val
##
parser = SafeConfigParser()
parser.read(out_config)

data_dir = ReqParse(parser, 'general', 'data_dir', 's')
out_dir = ReqParse(parser, 'general', 'out_dir', 's')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
##
results = Wavesolver_Results()
results.Set_Folder(data_dir)
all_data = results.Read_All_Audio()
N_points = all_data.shape[1]
N_steps  = all_data.shape[0]
sampfreq = ReqParse(parser, 'general', 'sampfreq', 'i')

padTime = OptParse(parser, 'general', 'pad_time', 'f', 0)
if padTime > 0:
    pad_amount = int(round(padTime * sampfreq))
    all_data = np.pad(all_data, [(pad_amount, 0), (0,0)], 'constant', constant_values=0)
    N_steps  = all_data.shape[0]

if ReqParse(parser, 'general', 'plot', 'b'):
    print '\n------ PLOTTING ------'
    xaxis_frame = ReqParse(parser, 'plot', 'xaxis_frame', 'b')
    xmax        = OptParse(parser, 'plot', 'xmax'       , 'f')
    plot_points = OptParse(parser, 'plot', 'plot_points', 'li')
    plt.figure()
    if (plot_points is None): ii_range = range(N_points)
    else                    : ii_range = plot_points
    for ii in ii_range:
        if xaxis_frame: t = np.array(range(len(all_data[:,ii])))
        else:           t = np.array(range(len(all_data[:,ii])))/float(sampfreq)
        y = all_data[:,ii]
        if xmax is not None:
            cut = -1
            for ii in range(len(t)):
                if t[ii] > xmax:
                    cut = ii
                    break
            if cut >= 0:
                t = t[:cut]
                y = y[:cut]
        plt.plot(t, y)
    plt.show()

if ReqParse(parser, 'general', 'write_wav', 'b'):
    print '\n------ WRITING WAV FILES ------'
    wavfreq  = ReqParse(parser, 'wav'    , 'wavfreq' , 'i')
    prefix   = ReqParse(parser, 'wav'    , 'prefix'  , 's')
    rateRatio = float(sampfreq) / float(wavfreq)
    for ii in range(N_points):
        print 'point %u' %(ii)
        outputdata = signal.resample(all_data[:,ii], int(float(N_steps)/rateRatio))
        normalization = np.absolute(outputdata).max()
        scipy.io.wavfile.write('%s/%s_%u.wav' %(out_dir, prefix, ii), wavfreq, outputdata/normalization)
