#!/usr/bin/env python
import sys
import numpy as np
import scipy.io.wavfile
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

class Normalize(colors.Normalize):
    """
    A class which, when called, can normalize data into
    the ``[0.0, 1.0]`` interval.

    """

    def __call__(self, value, clip=None):
        """
        Normalize *value* data in the ``[vmin, vmax]`` interval into
        the ``[0.0, 1.0]`` interval and return it.  *clip* defaults
        to *self.clip* (which defaults to *False*).  If not already
        initialized, *vmin* and *vmax* are initialized using
        *autoscale_None(value)*.
        Clip is in dB
        """
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        # Convert at least to float, without losing precision.
        (vmin,), _ = self.process_value(self.vmin)
        (vmax,), _ = self.process_value(self.vmax)
        if vmin == vmax:
            result.fill(0)   # Or should it be all masked?  Or 0.5?
        elif vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        else:
            #if clip:
            #    mask = np.ma.getmask(result)
            #    result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
            #                         mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data
            vmin = max(vmin, vmax - self.clip)
            resdat -= vmin
            resdat /= (vmax - vmin)
            #print(np.amax(resdat))
            #print(np.amin(resdat))
            #print(self.clip)
            #print(np.amax(resdat))
            #print(np.amin(resdat))

            result = np.ma.array(resdat, mask=result.mask)
        # Agg cannot handle float128.  We actually only need 32-bit of
        # precision, but on Windows, `np.dtype(np.longdouble) == np.float64`,
        # so casting to float32 would lose precision on float64s as well.
        if result.dtype == np.longdouble:
            result = result.astype(np.float64)
        if is_scalar:
            result = result[0]
        return result

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print '**Usage: %s <file_1> [file_2] ...' %(sys.argv[0])
        sys.exit()

    files = sys.argv[1:]
    print files

    Range = 40.0 #dB
    # Range = 80.0 #dB

    #cma = cm.get_cmap('gnuplot')
    # cma = cm.get_cmap('afmhot')
    # cma = cm.get_cmap('Reds')
    cma = cm.get_cmap('hot')
    # cma = cm.get_cmap('Greys')
    # cma = cm.get_cmap('OrRd')
    # cma = cm.get_cmap('BuGn')
    nfft = 1024
    # nfft = 2048
    # nover = round(0.5 * nfft)
    nover = round(0.9 * nfft)
    # nover = round(0.5 * nfft)

    plt.style.use('dark_background')

    plt.figure(figsize=[24,4])
    # plt.figure(figsize=[12,2*len(files)])
    # plt.figure(figsize=[24,6],dpi=300)
    for idx, a in enumerate(files):

        N_a, data_a = scipy.io.wavfile.read(a)

        max_a = np.max(np.abs(data_a))
        data_a = data_a.astype(float)/float(max_a)

        dt_a = 1./N_a
        t_a = np.linspace(0., len(data_a)*dt_a, len(data_a))

        myNorm = Normalize(None, None, Range)

        # plt.subplot(len(files),1,idx+1)
        plt.subplot(1,len(files),idx+1)
        plt.specgram(data_a, NFFT=nfft, noverlap=nover, Fs=N_a, cmap=cma, norm=myNorm)
        plt.xlim([0,max(t_a)])
        plt.xlim([0,2.25])
        # plt.ylim([0,14000])
        plt.ylim([0,1400])
        # plt.ylabel('Frequency (kHz)')
        # plt.plot(t_a, data_a)
        # plt.ylim([-1,1])

    plt.savefig('spectrogram_comparison.pdf',dpi=300)
    plt.savefig('spectrogram_comparison.jpg',dpi=300)
    plt.show()
    sys.exit()
