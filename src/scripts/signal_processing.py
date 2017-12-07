#/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
################################################################################
# Function ComputeAverageFFT
################################################################################
def ComputeAverageFFT(dt, y, N=1024, overlap_ratio=0.5): 
    # N = len(y)
    x = np.linspace(0.0, float(len(y))*dt, len(y))

    xf_accum = None
    yf_accum = np.zeros(N//2, dtype='complex128')
    magyf = np.zeros(N//2)

    ii = 0
    N_accum = 0
    while ii<len(y): 
        ii_m = ii+N
        if (ii_m >= len(y)-1):
            break
        N_windowed = ii_m - ii
        y_windowed = np.multiply(y[ii:ii_m], np.hanning(N))
        assert(len(y_windowed) == N_windowed)

        xf_accum  = np.linspace(0.0, 1.0/(2.0*dt), N/2)
        yf = np.fft.fft(y_windowed)[:N//2]
        yf_accum += yf
        magyf += np.abs(yf)

        ii += int(N*overlap_ratio)
        N_accum += 1

    magyf /= float(N_accum)
    yf_accum /= float(N_accum)
    return xf_accum, yf_accum, magyf

################################################################################
# Function Main
################################################################################
if __name__ == '__main__': 

    # testing for compute fft
    N = 600
    # sample spacing
    T = 1.0 / 800.0
    x = np.linspace(0.0, N*T, N)
    y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
    xf, yf, magyf = ComputeAverageFFT(T, y, N=512)
    plt.figure()
    plt.plot(xf, magyf)
    plt.show()
