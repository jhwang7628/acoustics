from pylab import *
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
 
def animate(nframe):
        plt.clf()
        print "Frame " + str(nframe)
        data = numpy.loadtxt("frames/frame" + str(nframe))
        # imshow(data, aspect='auto', cmap=get_cmap('seismic'),origin="lower")
        imshow(data, aspect='auto', cmap=get_cmap('seismic'), vmin=-0.1, vmax=0.1,origin="lower")
        # savefig("frames/frame" + str(i) + ".png",dpi=100,facecolor='gray')
 
def run():
        # Writer = animation.FFMpegWriter
        # writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        fig = plt.figure(figsize=(5,5))
        anim = animation.FuncAnimation(fig, animate, frames=200)
        anim.save('frames/wave.mp4', fps=30)
        return anim
 
if __name__ == '__main__':
        run()