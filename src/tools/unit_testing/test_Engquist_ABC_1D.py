#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import animation

float_formatter = lambda x: "% 2.2f" % x
np.set_printoptions(linewidth=999, formatter={'float_kind':float_formatter})

def Gaussian(x, x0, h): 
    return np.exp(-(x-x0)**2/2./h**2)


N = 200
h = 1./float(N)
t = h
l = t/h
l2 = l**2
c = 1.
Nstep = 1200
btype = 'Dirichlet'
# btype = 'Engquist Absorbing'

pim1 = np.zeros(N)
pip1 = np.zeros(N)
pi00 = np.zeros(N)

# initialization: left-traveling wave
for i in range(N): 
    x = h*float(i)
    pim1[i] = Gaussian(x, float(N)*3./4.*h+c*t, 0.025)
    pi00[i] = Gaussian(x, float(N)*3./4.*h    , 0.025)

# updates
colors = [cm.jet(x) for x in np.linspace(0, 1, Nstep)]
fig = plt.figure()
ax = plt.axes(xlim=(0,200), ylim=(-1,1))
ax.set_title('%s (x=0), Neumann (x=200)' %(btype))
line, = ax.plot([],[],lw=2)

def init(): 
    line.set_data([],[])
    return line, 

def animate(step): 
    global pim1,pip1,pi00,btype,colors
    # reflecting boundary
    pip1[N-1] = 2.*pi00[N-1] - pim1[N-1] + l2*(pi00[N-2] - pi00[N-1])

    # absorbing boundary
    if btype=='Dirichlet':
        pip1[1] = 2.*pi00[1] - pim1[1] + l2*(pi00[2] - 2.*pi00[1])
    elif btype=='Engquist Absorbing': 
        pip1[1] = l2/(1.+l)*((2./l2-2.)*pi00[1] + (l-1.)/l2*pim1[1] + 2.*pi00[2])

    for i in range(2,N-1): 
        pip1[i] = 2.*pi00[i] - pim1[i] + l2*(pi00[i+1] + pi00[i-1] - 2.*pi00[i])

    # update arrays
    pim1 = pi00.copy()
    pi00 = pip1.copy()
    print 'step %u' %step, pi00[:10], pi00[-10:]

    line.set_data(range(len(pi00)), pi00)
    line.set_color(colors[step])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Nstep, interval=20, blit=True)
anim.save('%s.mp4' %(btype), fps=60, extra_args=['-vcodec', 'libx264'])
plt.show()
