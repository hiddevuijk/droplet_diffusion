import numpy as np
import matplotlib.pyplot as plt
from sys import exit


c = np.loadtxt("c.dat").T
drops = np.loadtxt("drops.dat")



fig, axes = plt.subplots(2)

ax = axes[0]
ax.imshow(c, interpolation='none')
#ax.colorbar()
ax.set_aspect( 1 )




ax = axes[1]

ax.set_xlim([0,50])
ax.set_ylim([0,25])
ax.set_aspect( 1 )

for i in range(drops.shape[0]):
    x = drops[i][0]
    y = drops[i][1]
    R = drops[i][2]

    circle = plt.Circle(( x, y), R) 
    circle.set_alpha(0.5) 
    ax.add_artist(circle)

plt.show()

