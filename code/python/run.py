import coil_class as cl
import helix_class as hx
import time

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

def main():
    i = 1.e6
    center = (.15,.14,.13)
    tilt = (.2,.3)
    rad = .35

    c1 = cl.coil(i,rad,center,tilt)

    print(c1)

    print(c1.Bxyz((.15,.14,.13)))

    coil_data = "../data/CoilData.csv"
    subcoil_data = "../data/SubcoilData.csv"

    hlx = hx.helix(coil_data,subcoil_data)
    hlx.printcoil()
    hlx.printsubcoil()

    t = time.time()
    B = hlx.B((0.1,0.1,0.1))
    t = time.time() - t

    print("[Bx By Bz] = ", B)

    print("Time to compute: {} seconds".format(t))

    # hlx.PlotB(10,10,10,-.4,.4,-.4,.4,-.4,.4)

    # Plot B field

    xmin = -.4
    xmax = 0.4
    ymin = xmin
    ymax = xmax
    zmin = xmin
    zmax = xmax
    Nx = Ny = Nz = 5

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x,y,z = np.meshgrid(np.arange(xmin,xmax,(xmax-xmin)/Nx),\
                        np.arange(ymin,ymax,(ymax-ymin)/Ny),\
                        np.arange(zmin,zmax,(zmax-zmin)/Nz))

    x=x.reshape(1,len(x)**3)[0]
    y=y.reshape(1,len(y)**3)[0]
    z=z.reshape(1,len(z)**3)[0]
    Bb = list(zip(*list(map(hlx.B,zip(x,y,z)))))

    ax.quiver(x, y, z, Bb[0],Bb[1],Bb[2], length=0.1, normalize=True)
    plt.show()






main()
