import helix_class as hx
import time
import numpy as np
import math


def main():
    # Location of spreadsheets
    generaldata = '../../data/GeneralData.csv'
    coildata = '../../data/CoilData.csv'
    subcoildata = '../../data/SubcoilData.csv'
    # Make helix_class object
    helix_magnet = hx.helix(generaldata,coildata,subcoildata)

    # Check consistency with fortran routine
    ## run ./../fortran/helixtest to compare
    coord = [[0.1,0.1,0.05]]
    print('Magnetic field at',coord[0])
    print(helix_magnet.B(coord))
    print('Phi rotation with no theta rotation should not affect field due to axisymmetry')
    thetaphi = [0.,0.1]
    helix_magnet.set_coil_angle(0,thetaphi)
    print(helix_magnet.B(coord))
    print('Theta = pi rotation shoud effectively reverse current')
    thetaphi = [math.pi,0.]
    helix_magnet.set_coil_angle(0,thetaphi)
    helix_magnet.set_coil_angle(1,thetaphi)
    print(helix_magnet.B(coord))
    thetaphi = [0.,0.] # Reset angle
    
    
    # Example 1: how to get large number of field points at once
    ## Make some coordinates
    x = np.linspace(-0.1,0.1,num=20)
    y = np.linspace(-0.1,0.1,num=20)
    z = np.linspace(-0.1,0.1,num=20)
    coord = list()
    for cx in x:
        for cy in y:
            for cz in z:
                coord.append((cx,cy,cz))
    ## Calculate field
    B = helix_magnet.B(coord)
    ## Save to file
    filename = 'field.txt'
    with open(filename,"w") as fieldfile:
        index = 0
        for value in B:
            print("{}\t{}".format(coord[index],value),file=fieldfile)
            index+=1
        print('Wrote example 1 to',filename)
    


main()
