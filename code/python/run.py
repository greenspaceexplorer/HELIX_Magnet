import helix_class as hx
import time
import numpy as np
import math


def main():
    # Initialize data
    current = 93.0
    contraction = [0.99674,0.99585]
    width = [0.075818,0.075818]
    position_a = [0.,0.,-0.359747]
    position_b = [0.,0.,0.359747]
    position = [position_a,position_b]
    theta = [0.,0.]
    phi = [0.,0.]
    # Subcoil data
    sc_a1 = [.203517,.215011,1995.9,8,52]
    sc_a2 = [.215011,.234632,5150.,8,32]
    sc_a3 = [.234632,.255224,5489.5,8,32]
    sc_b1 = [.203492,.214731,1976.,8,52]
    sc_b2 = [.214732,.234709,5110.,8,32]
    sc_b3 = [.234709,.256222,5486.7,8,32]
    sclist_a = list(zip(sc_a1,sc_a2,sc_a3))
    sclist_b = list(zip(sc_b1,sc_b2,sc_b3))
    turns = [sclist_a[2],sclist_b[2]]
    inner_radius = [sclist_a[0],sclist_b[0]]
    outer_radius = [sclist_a[1],sclist_b[1]]
    nrho = [sclist_a[3],sclist_b[3]]
    nz = [sclist_a[4],sclist_b[4]]

    helix_magnet = hx.helix(current,contraction,width,position,theta,phi,turns,inner_radius,outer_radius,nrho,nz)


    # Example 1: change rotation and position of coils
    coord = [[0.1,0.1,0.05]]
    print('---Magnetic field at',coord[0])
    print(helix_magnet.B(coord))

    print('---Phi rotation with no theta rotation should not affect field due to axisymmetry')
    thetaphi = [0.,0.1]
    helix_magnet.set_coil_angle(0,thetaphi)
    print(helix_magnet.B(coord))
    
    print('---Theta = pi rotation shoud effectively reverse current')
    thetaphi = [math.pi,0.]
    helix_magnet.set_coil_angle(0,thetaphi)
    helix_magnet.set_coil_angle(1,thetaphi)
    print(helix_magnet.B(coord))
    
    print("---Reset coil rotation")
    thetaphi = [0.,0.] # Reset angle
    helix_magnet.set_coil_angle(0,thetaphi)
    helix_magnet.set_coil_angle(1,thetaphi)
    print(helix_magnet.B(coord))
    
    print("---Swapping coil positions should leave the field unchanged")
    helix_magnet.set_coil_position(0,position_b)
    helix_magnet.set_coil_position(1,position_a)
    print(helix_magnet.B(coord))

    print("---Moving the coils 50% closer should increase the field")
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    close_a = 0.5*position_a
    close_b = 0.5*position_b
    helix_magnet.set_coil_position(0,close_a)
    helix_magnet.set_coil_position(1,close_b)
    print(helix_magnet.B(coord))

    print("---Reset coil position")
    helix_magnet.set_coil_position(0,position_a)
    helix_magnet.set_coil_position(1,position_b)
    print(helix_magnet.B(coord))

    
    # Example 2: how to get large number of field points at once
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
        print('Wrote example 2 to',filename)
    


main()
