import helix_class as hx
import time
import numpy as np


def main():
    # Location of spreadsheets
    generaldata = '/Users/nwgreen/Google Drive/git-repo/HELIX_Magnet/data/GeneralData.csv'
    coildata = '/Users/nwgreen/Google Drive/git-repo/HELIX_Magnet/data/CoilData.csv'
    subcoildata = '/Users/nwgreen/Google Drive/git-repo/HELIX_Magnet/data/SubcoilData.csv'
    # Make helix_class object
    helix_magnet = hx.helix(generaldata,coildata,subcoildata)

    # Make some coordinates
    x = np.linspace(-0.1,0.1,num=20)
    y = np.linspace(-0.1,0.1,num=20)
    z = np.linspace(-0.1,0.1,num=20)
    coord = list()
    for cx in x:
        for cy in y:
            for cz in z:
                coord.append((cx,cy,cz))
    B = helix_magnet.B(coord)
    for value in B:
        print(value)
    


main()
