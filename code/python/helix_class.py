import coil_class as cl
import math
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

class helix():
    def __init__(self,cldat,scldat):

        self.__coil_data = pd.read_csv(cldat)
        # establish indices
        self.__coil_data.set_index(list(self.__coil_data)[0],inplace=True)

        self.__subcoil_data = pd.read_csv(scldat)
        # establish indices
        sc_names=list(self.__subcoil_data)[:2]
        sc_index = list(zip(self.__subcoil_data[sc_names[0]],self.__subcoil_data[sc_names[1]]))
        sc_index = pd.MultiIndex.from_tuples(sc_index)
        self.__subcoil_data.drop(sc_names[0],axis=1,inplace=True)
        self.__subcoil_data.drop(sc_names[1],axis=1,inplace=True)
        self.__subcoil_data.index = sc_index
        self.__subcoil_data.index.names=sc_names

        self.__rot = dict()
        for index in self.__coil_data.index:
            self.__rot[index] = self.__rotation(index)

################################################################################

    def printcoil(self):
        print(self.__coil_data)

################################################################################

    def printsubcoil(self):
        print(self.__subcoil_data)

################################################################################

    def origin(self,coil_index):
        """
        Returns location of center of coil

        Input: pandas array index
        Returns: 1x3 numpy array
        """
        return np.array(self.__subcoil_data.loc[coil_index][2:5])

################################################################################

    def __rotation(self,coil_index):
        """
        Returns rotation matrix for given coil

        Input: pandas array index
        Returns: 3x3 numpy array
        """

        tilt = np.array(self.__coil_data.loc[coil_index][-2:])

        cp = math.cos(tilt[1])
        sp = math.sin(tilt[1])
        ct = math.cos(tilt[0])
        st = math.sin(tilt[0])

        rotz = np.array([[cp,sp,0.],[-sp,cp,0.],[0.,0.,1.]])
        roty = np.array([[ct,0,-st],[0.,1.,0.],[st,0.,ct]])

        return roty.dot(rotz)

################################################################################

    def __subcoil_field(self,coil_index,subcoil_index,coords):


        #translate global origin to coil origin
        coords -= self.__coil_data.loc[coil_index][2:5]

        #rotate coordinates so z-axis corresponds with coil axis
        coords = self.__rot[coil_index].dot(coords)

        rho_out = self.__subcoil_data.loc[coil_index,subcoil_index][1]
        rho_in = self.__subcoil_data.loc[coil_index,subcoil_index][0]
        width_rho = rho_out - rho_in
        nrho = self.__subcoil_data.loc[coil_index,subcoil_index][3]
        drho = width_rho/float(nrho)

        width_z = self.__coil_data.loc[coil_index][1]
        nz = self.__subcoil_data.loc[coil_index,subcoil_index][4]
        dz = width_z/float(nz)

        dturns = self.__subcoil_data.loc[coil_index,subcoil_index][2]*dz*drho/(width_z*width_rho) #turns in ideal coil
        di = dturns*self.__coil_data.loc[coil_index][0]
        B = np.zeros(3)
        #break subcoil into ideal coils
        for i in range(int(nrho)):
            idc_radius = rho_in + (i+0.5)*drho #ideal coil radius
            for j in range(int(nz)):
                idc_z = (j+0.5)*dz - width_z/2. #ideal coil z-offset

                idc = cl.coil(current = di, radius = idc_radius, center = (0.,0.,idc_z))
                B += idc.Bxyz(coords)

        return B

################################################################################

    def B(self, coords):

        if len(coords)!=3:
            print("Coordinates have length ",len(coords))
            raise Exception("Error: invalid number of coordinates")
        coords = np.array(coords)
        B = np.zeros(3)
        for subcoil in self.__subcoil_data.index:
            B += self.__subcoil_field(subcoil[0],subcoil[1],coords)
        return B

################################################################################

    # def PlotB(self,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax):
    #     fig = plt.figure()
    #     ax = fig.gca(projection='3d')
    #
    #     x,y,z=np.meshgrid(np.arange(xmin,xmax,(xmax-xmin)/Nx),\
    #                     np.arange(ymin,ymax,(ymax-ymin)/Ny),\
    #                     np.arange(zmin,zmax,(zmax-zmin)/Nz))
    #     map(self.B,list(zip(x,y,z)))
    #
    #     ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)
    #
    #     plt.show()
