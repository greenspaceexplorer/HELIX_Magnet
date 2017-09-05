import math
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import helix_mag as hlx

class helix():
    def __init__(self,gdat,cldat,scldat):

        # Read in general data spreadsheet
        self.__general_data = pd.read_csv(gdat)
        # Get general data column names
        g_names = list(self.__general_data)
        self.current = self.__general_data[g_names[1]][0]
        self.contraction = [self.__general_data[g_names[2]][0],self.__general_data[g_names[3]][0]]

        # Read in coil data spreadsheet
        self.__coil_data = pd.read_csv(cldat)

        # Get coil data column names
        c_names = list(self.__coil_data)

        # Make coil data lists to pass to fortran routines
        self.width = self.__coil_data[c_names[1]].tolist()
        self.clx = self.__coil_data[c_names[2]].tolist()
        self.cly = self.__coil_data[c_names[3]].tolist()
        self.clz = self.__coil_data[c_names[4]].tolist()
        self.clposition = list(zip(self.clx,self.cly,self.clz))
        self.theta = self.__coil_data[c_names[5]].tolist()
        self.phi = self.__coil_data[c_names[6]].tolist()

        # Get coil names for subcoil list creation
        coil_index = self.__coil_data[c_names[0]]

        # Read in subcoil data spreadsheet
        self.__subcoil_data = pd.read_csv(scldat)
        # Get subcoil data column names
        sc_names = list(self.__subcoil_data)

        # Make subcoil lists
        self.ri = list()
        self.ro = list()
        self.turns = list()
        self.nrho = list()
        self.nz = list()

        # Populate subcoil lists with a list for each primary coil
        for i in range(len(coil_index)):
            self.ri.append(list())
            self.ro.append(list())
            self.turns.append(list())
            self.nrho.append(list())
            self.nz.append(list())

        # Extract data from spreadsheet to subcoil lists
        for j in range(len(self.__subcoil_data[sc_names[0]])):
            index = 0
            for coil_name in coil_index:
                if self.__subcoil_data[sc_names[0]][j] == coil_name:
                    self.ri[index].append(self.__subcoil_data[sc_names[2]][j])
                    self.ro[index].append(self.__subcoil_data[sc_names[3]][j])
                    self.turns[index].append(self.__subcoil_data[sc_names[4]][j])
                    self.nrho[index].append(self.__subcoil_data[sc_names[5]][j])
                    self.nz[index].append(self.__subcoil_data[sc_names[6]][j])
                index+=1

        # Send coil data to fortran code        
        hlx.helix_data(self.current,self.contraction,self.width,self.turns,self.ri,self.ro,self.nrho,self.nz,self.clposition,self.theta,self.phi)
#    def __check_data(self):
#    def __set_data(self):
        
    def print_coil(self):
        print(self.__coil_data)
    def print_subcoil(self):
        print(self.__subcoil_data)
    def ideal(self,i,r,x):
        return hlx.ideal(i,r,x)
    def B(self,x):
        return hlx.helix_magnet(x)
        


 
