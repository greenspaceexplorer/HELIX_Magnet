import math
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import helix_mag as hlx # This is the fortran library

class helix():
    def __init__(self,gdat,cldat,scldat):
        """
        Initializes helix class object.

        Input: 
         gdat: General data spreadsheet
         cldat: Coil data spreadsheet
         scldat: Subcoil data spreadsheet
         -Note: see ../../data for spreadsheet format
        """

        # Filenames
        self.gdat = gdat
        self.cldat = cldat
        self.scldat = scldat

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

    def print_general(self):
        print(self.__general_data)

    def print_coil(self):
        print(self.__coil_data)

    def print_subcoil(self):
        print(self.__subcoil_data)

    def set_coil_position(self,coil,position):
        """
        Changes position of coil.
        
        Input:
         coil: integer, coil index in spreadsheet
         position: list of doubles, (x,y,z) origin of coil
        """
        self.clx[coil] = position[0]
        self.cly[coil] = position[1]
        self.clz[coil] = position[2]
        self.clposition[coil] = position
        hlx.helix_data(self.current,self.contraction,self.width,self.turns,self.ri,self.ro,self.nrho,self.nz,self.clposition,self.theta,self.phi)
    def set_coil_angle(self,coil,angle):
        """
        Changes angle of coil.

        Input:
         coil: integer, coil index in spreadsheet
         angle: list of doubles, (theta,phi) angle of coil
        """
        self.theta[coil] = angle[0]
        self.phi[coil] = angle[1]
        hlx.helix_data(self.current,self.contraction,self.width,self.turns,self.ri,self.ro,self.nrho,self.nz,self.clposition,self.theta,self.phi)
    def ideal(self,i,r,x):
        """
        Calculates field of ideal coil at position x
        
        Input: 
        """
        return hlx.ideal(i,r,x)
    def B(self,x):
        """
        Calculate the full HELIX magnetic field at points x.

        Input:
         x: list of list of doubles, xyz coordinates i.e. ((x1,y1,z1),...,(xn,yn,zn))
        Returns: list of list of doubles, field values in cartesian coordinates i.e. ((Bx1,By1,Bz1),...,(Bxn,Byn,Bzn))         
        """
        return hlx.helix_magnet(x)
#    def __check_data(self):
#    def __set_data(self):
    def WriteData():
        """
        Writes current dataframes to file
        """
        # Make general dataframe
        general_names = 'Current CuComp AlComp'
        general_attr = [93.0,0.99674,0.99585]
        general_df = pd.DataFrame([general_attr],columns=general_names.split())
        general_df.to_csv('GeneralData.csv')

        # Make dataframe for coils
        coil_names = 'A B'
        coil_attr = 'Width X Y Z Theta Phi'
        coil_a = [.075818,0.,0.,-.359747,0.,0.]
        coil_b = [.075818,0.,0.,.359747,0.,0.]
        coil_df = pd.DataFrame([coil_a,coil_b],index=coil_names.split(),columns=coil_attr.split())
        coil_df.index.name = 'Coil'
        coil_df.to_csv('CoilData.csv')

        # Make dataframe for subcoils
        coil_hier1 = 'A A A B B B'.split()
        coil_hier2 = [1,2,3,1,2,3]
        subcoil_index = list(zip(coil_hier1,coil_hier2))
        subcoil_index = pd.MultiIndex.from_tuples(subcoil_index)
        subcoil_columns = ["Inner Radius","Outer Radius","Coil Turns","Nrho","Nz"]
        sc_a1 = [.203517,.215011,1995.9,8,52]
        sc_a2 = [.215011,.234632,5150.,8,32]
        sc_a3 = [.234632,.255224,5489.5,8,32]
        sc_b1 = [.203492,.214731,1976.,8,52]
        sc_b2 = [.214732,.234709,5110.,8,32]
        sc_b3 = [.234709,.256222,5486.7,8,32]
        subcoil_df = pd.DataFrame([sc_a1,sc_a2,sc_a3,sc_b1,sc_b2,sc_b3],index=subcoil_index,columns=subcoil_columns)
        subcoil_df.index.names = ['Coil','Subcoil']
        subcoil_df.to_csv('SubcoilData.csv')

        


 
