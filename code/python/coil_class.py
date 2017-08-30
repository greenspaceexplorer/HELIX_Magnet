import numpy as np
import math
from scipy import special

################################################################################
################################################################################

class coil:
    def __init__(self,current = 100.,radius = 1.,center = (0.,0.,0.),tilt = (0.,0.)):
        """
        Instantiates coil class

        Input:
            radius: float
                    Radius of coil
            center: (float,float,float)
                    (x,y,z) coordinates of origin
            tilt:   (float,float)
                    (theta,phi) orientation of normal vector wrt origin
        """

        self.__i = current
        self.__r = radius
        self.__center = center
        self.__tilt = tilt
        self.__isvalid()
        self.__brot = False
        if tilt[0] != 0. or tilt[1]!=0.:
            self.__brot = True
        self.__rotate = self.__make_rotation()
        self.__derotate = self.__rotate.transpose()

################################################################################

    def __isvalid(self):
        """
        Ensures valid coil is entered
        """
        self.__i = float(self.__i)
        self.__r = float(self.__r)

        if len(self.__center) != 3:
            raise Exception("Error: center tuple should have length 3")
        if len(self.__tilt) != 2:
            raise Exception("Error: Error: tilt tuple should have length 2")

        self.__tilt = np.array(list(map(float,self.__tilt)))
        self.__tilt[0]%=math.pi
        self.__tilt[1]%=(2.*math.pi)
        self.__center = np.array(list(map(float,self.__center)))

################################################################################

    def __str__(self):
        return "Current = {} amps".format(self.__i) + \
            "\nRadius = {} meters".format(self.__r) + \
            "\nCenter [x,y,z] = {}".format(self.__center) + \
            "\nNormal vec [theta,phi] = {}".format(self.__tilt)

################################################################################

    def __make_rotation(self):
        if self.__brot:
            cp = math.cos(self.__tilt[1])
            sp = math.sin(self.__tilt[1])
            ct = math.cos(self.__tilt[0])
            st = math.sin(self.__tilt[0])
            rotz = np.array([[cp,sp,0.],[-sp,cp,0.],[0.,0.,1.]])
            roty = np.array([[ct,0,-st],[0.,1.,0.],[st,0.,ct]])
            return roty.dot(rotz)
        else:
            return np.eye(3)

################################################################################

    def set_r(self,r):
        self.__r = float(r)

    def r(self):
        return self.__r

################################################################################

    def set_i(self,i):
        self.__i = float(i)

    def i(self):
        return self.__i

################################################################################

    def set_theta(self,theta):
        theta = float(theta)%math.pi
        self.__tilt[0] = theta
        if self.__tilt[0]!=0. or self.__tilt[1]!=0.:
            self.__brot = True
            self.__rotate = self.__make_rotation()

    def theta(self):
        return self.__tilt[0]

################################################################################

    def set_phi(self,phi):
        phi = float(phi)%(2.*math.pi)
        self.__tilt[1] = phi
        if self.__tilt[0]!=0. or self.__tilt[1]!=0.:
            self.__brot = True
            self.__rotate = self.__make_rotation()

    def phi(self):
        return self.__tilt[1]

################################################################################

    def __Brhophi0(self,rho,z):
        """
        Calculates components of B-field from coil in xy plane at origin at\\
         (rho,z)
        See Smythe 1950 pp 271 for details

        Input:
            rho: float
            z: float
        Returns: (float,float)
                 (Brho,Bz) components of magnetic field
        """
        rho = float(rho)
        z = float(z)
        d1 = (self.__r+rho)**2+z**2
        B = 2.e-7*self.__i/math.sqrt(d1)
        if rho == 0:
            return (0.,B*math.pi*self.__r**2/d1)
        else:
            d2 = (self.__r-rho)**2+z**2
            n1 = self.__r**2+rho**2+z**2
            n2 = self.__r**2-rho**2-z**2
            ksq = 4.*self.__r*rho/d1
            ell_k = special.ellipk(ksq)
            ell_e = special.ellipe(ksq)
            return (B*z/rho*(n1*ell_e/d2-ell_k),B*(ell_k+n2*ell_e/d1))

################################################################################

    def Bxyz(self,coords):
        """

        """
        if len(coords)!=3:
            raise Exception("Error: invalid number of coordinates")
        coords = np.array(coords)
        coords -= self.__center
        if self.__brot:
            coords = self.__rotate.dot(coords)


        rho = math.sqrt(coords[0]**2+coords[1]**2)
        Brho,Bz = self.__Brhophi0(rho,coords[2])
        B = np.array([Brho,0.,Bz])
        cphi = 1.
        sphi = 0.
        if (coords[0] != 0. or coords[1]!=0):
            cphi = coords[0]/math.sqrt(coords[0]**2+coords[1]**2)
            sphi = coords[1]/math.sqrt(coords[0]**2+coords[1]**2)
        m_cart = np.array([[cphi,-1.*sphi,0.],[sphi,cphi,0.],[0.,0.,1.]])
        return self.__derotate.dot(m_cart.dot(B))
