import helix_mag as hlx # This is the fortran library

class helix():
    def __init__(self,current,contraction,width,position
            ,theta,phi,turns,inner_radius,outer_radius,nrho,nz):
        """
        Initializes helix class object.

        Let nc be the number of coils and ns be the number of subcoils.

        Input:
            current = float
            contraction = list, len=2, (Cu, Al) for thermal contraction of metals
            width = list, len=nc, width of coils in meters
            position = list1(list2), len1 = nc, len2 = 3, origins of coils in meters
            theta = list, len=nc, theta tilt of coils in radians
            phi = list, len=nc, phi rotation of coils in radians
            turns = list1(list2), len1=nc, len2=ns, number of turns in each subcoil
            inner_radius = list1(list2), len1=nc, len2=ns, inner radius of subcoils
            outer_radius = list1(list2), len1=nc, len2=ns, outer radius of subcoils
            nrho = list1(list2), len1=nc, len2=ns, number radial FEA divisions in subcoils
            nz = list1(list2), len1=nc, len2=ns, number of FEA divisions in subcoils in z
        """
        self.i = current
        self.c = contraction
        self.w = width
        self.p = position
        self.th = theta
        self.ph = phi
        self.tn = turns
        self.ri = inner_radius
        self.ro = outer_radius
        self.nrho = nrho
        self.nz = nz

        # TODO: Check validity of input

        # Send coil data to fortran code        
        hlx.helix_data(self.i,self.c,self.w,self.tn,self.ri,self.ro,self.nrho,self.nz,self.p,self.th,self.ph)

#   TODO: make the following functions
#    def print_general(self):
#    def print_coil(self):
#    def print_subcoil(self):
#    def __check_data(self):
#    def __set_data(self):

    def set_coil_position(self,coil,position):
        """
        Changes position of coil.
        
        Input:
         coil: integer, coil index in spreadsheet
         position: list of doubles, (x,y,z) origin of coil
        """
        self.p[coil] = position
        hlx.helix_data(self.i,self.c,self.w,self.tn,self.ri,self.ro,self.nrho,self.nz,self.p,self.th,self.ph)
    def set_coil_angle(self,coil,angle):
        """
        Changes angle of coil.

        Input:
         coil: integer, coil index in spreadsheet
         angle: list of doubles, (theta,phi) angle of coil
        """
        self.th[coil] = angle[0]
        self.ph[coil] = angle[1]
        hlx.helix_data(self.i,self.c,self.w,self.tn,self.ri,self.ro,self.nrho,self.nz,self.p,self.th,self.ph)
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
