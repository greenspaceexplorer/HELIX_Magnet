c used in calculating HEAT magnetic field 5/94 SLN

c	OPTIONS /EXTEND_SOURCE/CHECK=(ALL)
	SUBROUTINE HEAT_MAG(XIN,YIN,ZIN,BRHO,BZ,B)
c
c  originally subroutine SMILI_MAG2
C 
C FUNCTIONAL DESCRIPTION:	
C 
C    Returns the value of the HEAT magnetic field at any point in space.
C 
C DUMMY ARGUMENTS:
C 
C    XIN,YIN,ZIN    The coordinates, in cm., of the point to get the field 
C 			for.  
C		    The origin is in the center of the bore, and the coils 
c			are centered about the Z axis.
C
C    BRHO,BZ	    The component values of the B-field at the above point in 
C			gauss.
C
C    B		    The magnitude of the B-field at the above point in gauss.
C 
C IMPLICIT INPUTS:
C 
C    none
C 
C IMPLICIT OUTPUTS:
C 
C    none
C 
C 
C SIDE EFFECTS:
C 
C    none
C 
C 
*********************************************************************
*  8/18/93:  -ZA changed from 23.304 to 23.360
*     
*            -RAD_IN_A, RAD_IN_B changed from 18.6817 to 18.682 
*
*            -"DATA NRHO,NZ /16,8/" changed to
*             "DATA NRHO,NZ /8,16/"
*
*            -Radial dimensions multiplied by thermal contraction
*             factor for copper (THERM_COMP_CU); axial dimensions
*             multiplied by thermal contraction factor for
*             aluminum (THERM_COMP_AL).  These factors describe
*             the fractional contraction of the material when it
*             is cooled from room temperature down to liquid helium
*             temp, and were taken from "A Compendium of the Properties 
*             of Materials at Low Temperature (Phase I)", Wright Air
*             Development Division Technical Report 60-56,
*
*  These changes were made by N. Greene, Boston University
*  (nrg@alvarez.bu.edu).  
*
*********************************************************************
c
c    Modified for use in calculating the magnetic field of the HEAT 
c	magnet. 
c    Note that the wrapping density of the coils differs along their 
c       width, requiring the superposition of 6 "coils" to emulate 
c	the actual two.
c
c    The "1" coil is the inside (close to coil center) part of the coil,
c	the "2" coil is the middle part, and "3" the outer part.
c
c    Changed W contraction factor from Greene's code from aluminum to copper
c
c    Magnet coil divided into 3 winding densities, A1,A2,A3,B1,B2,B3
c    Coils A and B have been split in to A1, A2, A3 and B1, B2, B3 each with 
c	their own location and wrapping parameters.
c
c    Chose segmentation of 
c	NRHO1,NZ1 /8,52/		!approximates square sections
c	NRHO2,NZ2 /8,32/		!approximates square sections
c	NRHO3,NZ3 /8,32/		!approximates square sections
c    after studying effect of increasing segmentation.
c
c    Note B coil is closest to stack, Bfield points from A to B.
c      Cryostat fixed end is stack end (thermal contractions do not effect 
c      position of B coil but DO effect A coil position).
c
c    HEAT coordinates map onto magnet coordinates like this: xHEAT = zMAG.
c       zHEAT is up.
c
c    5/94   Scott Nutter, Univ. of Michigan 
c	nutter@tarle3.physics.lsa.umich.edu
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	IMPLICIT NONE
	REAL BZ,B,R,XIN,YIN,ZIN
	REAL BRHO,X,Y,Z
	REAL I				!current
	real W				!coil width (in Z)
	real RAD_OUT_A1,RAD_IN_A1,RAD_OUT_B1,RAD_IN_B1
	real RAD_OUT_A2,RAD_IN_A2,RAD_OUT_B2,RAD_IN_B2
	real RAD_OUT_A3,RAD_IN_A3,RAD_OUT_B3,RAD_IN_B3
	real coilsep
	REAL BRHO_A1,BZ_A1
	REAL BRHO_A2,BZ_A2
	REAL BRHO_A3,BZ_A3
	real BRHO_B1,BZ_B1
	real BRHO_B2,BZ_B2
	real BRHO_B3,BZ_B3
	real NTURN_A1,B_A1,Z_A1,X_A1,Y_A1
	real NTURN_A2,B_A2,Z_A2,X_A2,Y_A2
	real NTURN_A3,B_A3,Z_A3,X_A3,Y_A3
	real NTURN_B1,B_B1,Z_B1,X_B1,Y_B1
	real NTURN_B2,B_B2,Z_B2,X_B2,Y_B2
	real NTURN_B3,B_B3,Z_B3,X_B3,Y_B3
        REAL THERM_COMP_CU, THERM_COMP_AL 
	DATA Z_A1,Z_B1 /-35.9747,35.9747/
	DATA Z_A2,Z_B2 /-35.9747,35.9747/
	DATA Z_A3,Z_B3 /-35.9747,35.9747/
	DATA X_A1,X_B1,Y_A1,Y_B1 /0.0,0.0,0.0,0.0/
	DATA X_A2,X_B2,Y_A2,Y_B2 /0.0,0.0,0.0,0.0/
	DATA X_A3,X_B3,Y_A3,Y_B3 /0.0,0.0,0.0,0.0/
	INTEGER*4 NRHO1,NZ1
	INTEGER*4 NRHO2,NZ2
	INTEGER*4 NRHO3,NZ3
	DATA RAD_OUT_A1,RAD_IN_A1,RAD_OUT_B1,RAD_IN_B1 	!outer/inner radii
     &       /21.5011,20.3517,21.4731,20.3492/
	DATA RAD_OUT_A2,RAD_IN_A2,RAD_OUT_B2,RAD_IN_B2 	!outer/inner radii
     &       /23.4632,21.5011,23.4709,21.4732/
	DATA RAD_OUT_A3,RAD_IN_A3,RAD_OUT_B3,RAD_IN_B3 	!outer/inner radii
     &       /25.5524,23.4632,25.6222,23.4709/
	DATA W /7.5818/
        DATA THERM_COMP_CU, THERM_COMP_AL /0.99674, 0.99585/ 
	DATA NTURN_A1,NTURN_B1 /1995.9,1976.0/
	DATA NTURN_A2,NTURN_B2 /5150.0,5110.0/
	DATA NTURN_A3,NTURN_B3 /5489.5,5486.7/
	DATA I /90.09/
	DATA NRHO1,NZ1 /8,52/		!approximates square sections
	DATA NRHO2,NZ2 /8,32/		!approximates square sections
	DATA NRHO3,NZ3 /8,32/		!approximates square sections
c	common/nloops/nrho1,nz1,nrho2,nz2,nrho3,nz3
C
c Do B coil first. Recall this coil's location is fixed in cryostat.
c
C Do "B1" coil...
C
	X = XIN - X_B1
	Y = YIN - Y_B1
	Z = ZIN - Z_B1
	CALL EXCOIL(RAD_OUT_B1*THERM_COMP_CU,RAD_IN_B1*THERM_COMP_CU,
     /          W*THERM_COMP_CU,I,NTURN_B1,NRHO1,NZ1,
     /          X,Y,Z,BRHO_B1,BZ_B1,B_B1)
C
C Do "B2" coil...
C
	X = XIN - X_B2
	Y = YIN - Y_B2
	Z = ZIN - Z_B2
	CALL EXCOIL(RAD_OUT_B2*THERM_COMP_CU,RAD_IN_B2*THERM_COMP_CU,
     /          W*THERM_COMP_CU,I,NTURN_B2,NRHO2,NZ2,
     /          X,Y,Z,BRHO_B2,BZ_B2,B_B2)
C
C Do "B3" coil...
C
	X = XIN - X_B3
	Y = YIN - Y_B3
	Z = ZIN - Z_B3
	CALL EXCOIL(RAD_OUT_B3*THERM_COMP_CU,RAD_IN_B3*THERM_COMP_CU,
     /          W*THERM_COMP_CU,I,NTURN_B3,NRHO3,NZ3,
     /          X,Y,Z,BRHO_B3,BZ_B3,B_B3)
C
c  Now do A coil, taking into effect change in position of coil 
c        due to thermal contraction. Contraction length is coil separation.
C
C Do "A1" coil...
C
 	COILSEP=Z_B1-Z_A1
c
C	TYPE *,' X,Y,Z:',XIN,YIN,ZIN
	X = XIN - X_A1
	Y = YIN - Y_A1
	Z = ZIN - (Z_A1+coilsep*(1.-THERM_COMP_AL))
	CALL EXCOIL(RAD_OUT_A1*THERM_COMP_CU,RAD_IN_A1*THERM_COMP_CU,
     /          W*THERM_COMP_CU,I,NTURN_A1,NRHO1,NZ1,
     /          X,Y,Z,BRHO_A1,BZ_A1,B_A1)
C
C Do "A2" coil...
C
 	COILSEP=ABS(Z_B2-Z_A2)
c
C	TYPE *,' X,Y,Z:',XIN,YIN,ZIN
	X = XIN - X_A2
	Y = YIN - Y_A2
	Z = ZIN - (Z_A2+coilsep*(1.-THERM_COMP_AL))
	CALL EXCOIL(RAD_OUT_A2*THERM_COMP_CU,RAD_IN_A2*THERM_COMP_CU,
     /          W*THERM_COMP_CU,I,NTURN_A2,NRHO2,NZ2,
     /          X,Y,Z,BRHO_A2,BZ_A2,B_A2)
C
C Do "A3" coil...
C
 	COILSEP=ABS(Z_B3-Z_A3)
c
C	TYPE *,' X,Y,Z:',XIN,YIN,ZIN
	X = XIN - X_A3
	Y = YIN - Y_A3
	Z = ZIN - (Z_A3+coilsep*(1.-THERM_COMP_AL))
	CALL EXCOIL(RAD_OUT_A3*THERM_COMP_CU,RAD_IN_A3*THERM_COMP_CU,
     /          W*THERM_COMP_CU,I,NTURN_A3,NRHO3,NZ3,
     /          X,Y,Z,BRHO_A3,BZ_A3,B_A3)
c
c now calculate field....
	BRHO = BRHO_A1 + BRHO_A2 + BRHO_A3 + BRHO_B1 + BRHO_B2 + BRHO_B3 
	BZ = BZ_A1 + BZ_A2 + BZ_A3 + BZ_B1 + BZ_B2 + BZ_B3
C	TYPE *,' BZ_A1:',BZ_A1,' BZ_B1:',BZ_B1
C	TYPE *,' BZ_A2:',BZ_A2,' BZ_B2:',BZ_B2
C	TYPE *,' BZ_A3:',BZ_A3,' BZ_B3:',BZ_B3
C	TYPE *,' BRHO_A1:',BRHO_A1,' BRHO_B1:',BRHO_B1
C	TYPE *,' BRHO_A2:',BRHO_A2,' BRHO_B2:',BRHO_B2
C	TYPE *,' BRHO_A3:',BRHO_A3,' BRHO_B3:',BRHO_B3
	B = SQRT(BRHO**2 + BZ**2)
	RETURN
	END
C
C----------------------------------------------------------------------
C
	SUBROUTINE EXCOIL(RHOOUT,RHOIN,W,I,NTURNS,NRHO,NZ,X,Y,Z,BRHOEC,BZEC,
     1  MODBEC)
C
C	Subroutine to break a magnet coil of rectangular cross section into
C	finite elements.
C
C		RHOOUT= Outer radius of coil (cm)
C
C		RHOIN= Inner radius of coil (cm)
C
C		W= Width of coil (Z extent) (cm)
C
C		I= Current in coil (Amperes)
C
C		NTURNS= Number of turns on coil
C
C		NRHO= Number of slices in radius (integer)
C
C		NZ= Number of slices in Z (integer)
C
C		X,Y,Z= Field point (cm)
C
C		BRHOEC,BZEC= Field components for Extended Coil (gauss)
C
C		MODBEC= Total field magnitude for Extended Coil (gauss)
C
C		Total number of elements (simple loop coils)=NRHO*NZ
C
C		Andrew Tomasch, November 1989
C
C----------------------------------------------------------------------
C
	IMPLICIT NONE
C
	REAL RHOOUT,RHOIN,W,I,NTURNS,X,Y,Z,BRHOEC,BZEC,MODBEC
	REAL H,BRHO,BZ,MODB,XNRHO,XNZ,XNEL,XN,XM
	REAL AEL,ZEL
C
	LOGICAL*4 LCK
	INTEGER*4 NRHO,NZ,NEL,N,M
CPAR$	PRIVATE N,M,BRHO,BZ,AEL,ZEL,XN,XM
CPAR$	CONTEXT_SHARED LCK
C
C	Calculate derived parameters
C
	LCK = .FALSE.
C
	H=RHOOUT-RHOIN
	NEL=NRHO*NZ
	XNRHO=DFLOAT(NRHO)
	XNZ=DFLOAT(NZ)
	XNEL=DFLOAT(NEL)
C
C	Initialize field components
C
	BRHOEC=0.
	BZEC=0.
	MODBEC=0.
C
C	Superpose field components for all elements (N,M)
C
CPAR$	DO_PARALLEL NRHO/NWORKERS()
	DO N=1,NRHO
C
	  DO M=1,NZ
C
	    XN = DFLOAT(N)
	    XM = DFLOAT(M)
C
	    AEL = RHOOUT-(0.5*(H/XNRHO)*(2.*XN - 1.))
	    ZEL = -0.5*W+(0.5*(W/XNZ)*(2.*XM - 1.))
C
	    CALL COILB(I,AEL,(NTURNS/XNEL),X,Y,(Z-ZEL),BRHO,BZ,MODB)
C
CPAR$	    LOCKON LCK
	    BRHOEC=BRHOEC+BRHO
	    BZEC=BZEC+BZ
CPAR$	    LOCKOFF LCK
C
	  ENDDO
C
	ENDDO
C
C	Calculate MODBEC for extended coil
C
	MODBEC=(BRHOEC*BRHOEC)+(BZEC*BZEC)
	MODBEC=SQRT(MODBEC)
C
	RETURN
	END
C
C
C
	SUBROUTINE COILB(I,A,NTURNS,X,Y,Z,BRHO,BZ,MODB)
C
C--------------------------------------------------------------------------
C
C	Subroutine to evaluate the magnetic field of an ideal coil of radius
C	A (cm) carrying current I (amperes) with NTURNS total turns.  The
C	coil lies in the X-Y plane with the coil center at the origin.  The
C	field is calculated in cylindrical coordinates using the formulation
C	of Smythe.  The components BRHO and BZ are evaluated at the point
C	X,Y,Z (cm) and are reported in Gauss.  The magnitude of B is also
C	reported as MODB.  The complete elliptic integrals of the first (K)
C	and second (E) kind are evaluated using the generalized complete
C	elliptic integral function CEL(QQC,PP,AA,BB), where
C	K(k)=CEL(KC,1.,1.,1.) and E(k)=CEL(KC,1.,1.,KC_2), where 
C	KC_2=k_c**2=1-k**2, and KC=k_c=SQRT(KC_2).  KARG_2=k**2.  The
C	CEL function routine is from Numerical Recipies, page 187.
C
C--------------------------------------------------------------------------
C
	IMPLICIT NONE
	REAL I,A,NTURNS,X,Y,Z,RHO,KARG_2,CEL,C,DENOM,K,E,BRHO,BZ,MODB
	REAL KC,KC_2
C
C	C= speed of light in cm/s
	C=2.99792458E10
C
	RHO=SQRT((X*X)+(Y*Y))
	KARG_2=(4.*A*RHO)/((A+RHO)*(A+RHO)+(Z*Z))
C
	KC_2=1.-KARG_2
	KC=SQRT(KC_2)
C
	K=CEL(KC,1.,1.,1.)
	E=CEL(KC,1.,1.,KC_2)
c+SLN modifications made to track down error in CEL during execution
        if(k.eq.0.) write(*,*) 'Possible failure in CEL, k=0'
        if(e.eq.0.) write(*,*) 'Possible failure in CEL, e=0'
 	if(k.eq.0. .or. e.eq.0.) then
	  write(*,*) 'x,y,z=',x,y,z
	  write(*,*) 'i,a,nturns=',i,a,nturns
	endif
c-SLN
C
	DENOM=SQRT((A+RHO)*(A+RHO)+(Z*Z))
C
C	Factor of 3.E9 converts amperes to statamperes
C
	IF(RHO.EQ.0.)THEN
		BRHO=0.
		BZ=((A*A)-(RHO*RHO)-(Z*Z))/((A-RHO)*(A-RHO)+(Z*Z))
		BZ=E*BZ
		BZ=BZ+K
		BZ=(I*NTURNS)*(2.99792458E9/C)*(2./DENOM)*BZ
		MODB=BZ
		RETURN
C
	ELSE
C
		BRHO=((A*A)+(RHO*RHO)+(Z*Z))/((A-RHO)*(A-RHO)+(Z*Z))
		BRHO=E*BRHO
		BRHO=BRHO-K
		BRHO=(I*NTURNS)*(2.99792458E9/C)*(2.*Z/(DENOM*RHO))*BRHO
C
		BZ=((A*A)-(RHO*RHO)-(Z*Z))/((A-RHO)*(A-RHO)+(Z*Z))
		BZ=E*BZ
		BZ=BZ+K
		BZ=(I*NTURNS)*(2.99792458E9/C)*(2./DENOM)*BZ
C
		MODB=SQRT((BRHO*BRHO)+(BZ*BZ))
C
		RETURN
	ENDIF
C
	END
C
C
C
	REAL FUNCTION CEL(QQC,PP,AA,BB)
C
C---------------------------------------------------------------------------
C
C	Function to evaluate the generalized complete elliptic integral
C	cel(k_c,p,a,b) with QQC=k_c, PP=p, AA=a, and BB=b.  The complete
C	elliptic integrals of the first (K) and second (E) kind are given by
C	K(k)=cel(k_c,1.,1.,1.) and E(k)=cel(k_c,1.,1.,k_c**2), where
C	k_c=sqrt(1-k*k).  Code is taken from
C	Numerical Recipies, Cambridge University Press, 1988, p. 187, and
C	modified to run in double precision.
C
C----------------------------------------------------------------------------
C
	IMPLICIT NONE
	REAL QQC,PP,AA,BB,A,B,P,E,EM,QC,F,Q,G,CA,PIO2
C
	PARAMETER(CA=.0003, PIO2=1.5707963268)
C	The desired accuracy is CA*CA
C
	IF(QQC.EQ.0.)then
c+SLN modifications made to track down error in CEL during execution
          write(*,*) 'Failure in CEL'
          cel=0.
          return
        endif
c-SLN 
	QC=ABS(QQC)
	A=AA
	B=BB
	P=PP
	E=QC
	EM=1.
C
	IF(P.GT.0.)THEN
		P=SQRT(P)
		B=B/P
	ELSE
		F=QC*QC
		Q=1.-F
		G=1.-P
		F=F-P
		Q=Q*(B-A*P)
		P=SQRT(F/G)
		A=(A-B)/G
		B=-Q/(G*G*P)+A*P
	ENDIF
C
 1	F=A
	A=A+B/P
	G=E/P
	B=B+F*G
	B=B+B
	P=G+P
	G=EM
	EM=QC+EM
	IF(ABS(G-QC).GT.G*CA)THEN
		QC=SQRT(E)
		QC=QC+QC
		E=QC*EM
		GOTO 1
	ENDIF
C
	CEL=PIO2*(B+A*EM)/(EM*(EM+P))
C
	RETURN
	END

