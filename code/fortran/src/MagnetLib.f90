!-------------------------------------------------------------------------------
! Library of functions and subroutines for calculating B field of the HELIX
! magnet.
!
! Author: Noah Green
! Created: August 21, 2017
!-------------------------------------------------------------------------------
!
SUBROUTINE HELIX_DATA(INC,II,ITHERMCOMP,IW,INS,IM,IRI,IRO,INRHO,INZ,ID,ITHETA,IPHI)
  ! Wraps SET_COIL_DATA subroutine from COIL_DATA module.
  USE COIL_DATA
  IMPLICIT NONE
  INTEGER :: INC,INS
  INTEGER :: INRHO(INC,INS),INZ(INC,INS)
  REAL*8 :: II,ITHERMCOMP(2)
  REAL*8 :: IW(INC),ITHETA(INC),IPHI(INC)
  REAL*8 :: ID(INC,3)
  REAL*8 :: IM(INC,INS),IRI(INC,INS),IRO(INC,INS)
  CALL SET_COIL_DATA(INC,II,ITHERMCOMP,IW,INS,IM,IRI,IRO,INRHO,INZ,ID,ITHETA,IPHI)
  RETURN
END SUBROUTINE HELIX_DATA
!
!-------------------------------------------------------------------------------
!
!!$SUBROUTINE CLEAR_DATA()
!!$  ! Wraps DEALLOC_DATA subroutine from COIL_DATA module. 
!!$  USE COIL_DATA
!!$  IMPLICIT NONE
!!$  CALL DEALLOC_DATA()
!!$  RETURN
!!$END SUBROUTINE CLEAR_DATA
!
!-------------------------------------------------------------------------------
!
SUBROUTINE HELIX_MAGNET(N,X,B)
  !   Calculates the magnetic field of HELIX magnet at position X.
  !
  !   N: integer, number of points at which field is calculated
  !   X: double(3,N), array of input points at which field is calculated
  !   B: double(3,N), field at input points
  !
  !   Iterates over each primary coil and adds resulting magnetic fields.
  !   Rotates and translates to coordinate system where primary coil is at origin
  !   and oriented in xy plane. Calculates field and rotates components back to
  !   original coordinate system prior to exiting subroutine.
  !
  USE COIL_DATA
  USE ROTATION
  IMPLICIT NONE
  INTEGER :: J,K,L
  !   Input/output variables
  INTEGER,INTENT(IN) :: N
  REAL*8,INTENT(IN) :: X(N,3)
  REAL*8,INTENT(INOUT) :: B(N,3)
  !   Local variables
  INTEGER :: NRHOCOIL(NS),NZCOIL(NS)
  REAL*8 :: XCOIL(3),BCOIL(3),DCOIL(3)
  REAL*8 :: MCOIL(NS),RICOIL(NS),ROCOIL(NS)

  IF (.NOT.DATA_SET)THEN ! Ensures data for magnet has been set
     ERROR STOP "FORTRAN ERROR: MAGNET DATA NOT SET"
  ENDIF
  !
  DO J=1,NC ! Iterate over main coils
     DO L=1,NS ! Get subcoil data for a single primary coil
        MCOIL(L) = M(J,L)
        RICOIL(L) = RI(J,L)
        ROCOIL(L) = RO(J,L)
        NRHOCOIL(L) = NRHO(J,L)
        NZCOIL(L) = NZ(J,L)
     ENDDO
     DCOIL = (/D(J,1),D(J,2),D(J,3)/)
     !
     IF((THETA(J).EQ.0.D0).AND.(PHI(J).EQ.0.D0))THEN
        CONTINUE
     ELSE
        CALL RMATRIX(THETA(J),PHI(J)) ! Make rotation matrices
     ENDIF
     !
     DO K=1,N ! Iterate over coordinates
        BCOIL = (/0.D0,0.D0,0.D0/) ! Reset BCOIL for each coordinate
        XCOIL = (/X(K,1),X(K,2),X(K,3)/) ! Get field calculation position
        XCOIL = XCOIL - DCOIL ! Translate to coil origin
        IF((THETA(J).NE.0.D0).OR.(PHI(J).NE.0.D0))THEN
           CALL CW(XCOIL) ! Rotate to coil coordinates
           CALL COIL(I,W(J),NS,MCOIL,RICOIL,ROCOIL,NRHOCOIL,NZCOIL,XCOIL,BCOIL)
           CALL CCW(BCOIL) ! Rotate field to original coordinates
        ELSE
           CALL COIL(I,W(J),NS,MCOIL,RICOIL,ROCOIL,NRHOCOIL,NZCOIL,XCOIL,BCOIL)
        ENDIF
        !       Add up field contributions
        B(K,1) = B(K,1) + BCOIL(1)
        B(K,2) = B(K,2) + BCOIL(2)
        B(K,3) = B(K,3) + BCOIL(3)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE HELIX_MAGNET
!
!-------------------------------------------------------------------------------
!
SUBROUTINE COIL(I,W,NS,M,RI,RO,NRHO,NZ,X,B)
  !   Calculates the magnetic field of a compound coil positioned at origin.
  !
  !   I: double, coil current
  !   W: double, coil width
  !   NS: integer, number of subcoils
  !   M: double(NS), number of turns in each subcoil
  !   RI: double(NS), inner radius of subcoils
  !   RO: double(NS), outer radius of subcoils
  !   NRHO: integer(NS), number of element divisions in rho direction
  !   NZ: integer(NS), number of element divisions in z direction
  !   X: double(3), position where field is calculated
  !   B: double(3), field at position X
  !
  !   Primary coil is broken into N subcoils to account for different turn
  !   densities. Each subcoil is broken into NRHO(ITER)*NZ(ITER) ideal current
  !   loops with current equal (coil current)*(turns in subcoil)
  !   *(subcoil cross sectional area)/(NRHO(ITER)*NZ(ITER)). Magnetic field of
  !   each ideal current loop is calculated with IDEAL(I,R,X,B) after accounting
  !   for position of the loop in the primary coil. Sum of magnetic fields of all
  !   ideal current loops in all subcoils gives the field of the primary coil.
  !
  IMPLICIT NONE
  INTEGER :: J,K,L
  INTEGER :: NS
  INTEGER :: NRHO(NS),NZ(NS)
  REAL*8 :: I,W
  REAL*8 :: X(3),B(3),XIC(3)
  REAL*8 :: M(NS),RI(NS),RO(NS)
  REAL*8 :: BSUB(3),BTEMP(3)
  REAL*8 :: IIC,WRHO,DRHO,DZ,RIC,ZIC
  !
  DO J = 1,NS
     IIC = I*M(J)/(DBLE(NRHO(J)*NZ(J))) ! Current in ideal current loop
     WRHO = RO(J)-RI(J) ! Width of subcoil in rho direction
     DRHO = WRHO/(DBLE(NRHO(J))) ! Rho step size
     DZ = W/(DBLE(NZ(J)))  ! Z step size
     BSUB(1) = 0.D0  ! Initialize subcoil field
     BSUB(2) = 0.D0
     BSUB(3) = 0.D0
     DO K = 1,NRHO(J)
        RIC = RO(J)-DRHO*(DBLE(K) - 0.5D0) ! Ideal coil radius
        DO L = 1,NZ(J)
           ZIC = (DBLE(L) - 0.5D0)*DZ - W/2.D0 ! Ideal coil z-displacement
           XIC = (/X(1),X(2),X(3)-ZIC/)
           CALL IDEAL(IIC,RIC,XIC,BTEMP)
           BSUB = BSUB+BTEMP
        ENDDO
     ENDDO
     B = B+BSUB
  ENDDO
  RETURN
END SUBROUTINE COIL
!
!-------------------------------------------------------------------------------
!
SUBROUTINE IDEAL(I,R,X,B)
  !   Calculates the magnetic field of an ideal current loop positioned at origin.
  !
  !   I: double, current
  !   R: double, radius
  !   X: double(3), position where field is calculated
  !   B: double(3), field at position X
  !
  !   See Smythe 1950 pp 271 for details.
  !
  IMPLICIT NONE
  REAL*8 :: I,R,RHO,PI,X(3),B(3) ! RHO = distance from origin in xy plane
  REAL*8 :: D1,D2,N1,N2,KSQ,KCSQ,KC,ELK,ELE,SINT,COST,BRHO,BB
  REAL*8 :: ELLIPI
  !
  PI = 4.D0*DATAN(1.D0)
  RHO = SQRT(X(1)*X(1)+X(2)*X(2))
  !
  D1 = (R+RHO)*(R+RHO)+X(3)*X(3)
  BB = 2.D-7*I/SQRT(D1)

  IF(RHO.EQ.0.D0)THEN
     B(1) = 0.D0
     B(2) = 0.D0
     B(3) = BB*PI*R*R/D1
  ELSE
     D2 = (R-RHO)*(R-RHO)+X(3)*X(3)
     N1 = R*R+RHO*RHO+X(3)*X(3)
     N2 = R*R-RHO*RHO-X(3)*X(3)
     KSQ = 4.D0*R*RHO/D1
     KCSQ = 1.D0-KSQ
     KC = SQRT(KCSQ)
     ELK = ELLIPI(KC,1.D0,1.D0,1.D0)
     ELE = ELLIPI(KC,1.D0,1.D0,KCSQ)
     !
     SINT = X(2)/RHO
     COST = X(1)/RHO
     !
     BRHO = (BB*X(3))/RHO*(N1*ELE/D2-ELK)
     B(1) = BRHO*COST
     B(2) = BRHO*SINT
     B(3) = BB*(ELK+N2*ELE/D2)
  ENDIF
  RETURN
END SUBROUTINE IDEAL
!
!-------------------------------------------------------------------------------
!
REAL*8 FUNCTION ELLIPI(QQC,PP,AA,BB)
  !  Function to evaluate the generalized complete elliptic integral
  !  cel(k_c,p,a,b) with QQC=k_c, PP=p, AA=a, and BB=b.  The complete
  !  elliptic integrals of the first (K) and second (E) kind are given by
  !  K(k)=cel(k_c,1.,1.,1.) and E(k)=cel(k_c,1.,1.,k_c**2), where
  !  k_c=sqrt(1-k*k).  Code is taken from
  !  Numerical Recipies, Cambridge University Press, 1988, p. 187, and
  !  modified to run in double precision.
  !
  IMPLICIT NONE
  REAL*8 :: QQC,PP,AA,BB,A,B,P,E,EM,QC,F,Q,G,CA,PIO2
  PARAMETER(CA=.0003, PIO2=1.5707963268)
  !  	The desired accuracy is CA*CA
  !
  IF(QQC.EQ.0.)THEN
     !    +SLN modifications made to track down error in CEL during execution
     WRITE(*,*) 'Failure in CEL'
     ELLIPI=0.
     return
  endif
  !    -SLN
  QC=ABS(QQC)
  A=AA
  B=BB
  P=PP
  E=QC
  EM=1.
  !
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
1 F=A
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
  !
  ELLIPI=PIO2*(B+A*EM)/(EM*(EM+P))
  !
  RETURN
END FUNCTION ELLIPI
!-------------------------------------------------------------------------------
