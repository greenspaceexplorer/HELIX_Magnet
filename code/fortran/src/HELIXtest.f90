PROGRAM HELIXTEST
  IMPLICIT NONE
  INTEGER J
  ! COMMON VARIABLES
  INTEGER NC,NS,N
  PARAMETER (NC=2,NS=3,N=2)
  REAL*8 I  
  ! HEAT_MAG VARIABLES
  REAL*8 BHEATRHO,BHEATZ,BHEAT,XIN,YIN,ZIN
  REAL*8 RHO,SPHI,CPHI,BHEATX,BHEATY
  ! HELIX_DATA VARIABLES
  INTEGER NRHO(NC,NS),NZ(NC,NS)
  REAL*8 W(NC),D(NC,3),THETA(NC),PHI(NC)  
  REAL*8 M(NC,NS),RI(NC,NS),RO(NC,NS)
  REAL*8 THERMCOMP(2)
  ! HELIX_MAGNET VARIABLES
  REAL*8 X(NC,3),B(NC,3)
  ! COIL VARIABLES
  REAL*8 XC(3),BC(3)  
  INTEGER NRHOC(NS),NZC(NS)  
  REAL*8 MC(NS),RIC(NS),ROC(NS)
  !
  !-----INPUT MAGNET COIL DATA-----
  I = 93.D0
  W = (/0.075818,0.075818/)
  THETA = (/0.D0,0.D0/)
  PHI = (/0.D0,0.D0/)
  THERMCOMP = (/1.D0,1.D0/)
  !-----INPUT MAGNET SUBCOIL DATA-----
  DATA RI(1,1),RO(1,1),M(1,1),NRHO(1,1),NZ(1,1)/0.203517,0.215011,1995.9,8,52/
  DATA RI(1,2),RO(1,2),M(1,2),NRHO(1,2),NZ(1,2)/0.215011,0.234632,5150.0,8,32/
  DATA RI(1,3),RO(1,3),M(1,3),NRHO(1,3),NZ(1,3)/0.234632,0.255224,5489.5,8,32/
  DATA RI(2,1),RO(2,1),M(2,1),NRHO(2,1),NZ(2,1)/0.203492,0.214731,1976.0,8,52/
  DATA RI(2,2),RO(2,2),M(2,2),NRHO(2,2),NZ(2,2)/0.214732,0.234709,5110.0,8,32/
  DATA RI(2,3),RO(2,3),M(2,3),NRHO(2,3),NZ(2,3)/0.234709,0.256222,5486.7,8,32/
  DATA D(1,1),D(1,2),D(1,3)/0.0,0.0,-0.359747/
  DATA D(2,1),D(2,2),D(2,3)/0.0,0.0,0.359747/
  !-----INPUT COORDINATES WHERE FIELD WILL BE CALCULATED-----
  DATA X(1,1),X(1,2),X(1,3)/0.1,0.1,0.05/
  DATA X(2,1),X(2,2),X(2,3)/0.0,0.0,0.15/
  !-----COPY DATA FOR SINGLE COIL CALCULATION-----
  MC = (/M(1,1),M(1,2),M(1,3)/)
  RIC = (/RI(1,1),RI(1,2),RI(1,3)/)
  ROC = (/RO(1,1),RO(1,2),RO(1,3)/)
  NRHOC = (/NRHO(1,1),NRHO(1,2),NRHO(1,3)/)
  NZC = (/NZ(1,1),NZ(1,2),NZ(1,3)/)
  XC = (/X(1,1),X(1,2),X(1,3)/)
  BC = 0.D0
  !
  !------HEAT MAGNET CALCULATION-----
  PRINT *,RI(1,1)*100.D0,RI(1,2)*100.D0,RI(1,3)*100.D0
  PRINT *,RO(1,1)*100.D0,RO(1,2)*100.D0,RO(1,3)*100.D0
  ! CONVERT TO CM
  XIN = X(2,1)*100.D0
  YIN = X(2,2)*100.D0
  ZIN = X(2,3)*100.D0
  PRINT *,"X,Y,Z = ",XIN,YIN,ZIN
  ! CALL ORIGINAL HEAT_MAG SUBROUTINE
  CALL HEAT_MAG(XIN,YIN,ZIN,BHEATRHO,BHEATZ,BHEAT)
  ! CONVERT TO CARTESIAN COORDINATES
  PRINT *,"BRHO,BZ,B = ",BHEATRHO,BHEATZ,BHEAT
  RHO = SQRT(XIN*XIN+YIN*YIN)
  IF (RHO.NE.0.D0) THEN
     CPHI = XIN/RHO
     SPHI = YIN/RHO
     BHEATX = BHEAT*CPHI*1.D-4 ! FACTOR OF 10^-4 TO CONVERT TO TESLA
     BHEATY = BHEAT*SPHI*1.D-4
  ELSE
     BHEATX = BHEATRHO*1.D-4
     BHEATY = BHEATRHO*1.D-4
  ENDIF
  BHEATZ = BHEATZ*1.D-4
  !-----END HEAT MAGNET CALCULATION-----
  !
  !
  !
  !-----HELIX MAGNET CALCULATION: NO ROTATION-----
  B = 0.D0
  CALL HELIX_DATA(NC,I,THERMCOMP,W,NS,M,RI,RO,NRHO,NZ,D,THETA,PHI)
  CALL HELIX_MAGNET(N,X,B)

  WRITE(*,90)
  WRITE(*,100) X(2,1),X(2,2),X(2,3),BHEATX,BHEATY,BHEATZ
  DO J=1,N
     WRITE(*,100) X(J,1),X(J,2),X(J,3),B(J,1),B(J,2),B(J,3)
  ENDDO
  PRINT *,"-----TESTING COIL SUBROUTINE-----"
  CALL COIL(I,W(1),NS,MC,RIC,ROC,NRHOC,NZC,XC,BC)
  PRINT *,"BXYZ = ",BC(1),BC(2),BC(3)
90 FORMAT(19X,'X',17X,'Y',17X,'Z',16X,'BX',16X,'BY',16X,'BZ')
100 FORMAT(2X,6(6X,E12.6))

  BC = 0.D0
  BHEATRHO = 0.D0
  BHEATZ = 0.D0
  BHEAT = 0.D0
  CALL IDEAL(100.D0,.2D0,XC,BC)
  XC = XC*100.D0
  CALL COILB(100.D0,20.D0,1.D0,XC(1),XC(2),XC(3),BHEATRHO,BHEATZ,BHEAT)

  PRINT *,"COILB BZ = ",BHEATZ*1.D-4
  PRINT *,"IDEAL BZ = ",BC(3)


  STOP
END PROGRAM HELIXTEST
!
!----------------------------------------------
!
SUBROUTINE PRINT_ARRAY(ROW,COL,NAME)
  IMPLICIT NONE
  INTEGER ROW,COL,I,J
  REAL*8 :: NAME(ROW,COL)
  DO I=1,ROW
     WRITE(*,100) (NAME(I,J),J=1,COL)
  ENDDO
100 FORMAT(100g15.5)
END SUBROUTINE PRINT_ARRAY

