PROGRAM HELIXTEST
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER ::     &
      VOUT = "(5X,A14,3ES12.4E2,A1)",&
      BOUT1 = "(BX BY BZ) = (",     &
      XOUT1 = "( X  Y  Z) = (",     &
      VOUT2 = ")"
  REAL*8 :: I,A
  REAL*8 :: NTURNS
  REAL*8 :: X,Y,Z,RHO
  REAL*8 :: BRHO,BZ,MODB
  REAL*8 :: BX,BY

! Test individual coil subroutine
  WRITE(*,*) "Single coil test:"
  I = 93.D0
  A = 25.5224D0
  NTURNS = 1.D0
  
  WRITE(*,*)
  WRITE(*,*) "    Current = ",I,"A"
  WRITE(*,*) "    Radius = ",A,"cm" 
  ! First test with on-axis coordinates

  X = 0.D0
  Y = 0.D0
  Z = 10.D0

  RHO = SQRT(X*X+Y*Y)

  CALL COILB(I,A,NTURNS,X,Y,Z,BRHO,BZ,MODB)

  IF(RHO.EQ.0.D0)THEN
      BX = BRHO
      BY = BRHO
  ELSE
      BX = BRHO*X/RHO
      BY = BRHO*Y/RHO
  ENDIF
  WRITE(*,*)
  WRITE(*,VOUT) XOUT1,X,Y,Z,VOUT2
  WRITE(*,VOUT) BOUT1,BX,BY,BZ,VOUT2
  ! Next test with off axis coordinates

  X = 5.D0
  Y = 7.D0
  Z = 10.D0

  RHO = SQRT(X*X+Y*Y)

  CALL COILB(I,A,NTURNS,X,Y,Z,BRHO,BZ,MODB)

  IF(RHO.EQ.0.D0)THEN
      BX = BRHO
      BY = BRHO
  ELSE
      BX = BRHO*X/RHO
      BY = BRHO*Y/RHO
  ENDIF
  WRITE(*,*)
  WRITE(*,VOUT) XOUT1,X,Y,Z,VOUT2
  WRITE(*,VOUT) BOUT1,BX,BY,BZ,VOUT2

! Test compound coil subroutine


! Test field for magnet with no shifts


  STOP
END PROGRAM HELIXTEST
