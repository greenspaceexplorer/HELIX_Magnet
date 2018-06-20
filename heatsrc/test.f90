PROGRAM HELIXTEST
  IMPLICIT NONE
  REAL*8 I,A
  REAL*8 NTURNS
  REAL*8 X,Y,Z
  REAL*8 BRHO,BZ,MODB
  REAL*8 ST,SP,CT,CP
  REAL*8 BX,BY,BZ

! Test individual coil subroutine
  I = 93.D0
  A = 25.5224D0
  NTURNS = 1.D0
  X = 0.D0
  Y = 0.D0
  Z = 10.D0

  CALL COILB(I,A,NTURNS,X,Y,Z,BRHO,BZ,MODB)

  CT = 

  PRINT *,"Single coil test:"

  PRINT *,"X,Y,Z = ",X,Y,Z
  PRINT *,"BX,BY,BZ = ",BX,BY,BZ
 
! Test compound coil subroutine


! Test field for magnet with no shifts


  STOP
END PROGRAM HELIXTEST
