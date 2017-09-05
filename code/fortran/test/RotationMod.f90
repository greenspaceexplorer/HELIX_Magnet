MODULE ROTATION
  !-------------------------------------------------------------------------------
  ! Provides means to calculate, store, and utilize magnet rotation matrices.
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
  REAL*8 :: RCCW(3,3),RCW(3,3)
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE RMATRIX(TH,PH)
    !   Builds rotation matrices out of theta and phi angles. Stores rotation
    !   matrices to the 'ROTATIONS' common block.
    !
    !   TH: double, theta
    !   PH: double, phi
    !   RCCW: double(3,3), compound rotation matrix
    !   RCW: double(3,3), inverse rotation matrix
    !   RT: double(3,3), theta rotation matrix
    !   RP: double(3,3), phi rotation matrix
    !
    IMPLICIT NONE
    INTEGER :: I,J,K
    REAL*8,INTENT(IN) :: TH,PH
    REAL*8 :: RT(3,3),RP(3,3)
    REAL*8 :: CT,ST,CP,SP
    !
    CT = COS(TH)
    ST = SIN(TH)
    CP = COS(PH)
    SP = SIN(PH)
    !   Rotation about z-axis
    RP(1,1) = CP
    RP(1,2) = -1.D0*SP
    RP(1,3) = 0.D0
    RP(2,1) = SP
    RP(2,2) = CP
    RP(2,3) = 0.D0
    RP(3,1) = 0.D0
    RP(3,2) = 0.D0
    RP(3,3) = 1.D0
    !   Rotation about y-axis
    RT(1,1) = CT
    RT(1,2) = 0.D0
    RT(1,3) = ST
    RT(2,1) = 0.D0
    RT(2,2) = 1.D0
    RT(2,3) = 0.D0
    RT(3,1) = -1.D0*ST
    RT(3,2) = 0.D0
    RT(3,3) = CT
    !
    !   Left multiply RT by RP to get y,z axis rotation
    RCW = 0.D0
    RCCW = 0.D0
    DO I=1,3
       DO J=1,3
           DO K=1,3
             RCCW(I,J) = RCCW(I,J) + RP(I,K)*RT(K,J)
          ENDDO
          !   SO(3) rotation matrices -> RCW = RCCW^T
          RCW(J,I) = RCCW(I,J)
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE RMATRIX
  !
  !-------------------------------------------------------------------------------
  !
  SUBROUTINE CCW(X)
    !   Rotate vector counterclockwise first around y-axis, then around z-axis using
    !   rotation matrix RM. Run RMATRIX subroutine to set RM.
    !
    !   X: double(3), some arbitrary vector
    !
    IMPLICIT NONE
    INTEGER :: I,J
    REAL*8 :: X(3),XTEMP(3)
    !
    XTEMP = 0.D0
    DO I = 1,3
       DO J = 1,3
          XTEMP(I) = XTEMP(I) + RCCW(I,J)*X(J)
       ENDDO
    ENDDO
    X = XTEMP
    RETURN
  END SUBROUTINE CCW
  !
  !-------------------------------------------------------------------------------
  !
  SUBROUTINE CW(X)
    !   Rotate vector clockwise first around z-axis, then around y-axis using
    !   rotation matrix RI. Run RMATRIX subroutine to set RI.
    !
    !   X: double(3), some arbitrary vector
    !
    IMPLICIT NONE
    INTEGER :: I,J
    REAL*8 :: X(3),XTEMP(3)
    !
     XTEMP = 0.D0
    DO I = 1,3
       DO J = 1,3
          XTEMP(I) = XTEMP(I) + RCW(I,J)*X(J)
       ENDDO
    ENDDO
    X = XTEMP
    RETURN
  END SUBROUTINE CW
  !
  !-------------------------------------------------------------------------------
END MODULE ROTATION
