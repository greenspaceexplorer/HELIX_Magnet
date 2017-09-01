MODULE COIL_DATA
  !-------------------------------------------------------------------------------
  !   Contains data for magnet coils.
  !
  !   NC: integer, number of coils
  !   I: double, coil currents
  !   W: double(NC), width of coils
  !   NS: integer, number of subcoils in each coil
  !   M: double(NC, NS), number of turns in each subcoil
  !   RI: double(NC, NS), inner radius of subcoils
  !   RO: double(NC, NS), outer radius of subcoils
  !   NRHO: integer(NC, NS), number of element divisions in rho direction in
  !         subcoil (coil index, subcoil index)
  !   NZ: integer(NC, NS), number of element divisions in z direction in subcoil
  !       (coil index, subcoil index)
  !   D: double(NC, 3), displacement of center of coil from origin
  !   THETA: double(NC), coil tilt in theta
  !   PHI: double(NC), coil rotation in phi
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
  LOGICAL :: DATA_SET = .FALSE.
  INTEGER :: NC,NS
  INTEGER,ALLOCATABLE :: NRHO(:,:),NZ(:,:)
  REAL*8 :: I
  REAL*8,ALLOCATABLE :: D(:,:)
  REAL*8,ALLOCATABLE :: W(:),THETA(:),PHI(:)
  REAL*8,ALLOCATABLE :: M(:,:),RI(:,:),RO(:,:)
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE SET_COIL_DATA(INC,II,ITHERMCOMP,IW,INS,IM,IRI,IRO,INRHO,INZ,ID,ITHETA,IPHI)
    ! Sets variables and allocates memory. ITHERMCOMP is the only new variable.
    !
    ! ITHERMCOMP: double(2), thermal contraction factors for Cu and Al respectively
    IMPLICIT NONE
    INTEGER :: J
    INTEGER :: INC,INS
    INTEGER :: INRHO(INC,INS),INZ(INC,INS)
    REAL*8 :: DELTAZ
    REAL*8 :: II,ITHERMCOMP(2)
    REAL*8 :: IW(INC),ITHETA(INC),IPHI(INC)
    REAL*8 :: ID(INC,3)
    REAL*8 :: IM(INC,INS),IRI(INC,INS),IRO(INC,INS)
    !  Reset memory if data has previously been set
    IF(DATA_SET)THEN
       CALL DEALLOC_MEMORY()
    ENDIF
    !   Allocate memory for arrays
    ALLOCATE(D(INC,3))
    ALLOCATE(NRHO(INC,INS),NZ(INC,INS))
    ALLOCATE(W(INC),THETA(INC),PHI(INC))
    ALLOCATE(M(INC,INS),RI(INC,INS),RO(INC,INS))
    !   Save values to variables
    NC = INC
    NS = INS
    NRHO = INRHO
    NZ = INZ
    I = II
    ! Apply copper thermal contraction to width and radius of coils
    W = IW*ITHERMCOMP(1)
    RI = IRI*ITHERMCOMP(1)
    RO = IRO*ITHERMCOMP(1)
    THETA = ITHETA
    PHI = IPHI
    D = ID
    ! Last coil remains fixed. Additional coils move along z-axis due to thermal
    ! contraction of aluminum.
    DO J=1,NC
       IF(J.NE.NC)THEN
          DELTAZ = ABS(D(J,3)-D(NC,3)) ! Distance to fixed coil
          D(J,3) = ID(J,3)+DELTAZ*(1.D0-ITHERMCOMP(2))
       ENDIF
    ENDDO
    M = IM
    DATA_SET = .TRUE.
    !
    RETURN
  END SUBROUTINE SET_COIL_DATA
  !-------------------------------------------------------------------------------
  SUBROUTINE DEALLOC_MEMORY()
    ! Clears memory allocated for arrays
    IMPLICIT NONE
    DEALLOCATE (NRHO,NZ)
    DEALLOCATE (D)
    DEALLOCATE (W,THETA,PHI)
    DEALLOCATE (M,RI,RO)
    RETURN
  END SUBROUTINE DEALLOC_MEMORY
  !-------------------------------------------------------------------------------
END MODULE COIL_DATA
