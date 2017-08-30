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
  LOGICAL :: SET_DATA = FALSE
  INTEGER :: NC,NS
  INTEGER,ALLOCATABLE :: NRHO(:,:),NZ(:,:)
  REAL*8 :: I,THERMCOMP(2)
  REAL*8,ALLOCATABLE :: D(:,:)
  REAL*8,ALLOCATABLE :: W(:),THETA(:),PHI(:)
  REAL*8,ALLOCATABLE :: M(:,:),RI(:,:),RO(:,:)
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE SET_COIL_DATA(INC,II,ITHERMCOMP,IW,INS,IM,IRI,IRO,INRHO,INZ,ID,ITHETA,IPHI)
    INTEGER :: INC,INS
    INTEGER :: INRHO(INC,INS),INZ(INC,INS)
    REAL*8 :: II,ITHERMCOMP(2)
    REAL*8 :: IW(INC),ITHETA(INC),IPHI(INC)
    REAL*8 :: ID(INC,3)
    REAL*8 :: IM(INC,INS),IRI(INC,INS),IRO(INC,INS)
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
    THERMCOMP = ITHERMCOMP
    W = IW
    THETA = ITHETA
    PHI = IPHI
    D = ID
    M = IM
    RI = IRI
    RO = IRO
    SET_DATA = TRUE
    !
    RETURN
  END SUBROUTINE SET_COIL_DATA
  !-------------------------------------------------------------------------------
END MODULE COIL_DATA
