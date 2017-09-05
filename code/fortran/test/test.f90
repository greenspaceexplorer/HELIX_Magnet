PROGRAM TEST
  IMPLICIT NONE
  REAL*8 :: A(3,3)
  REAL*8 :: SQUARE(SHAPE(A))

  A = RESHAPE((/1,2,3,4,5,6,7,8,9/),(/3,3/))
  A = SQUARE(A)
  PRINT *,A
  PRINT *,SHAPE(A)
  
  STOP
END PROGRAM TEST

FUNCTION SQUARE(A)
  IMPLICIT NONE
  INTEGER :: K,L,S(2)
  REAL*8 :: A(:,:)
  REAL*8 :: SQUARE(SHAPE(A))
  S = SIZE(A)
  DO K=1,S(1)
     DO L=1,S(2)
        !        SQUARE(K,L) = A(K,L)*A(K,L)
        PRINT *,A
        PRINT *,S
     ENDDO
  ENDDO
  
  RETURN
END FUNCTION SQUARE