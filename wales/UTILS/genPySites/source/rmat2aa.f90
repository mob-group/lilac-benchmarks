      SUBROUTINE RMAT2AA(RMAT,P)

!     transforms a rotation matrix RMAT(1:3,1:3) to its equivalent angle axis variables P(1:3)

      IMPLICIT NONE

      DOUBLE PRECISION RMAT(1:3,1:3), P(1:3)
      DOUBLE PRECISION :: ANGLE, FCTR

      ANGLE = ACOS(0.5D0*(RMAT(1,1) + RMAT(2,2) + RMAT(3,3) - 1.D0))
      FCTR  = 1.D0/(SQRT((RMAT(3,2) - RMAT(2,3))**2 + (RMAT(1,3) - RMAT(3,1))**2 + (RMAT(2,1) - RMAT(1,2))**2))
      P(1)  = FCTR*(RMAT(3,2) - RMAT(2,3))
      P(2)  = FCTR*(RMAT(1,3) - RMAT(3,1))
      P(3)  = FCTR*(RMAT(2,1) - RMAT(1,2))
      P(:)  = ANGLE*P(:)
      
!check for conditions which give FCTR =1/0 (i.e. symmetric matrix - zero rotation)
      IF (abs(RMAT(3,2)-RMAT(2,3)).lt.0.000000001) THEN
        IF (abs(RMAT(1,3)-RMAT(3,1)).lt.0.000000001)  THEN
          IF (abs(RMAT(2,1)-RMAT(1,2)).lt.0.000000001) THEN
          
          P(1)=0
          P(2)=0
          P(3)=2*3.1415926
	 END IF
	END IF
      END IF

      return
      END SUBROUTINE RMAT2AA
