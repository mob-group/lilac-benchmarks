      SUBROUTINE AA2RMAT(P,RMAT)

!     transforms angle axis to a rotation matrix RMAT(1:3,1:3) 

      IMPLICIT NONE

      DOUBLE PRECISION RMAT(1:3,1:3), P(1:3)
      DOUBLE PRECISION :: THETA, K(1:3), KX(1:3,1:3), OUTERK(1:3,1:3)
      INTEGER :: i,j

      !compute size of p
      THETA =SQRT(P(1)**2+P(2)**2+P(3)**2)

      !set up RMAT as identity
      RMAT=0
      RMAT(1,1)=1.0D0
      RMAT(2,2)=1.0D0
      RMAT(3,3)=1.0D0
      IF (abs(THETA).GE.1D-9) THEN
          K=P/THETA
          
         
          KX(1,1)=0
          KX(1,2)=-K(3)
          KX(1,3)=K(2)
          KX(2,1)=K(3)
          KX(2,2)=0
          KX(2,3)=-K(1)
          KX(3,1)=-K(2)
          KX(3,2)=K(1)
          KX(3,3)=0
          RMAT=RMAT*COS(THETA)
          RMAT=RMAT+KX*sin(THETA)
          DO i=1,3
             DO j=1,3
                 RMAT(i,j)=RMAT(i,j)+(1-COS(THETA))*K(i)*K(j)
             END DO
          END DO
      END IF
      return
      END SUBROUTINE AA2RMAT
