! 
! One-dimensional anti-periodic XY model. ch558
!
SUBROUTINE Energy_1d_APBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  IMPLICIT NONE
  INTEGER N
  DOUBLE PRECISION, dimension(NONEDAPBC) :: theta, GRAD
  DOUBLE PRECISION :: Energy
  LOGICAL GTEST,SECT

  n=NONEDAPBC

  ENERGY = 0.0D0
 Energy = sum(cos(xyphi(1:n-1)+theta(2:n)-theta(1:n-1)))  
 Energy = Energy + cos(xyphi(n) - theta(1)-theta(n))
 Energy = 1 - (Energy/n)
 !Energy = n - Energy

 IF (.NOT.(GTEST.OR.SECT)) RETURN

 grad(1)=-sin(xyphi(1) + theta(2) - theta(1)) - sin( xyphi(n) - theta(1) - theta(n))

 grad(n)= sin(xyphi(n-1) + theta(n) - theta(n-1)) - sin(xyphi(n) - theta(1) -theta(n))

 grad(2:(n-1))= sin(xyphi(1:(n-2)) + theta(2:(n-1)) - theta(1:(n-2))) - sin( xyphi(2:(n-1)) + theta(3:n) - theta(2:n-1))

 IF (.NOT.SECT) RETURN

 CALL Hessian_1d_APBC(THETA)

END SUBROUTINE ENERGY_1d_APBC

SUBROUTINE Hessian_1d_APBC(THETA)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  USE MODHESS
  IMPLICIT NONE
  INTEGER :: n,i, BC
  DOUBLE PRECISION :: APBC
  DOUBLE PRECISION, INTENT(IN), DIMENSION(1:NONEDAPBC) :: theta
  
  N=NONEDAPBC
  
  !Initialise matrix to zeros
  HESS = 0.0D0
  
  ! Implement boundary conditions with BC function. We also need to use
  ! APBC(i,N) for minus one in theta[i+N] = -theta[i]
  ! Initialise the upper and diagonal components

  DO i=1,n
     
     HESS(i,i)   =  cos( xyphi(i) + APBC(i+1,N)*theta( BC(i+1,n) ) - theta(i) ) + cos( xyphi(BC(i-1,N) ) + APBC(i-1,N)*theta(i) - theta( BC(i-1, n ) )    )
     
     HESS(BC(i+1,N) ,i) = -1.0D0*APBC(i+1,N)*cos(xyphi(i) + APBC(i+1,N)*theta( BC(i+1,N) ) - theta(i))
  ENDDO
  
  ! Since Hessian is symmetric, use upper com. for lower. 
  DO i=1,n
     HESS(BC(i-1,N) ,i) = HESS(i, BC(i-1,N) )
  ENDDO
  
 
  !  HESS(:,:) = HESS(:,:)/n ! Need to Scale, leave out and scale later
  
END SUBROUTINE Hessian_1d_APBC

!APBC is a function to get minus sign for boundary conditions if needed
! theta[i+N] = -theta[i]

FUNCTION APBC(i,N) 
  IMPLICIT NONE
  DOUBLE PRECISION :: APBC
  INTEGER, INTENT(IN) :: i, N
  
  if ( (i.EQ.0).OR.(i.EQ.(N+1)) ) THEN
     APBC = -1
  else
     APBC = 1
  ENDIF
END FUNCTION APBC
