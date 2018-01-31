!
! XY model with Periodic Boundary conditions. ch558
!

SUBROUTINE Energy_1d_PBC_INVT(THETA,GRAD,ENERGY,GTEST,SECT)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  IMPLICIT NONE
  INTEGER ::  N
  DOUBLE PRECISION , INTENT(IN), DIMENSION(NONEDAPBC) :: THETA
  DOUBLE PRECISION , DIMENSION(NONEDAPBC) :: GRAD
  DOUBLE PRECISION :: ENERGY
  LOGICAL GTEST, SECT
  
  N = NONEDAPBC

 ! theta(n)=0
  energy = sum( cos(xyphi(1:n-1) + theta(2:n) - theta(1:n-1))  )
  energy = energy + cos( xyphi(n) + theta(1)  - theta(n))
  energy = -1.0D0 + (energy/n)
  !energy = n - energy

  IF (.NOT.(GTEST.OR.SECT)) RETURN
  
  grad(1) = sin(xyphi(1) + theta(2) - theta(1)) - sin( xyphi(n) + theta(1) - theta(n))
  
  grad(2:(n-1)) = -sin( xyphi(1:(n-2)) + theta(2:(n-1)) - theta(1:(n-2))) + sin( xyphi(2:(n-1)) + theta(3:n) - theta(2:(n-1)))
  
  grad(n)= -sin(xyphi(n-1) + theta(n) - theta(n-1)) + sin(xyphi(n) + theta(1) -theta(n))

 ! grad(n)=0

  IF(.NOT.SECT) RETURN

  CALL  Hessian_1d_PBC_INVT(THETA)

END SUBROUTINE Energy_1d_PBC_INVT




SUBROUTINE Hessian_1d_PBC_INVT(THETA)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  USE MODHESS
  IMPLICIT NONE
  INTEGER  :: i,n,BC
  DOUBLE PRECISION, INTENT(IN), DIMENSION(NONEDAPBC) :: theta
  ! With PBC with fix one lattice site and require it to be zero. 
  ! This means grad_n=0 and hess_n=0
  
  N=NONEDAPBC
  !  theta(n)=0 ! Implement PBC
  
  !Initialise matrix to zeros
  HESS = 0.0D0

  ! Initialise the upper and diagonal components
  DO i=1,n

     HESS(i,i)   =  -cos(xyphi(i) + theta( BC(i+1,n) ) - theta(i) ) - cos( xyphi(BC(i-1,N) ) + theta(i) - theta( BC(i-1, n ) )    )

     HESS(BC(i+1,N) ,i) = +cos(xyphi(i) + theta( BC(i+1,N) ) - theta(i))
  ENDDO

  ! Since Hessian is symmetric, use upper com. for lower. 
  DO i=1,n
     HESS(BC(i-1,N) ,i) = HESS(i, BC(i-1,N) )
  ENDDO

!  HESS(:,:) = HESS(:,:)/n !Scaling factor

END SUBROUTINE Hessian_1d_PBC_INVT

