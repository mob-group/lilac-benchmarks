SUBROUTINE example1Dpotential(XC, DX, ENERGY, GTEST,SECT)

    ! *****************************************************************************************
  USE MODHESS
  USE COMMONS, ONLY: NATOMS

    IMPLICIT NONE
    DOUBLE PRECISION :: XC(NATOMS), DX(NATOMS) !Cartesian coords, derivatives
    DOUBLE PRECISION :: ENERGY, ENERGYNN, ENERGYDUMMY1, ENERGYDUMMY2, DUMMYDX
    LOGICAL          :: GTEST , SECT
    INTEGER          :: i,j, N
    REAL, PARAMETER  :: LAMBDA = 0.6D0 ! LAMBDA = 3/5
    REAL, PARAMETER  :: MUSQUARE = 1.0D0
    REAL, PARAMETER  :: JPARAM = 1.0D0


    N=NATOMS
    
    ENERGY=0.0D0
    DO i = 1, N
       ENERGY = ENERGY + (-0.5D0*XC(i)**2 + 0.25D0*XC(i)**4)
    ENDDO
    ENERGYDUMMY1 = 0.0D0
    DO i = 1, N
       ENERGYDUMMY1 = ENERGYDUMMY1 + XC(i)
    ENDDO
    ENERGY = ENERGY - (JPARAM/(2.0D0 * N)) *(ENERGYDUMMY1**2)
    
    IF (.NOT.(GTEST.OR.SECT)) RETURN
    
    DUMMYDX = 0.0D0
    DO i = 1, N
       DUMMYDX = DUMMYDX + XC(i)
    ENDDO
    DUMMYDX = (JPARAM/N)*DUMMYDX

! Gradient of the mean-field phi^4 model
     
    DO i = 1, N
       DX(i) = - XC(i) + XC(i)**3 - DUMMYDX
    ENDDO

      
    
IF (.NOT.SECT) RETURN
! if there is no analytical Hessian available for your potential, then call MAKENUMHESS(XC, NATOMS/3) subroutine. Here, NATOMS/3 is due to the fact that the subroutine
! assumes that the potential is defined for 3 coordinates. For 1D potential, dividing the argument NATOMS by 3 in the subroutine makes the subroutine work fine for the
! 1D case. The limitation is that then one must have the no. of atoms in multiples of 3 only. 
! If the analytical hessian is known and coded, as below, then the above limitation doesn't apply and one can work with any NATOMS.
 
!       PRINT*, 'Warning: There is no analytical hessian implemented for the phi4 model yet.'
!       CALL MAKENUMHESS(XC,NATOMS/3)
       CALL Hessian_EX1D(XC)


!call writespec_xyz(17,XC)

END SUBROUTINE Example1Dpotential

!!!!!!!!!!!!!!!!!!!
! Hessian of the mean-field phi^4 model
!!!!!!!!!!!!!!!!!!!

SUBROUTINE Hessian_EX1D(XC)
  USE COMMONS, ONLY : NATOMS
  USE MODHESS
  IMPLICIT NONE
  INTEGER  :: i,j, N
  DOUBLE PRECISION, INTENT(IN), DIMENSION(NATOMS) :: XC
  REAL, PARAMETER  :: JPARAM = 1.0D0

  N=NATOMS
  
  !Initialise matrix to zeros
  HESS = 0.0D0

  ! Initialise the upper and diagonal components
  DO i=1,N

     HESS(i,i) = -1.0D0 + 3.0D0*XC(i)**2 - (JPARAM/N)

  ENDDO

  DO i=1,N
     DO j= i+1,N    
       HESS(i,j) = - (JPARAM/N)
       HESS(j,i) = HESS(i,j)
     ENDDO 
    ENDDO

!  HESS(:,:) = HESS(:,:)/n !Scaling factor

END SUBROUTINE Hessian_EX1D


! The below subroutine must only be called for the case when there is at least one trivial zero (i.e., continuous symmetry) mode in the system.
! To call this subroutine, one needs to put the following line in two files: efol.f and bfgsts.f90
! IF (VARIABLES.AND.EX1DT) CALL SHIFTEX1D


SUBROUTINE SHIFTEX1D
  USE KEY, ONLY : SHIFTV
  USE COMMONS, ONLY : NATOMS
  USE MODHESS
  IMPLICIT NONE
  INTEGER :: N, J1, J2

  N=NATOMS
  DO J1=1,N
     DO J2=1,N
        HESS(J2, J1) = HESS(J2,J1) + SHIFTV
     ENDDO
  ENDDO

END SUBROUTINE SHIFTEX1D
