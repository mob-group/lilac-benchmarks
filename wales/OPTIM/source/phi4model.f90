SUBROUTINE Phi4Model (XC, DX, ENERGY, GTEST,STEST)

    ! *****************************************************************************************


    USE COMMONS, ONLY: NOPT
    USE KEY, ONLY: JPARAM
    USE MODHESS

    IMPLICIT NONE
    DOUBLE PRECISION :: XC(NOPT), DX(NOPT) !Cartesian coords, derivatives
    DOUBLE PRECISION :: ENERGY, ENERGYNN, ENERGYDUMMY1, ENERGYDUMMY2, DUMMYDX
    LOGICAL          :: GTEST , STEST
    INTEGER          :: i,j, N, J1, J2
    REAL, PARAMETER  :: LAMBDA = 0.6D0 ! LAMBDA = 3/5
    REAL, PARAMETER  :: MUSQUARE = 1.00D0
!   REAL, PARAMETER  :: JPARAM = 0.001D0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D Nearest-Neighbour Phi^ model

!    ENERGY = 0.00    
!    DO i = 1, NATOMS
!      ENERGY = ENERGY + (LAMBDA/(24.D0))* XC(i)**4.D0 - (MUSQUARE/(2.D0))*XC(i)**2.D0

!       ENERGYNN = 0.00
!         IF (i==1) THEN
!            ENERGYNN = ENERGYNN + (JPARAM/(4.D0))*(XC(i) - XC(NATOMS))**2.D0 + (JPARAM/(4.D0))*(XC(i) - XC(i+1))**2.D0
!             ELSE IF (i==NATOMS) THEN
!                ENERGYNN = ENERGYNN + (JPARAM/(4.D0))*(XC(i) - XC(i-1))**2.D0 + (JPARAM/(4.D0))*(XC(i) - XC(1))**2.D0
!             ELSE
!          ENERGYNN = ENERGYNN + (JPARAM/(4.D0))*(XC(i) - XC(i-1))**2.D0 + (JPARAM/(4.D0))*(XC(i) - XC(i+1))**2.D0
!          END IF
!       ENERGY = ENERGY + ENERGYNN
!    ENDDO
       


!    IF (GTEST) THEN
!       DO i = 1, NATOMS
!          DX(i) = (LAMBDA/6.D0)* XC(i)**3.D0 + ((4.D0) * JPARAM - MUSQUARE)*XC(i)
!          DUMMYDX = 0.00
!          IF (i==1) THEN
!             DUMMYDX = DUMMYDX + XC(NATOMS) + XC(i+1)
!             ELSE IF (i==NATOMS) THEN
!                DUMMYDX = DUMMYDX + XC(i-1) + XC(1)
!                ELSE
!                   DUMMYDX = DUMMYDX + XC(i-1) + XC(i+1)
!                END IF
!                DUMMYDX = JPARAM*DUMMYDX

!                DX(i) = DX(i) - DUMMYDX
!             ENDDO
    
!          ELSE


!          ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mean-field phi^4 model

N=NOPT

 ENERGY=0.0D0
 DO i = 1, N
    ENERGY = ENERGY + (-0.5D0*XC(i)**2 + 0.25D0*XC(i)**4)
 ENDDO
 ENERGYDUMMY1 = 0.0D0
 DO i = 1, N
    ENERGYDUMMY1 = ENERGYDUMMY1 + XC(i)
 ENDDO
 ENERGY = ENERGY - (JPARAM/(2.0D0*N)) *(ENERGYDUMMY1**2)
 
 IF (.NOT. GTEST) RETURN
 
 DUMMYDX = 0.0D0
 DO i = 1, N
    DUMMYDX=DUMMYDX+XC(i)
 ENDDO
 DUMMYDX=DUMMYDX*(JPARAM/(1.0D0*N))
 
!
! gradient
! 
 DO I=1,N
    DX(I)=-XC(I)+XC(I)**3-DUMMYDX
 ENDDO

!
! Hessian. DJW
!
 IF (STEST) THEN
    DO J1=1,N
       HESS(J1,J1)=-1.0D0+3.0D0*XC(J1)**2-(JPARAM/N)*(XC(J1)**2+DUMMYDX)
       DO J2=J1+1,N
          HESS(J1,J2)=-JPARAM*XC(J1)*XC(J2)/(1.0D0*N)
          HESS(J2,J1)=HESS(J1,J2)
       ENDDO
    ENDDO
!   PRINT '(A)','Coordinates:'
!   PRINT '(3G20.10)',XC(1:NOPT)
!   PRINT '(A)','Hessian elements'
!   PRINT '(3G20.10)',HESS(1:NOPT,1:NOPT)
 ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simple example

!ENERGY = XC(i)**3 + XC(i)**2 + 2.0

!IF (.NOT. GTEST) RETURN 
!DO i = 1, NATOMS
!DX(i) = 3.0 * XC(i)**2 + 2.0*XC(i)
!ENDDO

!ENERGY = XC(1)**4 + 4.0*XC(1)**2 * XC(2)**2 - 2.0 *XC(1)**2 + 2.0 * XC(2)**2

!IF (.NOT. GTEST) RETURN
!DX(1) = 4.0 * XC(1)**3 + 8.0 * XC(1) * XC(2)**2 - 4.0*XC(1) 
!DX(2) = 8.0*XC(1)**2 * XC(2) + 4.0*XC(2)
!DX(3:9) = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   IF (STEST) THEN
!       PRINT*, 'Warning: There is no analytical hessian implemented for the phi4 model yet.'
!
! This call is not allowed because it implies recursive calls to 
! non-recursive subroutine potential. DJW
!
!       CALL MAKENUMHESS(XC,NATOMS)
!   ENDIF

!call writespec_xyz(17,XC)

END SUBROUTINE Phi4Model

SUBROUTINE Phi4dist(COORDSB,COORDSA,NATOMS,DISTANCE)
    USE COMMONS, ONLY: NOPT
    IMPLICIT NONE
    double precision, intent(in) :: COORDSB(NOPT),COORDSA(NOPT)
    double precision, intent(out) :: DISTANCE
    integer :: icount, NATOMS

    DO ICOUNT = 1,NOPT
        DISTANCE = DISTANCE + (COORDSA(ICOUNT) - COORDSB(ICOUNT))**2
    ENDDO
    DISTANCE = DSQRT(DISTANCE)

END SUBROUTINE Phi4dist
