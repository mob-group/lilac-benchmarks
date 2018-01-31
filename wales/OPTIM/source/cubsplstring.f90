MODULE CUBSPLSTRING
  ! functions for dealing with a string interpolated via cubic splines
  ! WARNING: ARC LENGTH DERIVATIVES ARE NOT ALWAYS QUITE RIGHT!!!!
  ! THIS CAN OCCASSIONALLY MAKE FOLLOWARCNEWT FAIL TO CONVERGE

  USE GSDATA

  IMPLICIT NONE

  INTEGER :: DQAGKEY

  DOUBLE PRECISION, ALLOCATABLE, TARGET :: COEFF(:,:,:), ABSC(:)
  DOUBLE PRECISION, TARGET :: ARCTOL
  COMPLEX, PARAMETER :: CMPLXI=(0.0D0,1.0D0)

  INTEGER :: GLOBALB ! for use in calculating arclength with DQAG only

CONTAINS
  SUBROUTINE INTERPCUBSPL
    ! string with interpolate with cubic splines

    IMPLICIT NONE
    TYPE(IMGNODE), POINTER :: DUMMYP
    INTEGER, PARAMETER :: IBCBEG = 0, IBCEND = 0
    INTEGER :: CRD, IM

    DUMMYP => FIRST%PREV
    DO IM = 1, NIM+1 
       DUMMYP => DUMMYP%NEXT
       DUMMYP%DIFF(:) = DUMMYP%XYZ(:) - DUMMYP%PREV%XYZ(:)
       DUMMYP%CHORD = SQRT(DOT_PRODUCT(DUMMYP%DIFF, DUMMYP%DIFF))
       ABSC(IM+1) = ABSC(IM) + DUMMYP%CHORD ! abscissa
       IF(.NOT.CUBSPLT) DUMMYP%ARC=ABSC(IM+1)
       COEFF(1,IM+1,:) = DUMMYP%XYZ(:) ! cubic spline coefficients  
    ENDDO

    DO CRD = 1,NC
       ! set second derivatives at end points equal to zero
       COEFF(2,1,CRD) = 0.0D0; COEFF(2,NIM+2,CRD) = 0.0D0

       CALL CUBSPL(ABSC(1:NIM+2), COEFF(1:4,1:NIM+2,CRD), NIM+2, IBCBEG, IBCEND)       
       COEFF(3,1:NIM+2,CRD) = COEFF(3,1:NIM+2,CRD) / 2
       COEFF(4,1:NIM+2,CRD) = COEFF(4,1:NIM+2,CRD) / 6
    ENDDO    
    
    IF(CUBSPLT) THEN
    ! calculate the arc lengths at each image
    DUMMYP => FIRST%PREV
    DO IM = 1,NIM+1
       DUMMYP => DUMMYP%NEXT
       DUMMYP%ARCDIFF = ARCLENGTH(DUMMYP%CHORD, IM-1)
       DUMMYP%ARC = DUMMYP%PREV%ARC + DUMMYP%ARCDIFF
    ENDDO
    ENDIF
    STRINGLEN = DUMMYP%ARC

  END SUBROUTINE INTERPCUBSPL

  SUBROUTINE GETSPLVAL(X,B,N,ANS)
    ! Get value of Nth derivative of spline, at abscissa X over B, put into ANS
    ! WARNING: does NOT check for X to be within limits

    DOUBLE PRECISION, INTENT(IN) :: X
    INTEGER, INTENT(IN) :: B, N
    DOUBLE PRECISION, INTENT(OUT) :: ANS(NC)

    IF (N.EQ.0) THEN
       ANS(:) = COEFF(4,B+1,:)*X**3+COEFF(3,B+1,:)*X**2+COEFF(2,B+1,:)*X + COEFF(1,B+1,:)
    ELSE IF(N.EQ.1) THEN
       ANS(:) = 3*COEFF(4,B+1,:)*X**2+2*COEFF(3,B+1,:)*X+COEFF(2,B+1,:)
    ELSE IF(N.EQ.2) THEN
       ANS(:) = 6*COEFF(4,B+1,:)*X+2*COEFF(3,B+1,:)
    ELSE
       print*, 'N must be 0,1,or 2', N
       STOP
    ENDIF

  END SUBROUTINE GETSPLVAL

  DOUBLE PRECISION FUNCTION ARCDERV(X,B)
    ! get norm of tangent at spline parameter X over B
    IMPLICIT NONE

    DOUBLE PRECISION :: X, D(NC)
    INTEGER :: B

    CALL GETSPLVAL(X,B,1,D) ! get spline derivative
    ARCDERV = SQRT(DOT_PRODUCT(D,D))
!    IF(BMCURPAIR.EQ.10) print*, 'TESTXA: ', X, B, ARCDERV
    IF (ARCDERV.LT.TINY) THEN
       print*, 'ERROR: very small arcderv - bad parametrization', X, B, ARCDERV
       STOP
    ENDIF
  END FUNCTION ARCDERV  

  DOUBLE PRECISION FUNCTION ARCDERVF(X)
    ! function for calling with DQAG integration
    ! works just like function ARCDERV but uses the 
    ! variable GLOBALB instead of B
    ! get norm of tangent at spline parameter X over B
    IMPLICIT NONE

    DOUBLE PRECISION :: X
    DOUBLE PRECISION :: D(NC)

    CALL GETSPLVAL(X,GLOBALB,1,D) ! get spline derivative
    ARCDERVF = SQRT(DOT_PRODUCT(D,D))
    IF (ARCDERVF.LT.TINY) THEN
       print*, 'ERROR: very small arcderv - bad parametrization', &
            & X, GLOBALB, ARCDERVF
       STOP
    ENDIF
  END FUNCTION ARCDERVF

  DOUBLE PRECISION FUNCTION ARCLENGTH(X,B)
    ! calculate the arc length from base image B to X above that image
    ! assume X is between 0 and absc(b+2)-absc(b+1)
    ! uses adaptive quadrature function from SLATEC (DQAG)

    IMPLICIT NONE
    
    DOUBLE PRECISION :: X
    INTEGER :: B
    
    INTEGER, PARAMETER :: LIMIT = 100 ! maximal number of subdivisions
    DOUBLE PRECISION :: ANS, ERR
    INTEGER :: NEVAL, IERR
    INTEGER :: LAST, IWORK(LIMIT)
    DOUBLE PRECISION :: WORK(LIMIT*4)

    GLOBALB = B
    CALL DQAG(ARCDERVF,0.0D0, X, ARCTOL*0.1, 1.0D-16, DQAGKEY, ANS, ERR, NEVAL, IERR, LIMIT, LIMIT*4, LAST, IWORK, WORK)
    
    ARCLENGTH = ANS

!    print*, 'TESTXA:', ANS, ERR, LAST
    IF(IERR.NE.0) THEN
       print*, 'Error in calculating arc length: ', X, B, ANS, ERR
       STOP
    ENDIF
    RETURN

  END FUNCTION ARCLENGTH    

  DOUBLE PRECISION FUNCTION FOLLOWARCNEWT(A, S, R)
    ! S is the image providing the starting point
    ! R is location of S relative to A: if R = 0, S is lower bound, if R=1, S is upper bound
    ! Use Newton-Raphston to find the parameter value within appropriate spline segment
    ! that gives an arc length of A within that segment

    IMPLICIT NONE

    DOUBLE PRECISION :: A
    INTEGER :: S, R

    INTEGER :: ITS, MAXITS, B
    DOUBLE PRECISION :: X, FX, FPX, STEP
    DOUBLE PRECISION :: LOWER, UPPER

    MAXITS = 20

    IF (R.EQ.0) THEN
       B = S
    ELSE IF (R.EQ.1) THEN
       B = S-1
    ELSE
       print*, 'R must be 0 or 1.'
       STOP
    ENDIF

    UPPER = ABSC(B+2)-ABSC(B+1)
    LOWER = 0.0D0

!    print*, 'TESTFOLLOWARcNEWT: ', A, S, B, R, NIM, TIM
    
    X = ABSC(S+1)-ABSC(B+1)
    STEP = ARCTOL + 1
    ITS = 0
    DO WHILE (STEP.GT.ARCTOL)
       ITS = ITS + 1
       IF (ITS.GT.MAXITS) THEN
          IF (PTEST) print*, 'FOLLOWARC NEWT FAILED! Taking bisection step', A, S, R, ITS          
          X = (UPPER+LOWER)/2
          ITS = 1          
       ENDIF       
       FX = ARCLENGTH(X,B) - A

       ! keeps track of a lower and upper bound in case have to call bisection method later

       IF (FX.GT.0.AND.X.LT.UPPER) THEN
          UPPER = X
       ELSE IF (FX.LT.0.AND.X.GT.LOWER) THEN
          LOWER = X
       ENDIF

       FPX = ARCDERV(X,B)
       STEP = ABS(FX/FPX)
       X = X-FX/FPX
    ENDDO

    FOLLOWARCNEWT = X
        
    RETURN
  END FUNCTION FOLLOWARCNEWT

  DOUBLE PRECISION FUNCTION FOLLOWARCBIS(A, B, L, U)
    ! Find parameter corresponding to given arclength by bisection
    ! less efficient than FOLLOWARCNEWT, only for use where that doesn't converge
    ! Returns parameter above base image B corresponding to arclength A above B
    ! L and U are starting upper and lower bounds

    DOUBLE PRECISION :: A, L, U
    INTEGER :: B, ITS
    DOUBLE PRECISION :: LOWER, UPPER, X, AX
    
    LOWER = L
    UPPER = U
    ITS = 0
    DO WHILE (ABS(AX-A).GT.ARCTOL.OR.ITS.EQ.0)
       ITS = ITS+1
       X = (UPPER+LOWER)/2
       AX = ARCLENGTH(X,B)
       print '(A,5F20.10)', 'TESTXB: ', X, AX, LOWER, UPPER, A
       IF (AX.GT.A) THEN
          UPPER = X
       ELSE
          LOWER = X
       ENDIF
       IF (ITS.GT.1.0D6) THEN
          print*, 'ERROR: BISECTION METHOD NOT CONVERGING! Something is wrong.'
          STOP
       ENDIF
    ENDDO
    FOLLOWARCBIS = X
    RETURN
  END FUNCTION FOLLOWARCBIS

END MODULE CUBSPLSTRING

