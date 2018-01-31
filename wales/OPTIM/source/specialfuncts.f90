MODULE SPFUNCTS
! EFK 11/01/07
! This file contains subroutines for computing various special functions
!
  IMPLICIT NONE

  COMPLEX (KIND=KIND(1.0D0)), PARAMETER :: CMPLXI = (0.0D0,1.0D0)

CONTAINS
  SUBROUTINE INITIATE_RANDOM(SEED)
    ! initiate random number generator with the given seed

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: SEED
    INTEGER, ALLOCATABLE :: FULLSEED(:)
    INTEGER :: SEEDSIZE
    DOUBLE PRECISION :: RANDSTART(1000)

    CALL RANDOM_SEED(SIZE=SEEDSIZE)
    ALLOCATE(FULLSEED(SEEDSIZE))
    FULLSEED(:) = SEED

    CALL RANDOM_SEED(PUT=FULLSEED(:))
    CALL RANDOM_NUMBER(RANDSTART(1:1000))

    DEALLOCATE(FULLSEED)
    RETURN
  END SUBROUTINE INITIATE_RANDOM

  SUBROUTINE BAND2SCS(BANDM,VAL,INDX,JNDX,KD,N,NNZ)
    ! convert a band matrix of dimension 2*KD+1 where KD is the number
    ! of superdiagonals
    ! to a matrix in sparse coordinate storage format 
    ! (such as required by NAG sparse Cholesky decomposition routine F11JAF)
    ! the entries in the SCS matrix are ordered by increasing row and then
    ! by increasing column within each row
    ! only the upper triangle part of BANDM is used
    ! SCS is transposed (so that row numbers are always greater than columns)
    ! N by N is the dimension of the full matrix
    ! returns NNZ - the number of nonzero elements
    ! NNZ <= (KD+1)*N
    ! The VAL, INDX, JNDX arrays have dimension >= 2*NNZ

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, KD
    DOUBLE PRECISION, INTENT(IN) :: BANDM(2*KD+1,N)
    DOUBLE PRECISION, INTENT(OUT) :: VAL(2*(KD+1)*N)
    INTEGER, INTENT(OUT) :: INDX(2*(KD+1)*N), JNDX(2*(KD+1)*N), NNZ

    INTEGER :: I, R, C
    DOUBLE PRECISION :: V
    DOUBLE PRECISION, PARAMETER :: TINY=1.0D-10

    I = 0
    DO C = 1,KD !column
       DO R = 1,C !row
          V = BANDM(KD+1+R-C,C)
          IF (ABS(V).GT.TINY) THEN
             I = I + 1
             VAL(I) = V
             INDX(I) = C
             JNDX(I) = R
          ENDIF
       ENDDO
    ENDDO

    DO C = KD+1,N !column
       DO R = C-KD,C !row
          V = BANDM(KD+1+R-C,C)
          IF (ABS(V).GT.TINY) THEN
             I = I + 1
             VAL(I) = V
             INDX(I) = C
             JNDX(I) = R
          ENDIF
       ENDDO
    ENDDO

    NNZ = I
    
!    print*, 'Fraction nonsparse: ', DBLE(NNZ)/((KD+1)*N-KD*(KD+1)/2)
  END SUBROUTINE BAND2SCS

  SUBROUTINE DUMPCOORDS(X, FNAME, APPEND)
    USE COMMONS, ONLY : NATOMS
    IMPLICIT NONE
    CHARACTER (*) :: FNAME
    LOGICAL :: APPEND
    DOUBLE PRECISION :: X(3*NATOMS)
    INTEGER :: A

    ! given coordinate array X for a single molecule, dump into file FNAME
    IF (APPEND) THEN
       OPEN (UNIT = 55, FILE = FNAME, STATUS = 'UNKNOWN', POSITION = 'APPEND')
    ELSE
       OPEN (UNIT = 55, FILE = FNAME, STATUS = 'UNKNOWN')
    ENDIF
    WRITE(55,'(I6)') NATOMS
    WRITE(55,'(A)') ' '
    WRITE(55,'(A3,3G25.15)') ('AX ',X(3*(A-1)+1:3*(A-1)+3),A=1,NATOMS)
    CLOSE(55)
  END SUBROUTINE DUMPCOORDS

  SUBROUTINE DUMPFRAMES(X, NFRAME,FNAME)
    ! dump out the frames from a 1D coordinate array
    ! NFRAME is the number of frames
    USE COMMONS, ONLY : NATOMS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NFRAME
    CHARACTER (*), INTENT(IN) :: FNAME
    DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS*NFRAME)
    INTEGER :: A, F

    ! given coordinate array X for a single molecule, dump into file FNAME
    OPEN (UNIT = 55, FILE = FNAME, STATUS = 'UNKNOWN')
    DO F = 1,NFRAME
       WRITE(55,'(I6)') NATOMS
       WRITE(55,'(A,I4)') ' frame ', F
       WRITE(55,'(A3,3G25.15)') ('AX ',X(3*NATOMS*(F-1)+3*(A-1)+1:3*NATOMS*(F-1)+3*(A-1)+3),A=1,NATOMS)
    ENDDO
    CLOSE(55)

  END SUBROUTINE DUMPFRAMES

    SUBROUTINE CROSS_PRODUCT(A, B, C)
      ! take the cross product of 3D vectors A and B; return result in C

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: A(3), B(3)
      DOUBLE PRECISION, INTENT(OUT) :: C(3)

      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1)-A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)

      RETURN
    END SUBROUTINE CROSS_PRODUCT

    SUBROUTINE NORMALIZE(X)
      ! normalize a 3 dimensional vector

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: X(3)
      DOUBLE PRECISION :: DX

      DX = SQRT(DOT_PRODUCT(X,X))
      X(:) = X(:)/DX

      RETURN
    END SUBROUTINE NORMALIZE

    SUBROUTINE IRANDOMSELECT(CHOICES, PICKED, NC, NS, SEED)
      ! From an array CHOICES of size NC, select NP elements randomly
      ! and put results in PICKED
      ! CHOICES and PICKED should be integer arrays
      ! uses simple algorithm S from Vitter, 1984, Commun. ACM
      ! at each step, if m elements remain to be selected from N possibilities
      ! the probability of selecting the next element is m/N

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NC, NS, CHOICES(NC), SEED
      INTEGER, INTENT(OUT) :: PICKED(NS)

      INTEGER :: SEEDSIZE, N, S
      INTEGER, ALLOCATABLE :: FULLSEED(:)
      DOUBLE PRECISION :: RANSTART(1000), R

      IF (NS > NC) THEN
         print*, 'ERROR: NS > NC in IRANDOMSELECT'
         STOP
      END IF

      CALL RANDOM_SEED(SIZE=SEEDSIZE)
      ALLOCATE(FULLSEED(SEEDSIZE))
      FULLSEED(:) = SEED
      CALL RANDOM_SEED(PUT=FULLSEED(:SEEDSIZE))
      CALL RANDOM_NUMBER(RANSTART(1:1000)) ! get the generator started

      N = NC; S = NS
      DO 
         CALL RANDOM_NUMBER(R)
         IF (R*N <= S) THEN
            PICKED(NS-S+1) = CHOICES(NC-N+1)
            S = S - 1
         END IF
         N = N - 1
         IF (S <= 0) THEN
            EXIT
         ENDIF
      ENDDO
      DEALLOCATE(FULLSEED)
      RETURN      
    END SUBROUTINE IRANDOMSELECT

    SUBROUTINE SHIFTZERO(COORDS,NATMS,ATM)
      ! Shift entire molecule so that atom ATM is at the origin
      ! NATMS is the total number of atoms

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATMS, ATM
      DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATMS)
      DOUBLE PRECISION :: SHIFTS(3)
      INTEGER :: A, C

      SHIFTS(1:3) = COORDS(3*(ATM-1)+1:3*(ATM-1)+3)
      DO A = 1,NATMS
         DO C = 1,3
            COORDS(3*(A-1)+C) = COORDS(3*(A-1)+C) - SHIFTS(C)
         ENDDO
      ENDDO
    END SUBROUTINE SHIFTZERO

    SUBROUTINE ROTATEZERO(COORDS, NATMS, ATM, RAXIS, ZAXIS)
      ! Rotate the molecule around the RAXIS of 
      ! atom number ATM so as to zero the ZAXIS dimension of ATM
      ! RAXIS and ZAXIS must be 1,2, or 3 and not equal
      ! NATMS is the total number of atoms
      ! COORDP is a pointer to an array of dimension 3*NATMS

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATMS, ATM, RAXIS, ZAXIS
      DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATMS)
      INTEGER :: XAXIS, A
      DOUBLE PRECISION :: CT, ST, THETA, XX, YY
      COMPLEX (KIND=KIND(1.0D0)) :: DUMMYC, DUMMY1, DUMMY2

      IF (RAXIS.EQ.ZAXIS.OR.(RAXIS.LT.1.OR.RAXIS.GT.3) &
           & .OR.(ZAXIS.LT.1.OR.ZAXIS.GT.3)) THEN
         print*, 'RAXIS and ZAXIS must be btwn 1 and 3 and &
              & cannot be the same in ROTATETOZ!', RAXIS, ZAXIS
         STOP
      ENDIF

      DO XAXIS = 1,3
         IF (XAXIS.NE.RAXIS.AND.XAXIS.NE.ZAXIS) EXIT
      ENDDO
      DUMMY1 = COORDS(3*(ATM-1)+XAXIS)
      DUMMY2 = COORDS(3*(ATM-1)+ZAXIS)
      DUMMYC = LOG(DUMMY1 + DUMMY2*CMPLXI)

      THETA = -AIMAG(DUMMYC)
      CT = COS(THETA); ST = SIN(THETA)

      DO A = 1,NATMS
         XX = COORDS(3*(A-1)+XAXIS); YY = COORDS(3*(A-1)+ZAXIS)
         COORDS(3*(A-1)+XAXIS) = CT*XX-ST*YY
         COORDS(3*(A-1)+ZAXIS) = ST*XX+CT*YY
      ENDDO
    END SUBROUTINE ROTATEZERO

    SUBROUTINE ROTATE(COORDS, THETA, AXIS)
      ! Rotate a given point in 3D space around a given axis (1,2,or 3)
      ! By an angle theta (in radians)
      ! Return result as outcoords

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3)
      DOUBLE PRECISION, INTENT(IN) :: THETA
      INTEGER, INTENT(IN) :: AXIS
      DOUBLE PRECISION :: CT, ST, X, Y, Z

      CT = DCOS(THETA)
      ST = DSIN(THETA)

      X = COORDS(1); Y = COORDS(2); Z = COORDS(3)
      IF (AXIS.EQ.1) THEN
         COORDS(2) = CT*Y-ST*Z
         COORDS(3) = ST*Y+CT*Z
      ELSE IF (AXIS.EQ.2) THEN
         COORDS(1) = CT*X-ST*Z
         COORDS(3) = ST*X+CT*Z
      ELSE IF (AXIS.EQ.3) THEN
         COORDS(1) = CT*X-ST*Y
         COORDS(2) = ST*X+CT*Y
      ELSE
         print*, 'ERROR: in ROTATE, AXIS must be 1,2,or3', AXIS
         STOP
      ENDIF
      RETURN
    END SUBROUTINE ROTATE

    SUBROUTINE ROTATEANGLAXIS(X,THETA,V)
      ! rotate 3d vector X by an angle theta around vector V
      ! use Rodrigues' rotation formula

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: V(3)
      DOUBLE PRECISION, INTENT(IN) :: THETA
      DOUBLE PRECISION, INTENT(INOUT) :: X(3)
      DOUBLE PRECISION :: DUMMY1(3), DUMMY2

      CALL NORMALIZE(V)

      CALL CROSS_PRODUCT(V,X,DUMMY1)
      DUMMY2 = DOT_PRODUCT(V,X)*(1-COS(THETA))

      X(:) = X(:)*COS(THETA) + DUMMY1(:)*SIN(THETA) + DUMMY2*V(:)
      
      RETURN
    END SUBROUTINE ROTATEANGLAXIS

    SUBROUTINE LEGZO(N,X,W)
      ! Copied from: http://jin.ece.uiuc.edu/, mlegzo.for
      ! Written by Jian-Ming Jin
      !
      !       =========================================================
      !       Purpose : Compute the zeros of Legendre polynomial Pn(x)
      !                 in the interval [-1,1], and the corresponding
      !                 weighting coefficients for Gauss-Legendre
      !                 integration
      !       Input :   n    --- Order of the Legendre polynomial
      !       Output:   X(n) --- Zeros of the Legendre polynomial
      !                 W(n) --- Corresponding weighting coefficients
      !       =========================================================
      !
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: N, N0,NR,I,J,K
      DOUBLE PRECISION ::  X(N),W(N)
      N0=(N+1)/2
      DO 45 NR=1,N0
         Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10       Z0=Z
         P=1.0D0
         DO 15 I=1,NR-1
15          P=P*(Z-X(I))
            F0=1.0D0
            IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
            F1=Z
            DO 20 K=2,N
               PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
               PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
               F0=F1
20             F1=PF
               IF (Z.EQ.0.0) GO TO 40
               FD=PF/P
               Q=0.0D0
               DO 35 I=1,NR
                  WP=1.0D0
                  DO 30 J=1,NR
                     IF (J.NE.I) WP=WP*(Z-X(J))
30                   CONTINUE
35                   Q=Q+WP
                     GD=(PD-Q*FD)/P
                     Z=Z-FD/GD
                     IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40                   X(NR)=Z
                     X(N+1-NR)=-Z
                     W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
45                   W(N+1-NR)=W(NR)
                     RETURN
                   END SUBROUTINE LEGZO

 END MODULE SPFUNCTS

