MODULE ORBITALS_MOD
IMPLICIT NONE

  CONTAINS

    SUBROUTINE ORBITALS_INIT(NORBS, NROTS)

        ! Initialise information needed for orbital optimisation.
        ! This requires reading in the integral arrays from file.

        USE COMMONS, ONLY: R2INTS, DIPINTS, ORBVAREXPONENT

        INTEGER, INTENT(OUT) :: NORBS, NROTS

        INTEGER :: ORB_FILE_UNIT, GETUNIT

        CHARACTER (LEN=10) :: WORD

        LOGICAL :: END

        ! First check we've initialised sensibly.
        IF (ORBVAREXPONENT.LT.1) THEN
            STOP 'Need to set penalty function exponent to a positive value.'
        END IF

        NORBS = -1

        CALL CHECK_FILE_EXISTS('INTEGRALS_INFO', EX=.TRUE.)

        ORB_FILE_UNIT=GETUNIT()
        OPEN(UNIT=ORB_FILE_UNIT,FILE='INTEGRALS_INFO',STATUS='OLD')
        READ(ORB_FILE_UNIT, *) NORBS

        NROTS = NORBS * (NORBS - 1) / 2
        !PRINT '(A,3I10)','NORBS,NROTS,ORBVAREXPONENT=',NORBS,NROTS,ORBVAREXPONENT

        IF (NORBS.LT.1) THEN
            STOP 'Could not successfully read the system information from file.'
        ELSE IF (NROTS.LT.1) THEN
            STOP 'We appear to have no possible orbital rotations.'
        END IF

        CLOSE(ORB_FILE_UNIT)

        ALLOCATE(R2INTS(NORBS,NORBS))
        ALLOCATE(DIPINTS(3,NORBS,NORBS))

        !PRINT '(A)','R2'
        CALL READ_INTEGRALS('INTEGRALS_R2.csv', NORBS, R2INTS(:,:))
        !PRINT '(A)','X'
        CALL READ_INTEGRALS('INTEGRALS_X.csv', NORBS, DIPINTS(1,:,:))
        !PRINT '(A)','Y'
        CALL READ_INTEGRALS('INTEGRALS_Y.csv', NORBS, DIPINTS(2,:,:))
        !PRINT '(A)','Z'
        CALL READ_INTEGRALS('INTEGRALS_Z.csv', NORBS, DIPINTS(3,:,:))

    END SUBROUTINE ORBITALS_INIT

    SUBROUTINE READ_INTEGRALS(INT_FILENAME, NORBS, INTS)

        ! Given filename of intergral files reads the integrals within
        ! it into the array INTS. Integrals are assumed to always have
        ! dimension NORBSxNORBS.

        CHARACTER (LEN=*), INTENT(IN) :: INT_FILENAME
        INTEGER, INTENT(IN) :: NORBS
        DOUBLE PRECISION, INTENT(OUT) :: INTS(NORBS, NORBS)
        CHARACTER (LEN=32) :: CALCID
        INTEGER :: INT_FILE_UNIT, I, GETUNIT
        LOGICAL :: END

        CALL CHECK_FILE_EXISTS(INT_FILENAME, EX=.TRUE.)
        INT_FILE_UNIT=GETUNIT()
        OPEN(UNIT=INT_FILE_UNIT,FILE=INT_FILENAME,STATUS='OLD')
!        CALL FILE_OPEN(INT_FILENAME,INT_FILE_UNIT,.FALSE.)
!        CALL INPUT(END, INT_FILE_UNIT)
        READ(INT_FILE_UNIT, '(A)') CALCID
        !CALL FILE_OPEN(INT_FILENAME,INT_FILE_UNIT,.FALSE.)
        !CALL INPUT(END, INT_FILE_UNIT)
        !CALL READU(CALCID)
        !PRINT '(A,I10,X,A)','INT_FILE_UNIT,CALCID=',INT_FILE_UNIT,CALCID
        !PRINT '(A,I10,A)','NORBS=',NORBS

        DO I = 1, NORBS
           READ(INT_FILE_UNIT,*) INTS(I,1:NORBS) 
           !PRINT *,INTS(I,1:NORBS)
        END DO

        CLOSE(INT_FILE_UNIT)
    END SUBROUTINE READ_INTEGRALS

    SUBROUTINE GET_ORBITAL_LOCALITY(X, GRAD, LOCALITY, GTEST, SECT)

        ! Obtains the Power of the Orbital Variance locality measure for
        ! a given position. If GTEST is set also returns the gradient.

        ! Note that for intermediates dependent upon coefficient of current
        ! (rotated) orbital set in the original MO basis we use indexing
        ! of (ORIG_ORB,NEW_ORB) in all cases.

        USE COMMONS, ONLY: NROTS, NORBS, R2INTS, DIPINTS, ORBVAREXPONENT

        DOUBLE PRECISION, INTENT(INOUT) :: X(NROTS)
        DOUBLE PRECISION, INTENT(OUT) :: GRAD(NROTS), LOCALITY
        LOGICAL, INTENT(IN) :: GTEST,SECT

        DOUBLE PRECISION :: ROTATION(NORBS,NORBS), ORBVAR(NORBS)
        DOUBLE PRECISION :: ORBDIPOLES(3,NORBS), IND_ORBVAR_DERIV(NORBS,NORBS)
        DOUBLE PRECISION :: PENALTY_DERIV(NORBS,NORBS), ROT_DERIV(NORBS,NORBS)
        DOUBLE PRECISION :: PROD_DERIVS(NORBS,NORBS)

        INTEGER :: I,J,K

        CALL CHECK_COORDINATES(X)

        CALL GET_ROTATION(X, ROTATION)

        CALL GET_CURRENT_ORB_DIPOLES(ROTATION, ORBDIPOLES)
        DO I = 1, NORBS
            ORBVAR(I) = -SUM(ORBDIPOLES(1:3,I)**2)
            !PRINT*,ORBVAR(I)
            DO J = 1, NORBS
                DO K = 1, NORBS
                    ORBVAR(I)=ORBVAR(I)+R2INTS(J,K)*ROTATION(J,I)*ROTATION(K,I)
                END DO
            END DO
            !PRINT *,'ORBITAL',I,',ORBVAR=',ORBVAR(I),',DIPOLES=',ORBDIPOLES(1:3,I)
        END DO

        !PRINT *,'ORBVARS=',ORBVAR
        !PRINT *,'ORBVARSEXP=',ORBVAR**ORBVAREXPONENT

        LOCALITY = SUM(ORBVAR ** ORBVAREXPONENT)

        !PRINT *,'LOCALITY=',LOCALITY

        IF (GTEST) THEN
            !PRINT *,'ORBDIPOLES='
            !DO J = 1,3
            !    PRINT *,ORBDIPOLES(J,1:NORBS)
            !END DO

            CALL GET_INDIVIDUAL_DERIVATIVES(ROTATION, ORBDIPOLES, IND_ORBVAR_DERIV)
            !PRINT *,'IND_ORBVAR_DERIV='
            !DO J = 1,NORBS
            !    PRINT *, IND_ORBVAR_DERIV(J,1:NORBS)
            !END DO

            CALL GET_PENALTY_DERIVATIVES(IND_ORBVAR_DERIV, ORBVAR, PENALTY_DERIV)

            !PRINT *,'PEN_DERIV='
            !DO J = 1,NORBS
            !    PRINT *, PENALTY_DERIV(J,1:NORBS)
            !END DO

            DO I = 1, NROTS
                CALL GET_ROTATION_DERIVATIVE(X, I, ROT_DERIV)
                !PRINT *,'ROT_DERIV='
                !DO J = 1,NORBS
                !    PRINT *, ROT_DERIV(J,1:NORBS)
                !END DO

                CALL ELEMENTWISE_MULTIPLICATION(ROT_DERIV,PENALTY_DERIV,PROD_DERIVS)
                !PRINT *,'OVERALL_PROD='
                !DO J = 1,NORBS
                !    PRINT *, PROD_DERIVS(J,:)
                !END DO
                GRAD(I) = SUM(PROD_DERIVS)
                !PRINT *,'GRAD',I,'=',GRAD(I)
            END DO

        END IF

        IF (SECT) THEN
            ! If haven't calculated gradient as well need to calculate additional intermediates.
            IF (.NOT.GTEST) THEN
                CALL GET_INDIVIDUAL_DERIVATIVES(ROTATION, ORBDIPOLES, IND_ORBVAR_DERIV)
                CALL GET_PENALTY_DERIVATIVES(IND_ORBVAR_DERIV, ORBVAR, PENALTY_DERIV)
            END IF

            CALL CALC_HESSIAN(X,ROTATION,ORBVAR,ORBDIPOLES,IND_ORBVAR_DERIV,PENALTY_DERIV)

        END IF
        !PRINT*,"COORDS=",X,",GRAD=",GRAD
    END SUBROUTINE GET_ORBITAL_LOCALITY

    SUBROUTINE CHECK_COORDINATES(COORDS)

        USE COMMONS, ONLY: NROTS

        DOUBLE PRECISION, INTENT(INOUT) :: COORDS(1:NROTS)

        INTEGER :: I, NDIFF

        DOUBLE PRECISION, PARAMETER :: PI = 3.141592654D0
        DOUBLE PRECISION :: COORDSSAVE(1:NROTS)

        COORDSSAVE = COORDS

        DO I=1,NROTS
            IF (COORDS(I) > PI) THEN
                NDIFF = CEILING((COORDS(I)+PI)/(2*PI))
                COORDS(I) = COORDS(I) - 2*PI*(NDIFF-1)
!                PRINT*,COORDSSAVE(I),'->',COORDS(I),SIN(COORDSSAVE(I))-SIN(COORDS(I))
            ELSE IF (COORDS(I) < -PI) THEN
                NDIFF = CEILING((-COORDS(I)+PI)/(2*PI))
                COORDS(I) = COORDS(I) + 2*PI*(NDIFF-1)
!                PRINT*,COORDSSAVE(I),'->',COORDS(I),SIN(COORDSSAVE(I))-SIN(COORDS(I))
            END IF
        END DO

    END SUBROUTINE CHECK_COORDINATES

    SUBROUTINE ELEMENTWISE_MULTIPLICATION(MATA,MATB,MATC)

        USE COMMONS, ONLY: NORBS

        DOUBLE PRECISION, INTENT(IN) :: MATA(NORBS,NORBS), MATB(NORBS,NORBS)

        DOUBLE PRECISION, INTENT(OUT) :: MATC(NORBS, NORBS)

        INTEGER :: I, J

        DO I = 1, NORBS
            DO J = 1, NORBS
                MATC(I,J) = MATA(I,J) * MATB(I,J)
            END DO
        END DO

    END SUBROUTINE ELEMENTWISE_MULTIPLICATION

    SUBROUTINE GET_ROTATION(COORDS, ROTATION)

        ! Obtains the overall rotation resulting from the direct product
        ! of the Givens rotations corresponding to the specified coordinates.

        USE COMMONS, ONLY: NROTS, NORBS

        DOUBLE PRECISION, INTENT(IN) :: COORDS(NROTS)

        DOUBLE PRECISION, INTENT(OUT) :: ROTATION(NORBS,NORBS)
        DOUBLE PRECISION :: NEXT_ROT(NORBS,NORBS)

        INTEGER :: I, J, TRIIND

        ROTATION(1:NORBS,1:NORBS) = 0.0

        DO TRIIND = 1, NORBS
            ROTATION(TRIIND, TRIIND) = 1.0
        END DO

        DO TRIIND = 1, NROTS
            CALL DECODE_TRIIND(TRIIND, I, J)
            CALL GET_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
            ROTATION = MATMUL(ROTATION,NEXT_ROT)
        END DO

    END SUBROUTINE GET_ROTATION

    SUBROUTINE GET_ROTATION_DERIVATIVE(COORDS, DERIVTRIIND, DERIV_ROTATION)

        ! Obtain the derivative of the rotation matrix with respect to a single
        ! rotation angle theta at the current coordinates.

        USE COMMONS, ONLY: NROTS, NORBS

        DOUBLE PRECISION, INTENT(IN) :: COORDS(NROTS)
        INTEGER, INTENT(IN) :: DERIVTRIIND

        DOUBLE PRECISION, INTENT(OUT) :: DERIV_ROTATION(NORBS,NORBS)

        DOUBLE PRECISION :: NEXT_ROT(NORBS,NORBS)
        INTEGER :: TRIIND, I, J

        DERIV_ROTATION(1:NORBS,1:NORBS) = 0.0

        DO TRIIND = 1, NORBS
            DERIV_ROTATION(TRIIND, TRIIND) = 1.0
        END DO

        DO TRIIND = 1, NROTS
            CALL DECODE_TRIIND(TRIIND, I, J)
            IF (TRIIND.EQ.DERIVTRIIND) THEN
                CALL GET_DERIVATIVE_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
            ELSE
                CALL GET_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
            END IF
            DERIV_ROTATION = MATMUL(DERIV_ROTATION,NEXT_ROT)
        END DO

        !DO TRIIND = 1, NORBS
        !    PRINT *,TRIIND,'DERIV_ROT=',DERIV_ROTATION(TRIIND,1:NORBS)
        !END DO

    END SUBROUTINE GET_ROTATION_DERIVATIVE

    SUBROUTINE GET_ROTATION_SECOND_DERIVATIVE(COORDS, DERIVTRIINDA, DERIVTRIINDB, SECDERIV_ROTATION)

        ! Obtain the derivative of the rotation matrix with respect to a single
        ! rotation angle theta at the current coordinates.

        USE COMMONS, ONLY: NROTS, NORBS

        DOUBLE PRECISION, INTENT(IN) :: COORDS(NROTS)
        INTEGER, INTENT(IN) :: DERIVTRIINDA, DERIVTRIINDB

        DOUBLE PRECISION, INTENT(OUT) :: SECDERIV_ROTATION(NORBS,NORBS)

        DOUBLE PRECISION :: NEXT_ROT(NORBS,NORBS)
        INTEGER :: TRIIND, I, J

        SECDERIV_ROTATION(1:NORBS,1:NORBS) = 0.0

        DO TRIIND = 1, NORBS
            SECDERIV_ROTATION(TRIIND, TRIIND) = 1.0
        END DO

        DO TRIIND = 1, NROTS
            CALL DECODE_TRIIND(TRIIND, I, J)
            IF (TRIIND.EQ.DERIVTRIINDA) THEN
                IF (DERIVTRIINDA.EQ.DERIVTRIINDB) THEN
                    CALL GET_SECOND_DERIVATIVE_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
                ELSE
                    CALL GET_DERIVATIVE_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
                END IF
            ELSE IF (TRIIND.EQ.DERIVTRIINDB) THEN
                CALL GET_DERIVATIVE_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
            ELSE
                CALL GET_GIVENS_ROTATION(I,J,COORDS(TRIIND),NEXT_ROT)
            END IF
            SECDERIV_ROTATION = MATMUL(SECDERIV_ROTATION,NEXT_ROT)
        END DO

        !DO TRIIND = 1, NORBS
        !    PRINT *,TRIIND,'DERIV_ROT=',SECDERIV_ROTATION(TRIIND,1:NORBS)
        !END DO

    END SUBROUTINE GET_ROTATION_SECOND_DERIVATIVE

    SUBROUTINE GET_GIVENS_ROTATION(I, J, THETA, MAT)

        ! Obtains Givens rotation between orbitals I and J with rotation angle theta.

        USE COMMONS, ONLY: NORBS

        INTEGER, INTENT(IN) :: I, J
        DOUBLE PRECISION, INTENT(IN) :: THETA
        DOUBLE PRECISION, INTENT(OUT) :: MAT(NORBS,NORBS)

        INTEGER :: VAL

        MAT(1:NORBS,1:NORBS) = 0.0

        DO VAL = 1, NORBS
            MAT(VAL, VAL) = 1.0
        END DO

        MAT(I,I) = COS(THETA)
        MAT(J,J) = COS(THETA)
        MAT(I,J) = -SIN(THETA)
        MAT(J,I) = SIN(THETA)

        !PRINT *, 'GIVENSROT',I,'-',J,'='
        !DO VAL=1,NORBS
        !    PRINT *,MAT(VAL,1:NORBS)
        !END DO

    END SUBROUTINE GET_GIVENS_ROTATION

    SUBROUTINE GET_DERIVATIVE_GIVENS_ROTATION(I, J, THETA, MAT)

        ! Obtains the first derivative of the Givens rotation between orbitals I and J
        !  wrt the rotation angle at angle theta as a matrix.

        USE COMMONS, ONLY: NORBS

        INTEGER, INTENT(IN) :: I, J
        DOUBLE PRECISION, INTENT(IN) :: THETA
        DOUBLE PRECISION, INTENT(OUT) :: MAT(NORBS,NORBS)

        INTEGER :: VAL

        MAT(1:NORBS,1:NORBS) = 0.0

        MAT(I,I) = -SIN(THETA)
        MAT(J,J) = -SIN(THETA)
        MAT(I,J) = -COS(THETA)
        MAT(J,I) = COS(THETA)
        !PRINT *, 'DERIVGIVENSROT',I,'-',J,'='
        !DO VAL=1,NORBS
        !    PRINT *,MAT(VAL,1:NORBS)
        !END DO

    END SUBROUTINE GET_DERIVATIVE_GIVENS_ROTATION

    SUBROUTINE GET_SECOND_DERIVATIVE_GIVENS_ROTATION(I, J, THETA, MAT)

        ! Obtains the second derivative of the Givens rotation between orbitals I and J
        !  wrt the rotation angle at angle theta as a matrix.

        USE COMMONS, ONLY: NORBS

        INTEGER, INTENT(IN) :: I, J
        DOUBLE PRECISION, INTENT(IN) :: THETA
        DOUBLE PRECISION, INTENT(OUT) :: MAT(NORBS,NORBS)

        INTEGER :: VAL

        MAT(1:NORBS,1:NORBS) = 0.0

        MAT(I,I) = -COS(THETA)
        MAT(J,J) = -COS(THETA)
        MAT(I,J) = SIN(THETA)
        MAT(J,I) = -SIN(THETA)
        !PRINT *, 'SECONDDERIVGIVENSROT',I,'-',J,'='
        !DO VAL=1,NORBS
        !    PRINT *,MAT(VAL,1:NORBS)
        !END DO

    END SUBROUTINE GET_SECOND_DERIVATIVE_GIVENS_ROTATION

    SUBROUTINE GET_CURRENT_ORB_DIPOLES(ROTATION, ORBDIPOLES)

        ! Obtain the dipole moments of the orbital set created by the current rotation of the starting MOs.

        USE COMMONS, ONLY: NORBS, DIPINTS

        DOUBLE PRECISION, INTENT(IN) :: ROTATION(NORBS, NORBS)
        DOUBLE PRECISION :: ORBDIPOLES(3,NORBS)
        ! DOUBLE PRECISION :: TEMPCURRENTDIPOLES(NORBS,NORBS)

        INTEGER :: CURRENT_ORB, ORIG_ORBA, ORIG_ORBB, DIRECTION

        ORBDIPOLES(1:3,1:NORBS) = 0.0

        DO DIRECTION = 1,3
            DO CURRENT_ORB = 1, NORBS
                DO ORIG_ORBA = 1, NORBS
                    DO ORIG_ORBB = 1, NORBS
                        ORBDIPOLES(DIRECTION,CURRENT_ORB) = ORBDIPOLES(DIRECTION,CURRENT_ORB) + &
                            ROTATION(ORIG_ORBA,CURRENT_ORB) * ROTATION(ORIG_ORBB,CURRENT_ORB) * &
                            DIPINTS(DIRECTION,ORIG_ORBA,ORIG_ORBB)
                    END DO
                END DO
            END DO
        END DO

        ! Comparison using matrix multiplication.
        !DO DIRECTION = 1,3
        !    PRINT *, 'Prev and matmul dipoles, direction', DIRECTION, '.'
        !    TEMPCURRENTDIPOLES = MATMUL(MATMUL(TRANSPOSE(ROTATION),DIPINTS(DIRECTION,1:NORBS,1:NORBS)),ROTATION)
        !    DO CURRENT_ORB = 1, NORBS
        !        PRINT *, ORBDIPOLES(DIRECTION, CURRENT_ORB), TEMPCURRENTDIPOLES(CURRENT_ORB, CURRENT_ORB)
        !    END DO
        !END DO

    END SUBROUTINE GET_CURRENT_ORB_DIPOLES

    SUBROUTINE DECODE_TRIIND(TRIIND, I, J)

        ! Given a triangular index obtains the 2D indices corresponding to it.

        USE COMMONS, ONLY: NORBS

        INTEGER, INTENT(IN) :: TRIIND
        INTEGER, INTENT(OUT) :: I, J

        INTEGER :: RUNNING_SUM, VAL, I_, J_

        RUNNING_SUM = 0

        DO VAL = 1, NORBS
            IF (RUNNING_SUM + VAL .GE. TRIIND) THEN
                I_ = VAL + 1
                J_ = TRIIND - RUNNING_SUM
                EXIT
            ELSE
                RUNNING_SUM = RUNNING_SUM + VAL
            END IF

        END DO

        IF (I_ < J_) THEN
            I = J_
            J = I_
        ELSE
            I = I_
            J = J_
        END IF

        !PRINT '(A,3I10)','TRIIND,I,J=',TRIIND,I,J

    END SUBROUTINE DECODE_TRIIND

    SUBROUTINE GET_INDIVIDUAL_DERIVATIVES(CURRENT_ROTATION, CURRENT_DIP, INDIVIDUAL_DERIVATIVES)

        ! Generates the matrix of all derivatives of different orbital variances w.r.t.
        ! the coefficients of their expansion in the original MO basis.

        USE COMMONS, ONLY: NORBS

        DOUBLE PRECISION, INTENT(IN) :: CURRENT_ROTATION(NORBS, NORBS), CURRENT_DIP(3,NORBS)

        DOUBLE PRECISION, INTENT(OUT) :: INDIVIDUAL_DERIVATIVES(NORBS,NORBS)

        INTEGER :: ORIGORB, NEWORB

        !PRINT *,'CURRENTROTATION=',CURRENT_ROTATION

        DO NEWORB = 1, NORBS
            DO ORIGORB = 1, NORBS
                CALL GET_SINGLE_ORBVAR_DERIVATIVE(CURRENT_ROTATION, CURRENT_DIP, NEWORB, ORIGORB, INDIVIDUAL_DERIVATIVES(ORIGORB,NEWORB))
            END DO
        END DO

    END SUBROUTINE GET_INDIVIDUAL_DERIVATIVES

    SUBROUTINE GET_SINGLE_ORBVAR_DERIVATIVE(CURRENT_ROTATION, CURRENT_DIP, NEW_ORB, ORIG_ORB, DERIVATIVE)

        ! Obtains the derivative of a single orbital at the current rotation's orbital
        ! variance w.r.t. to the coefficient an original (unrotated) orbital.

        USE COMMONS, ONLY: NORBS, R2INTS, DIPINTS

        DOUBLE PRECISION, INTENT(IN) :: CURRENT_ROTATION(NORBS, NORBS), CURRENT_DIP(3,NORBS)

        DOUBLE PRECISION, INTENT(OUT) :: DERIVATIVE

        INTEGER,INTENT(IN) :: NEW_ORB, ORIG_ORB

        INTEGER :: OTHER_ORIG_ORB

        DOUBLE PRECISION :: CONTRIB

        DERIVATIVE = 0.0

        DO OTHER_ORIG_ORB = 1, NORBS
            CONTRIB = 2 * CURRENT_ROTATION(OTHER_ORIG_ORB,NEW_ORB)*R2INTS(ORIG_ORB, OTHER_ORIG_ORB)

            DERIVATIVE = DERIVATIVE + CONTRIB


            CONTRIB = 4 * CURRENT_ROTATION(OTHER_ORIG_ORB,NEW_ORB)*&
                        (CURRENT_DIP(1,NEW_ORB)*DIPINTS(1,ORIG_ORB,OTHER_ORIG_ORB)&
                        +CURRENT_DIP(2,NEW_ORB)*DIPINTS(2,ORIG_ORB,OTHER_ORIG_ORB)&
                        +CURRENT_DIP(3,NEW_ORB)*DIPINTS(3,ORIG_ORB,OTHER_ORIG_ORB))

            DERIVATIVE = DERIVATIVE - CONTRIB

        END DO

    END SUBROUTINE GET_SINGLE_ORBVAR_DERIVATIVE

    SUBROUTINE GET_PENALTY_DERIVATIVES(IND_ORBVAR_DERIV, ORB_VAR, PENALTY_DERIVATIVES)

        ! Get the derivative of the overall penalty function with respect to each element
        ! of the rotation matrix.

        USE COMMONS, ONLY: NORBS, ORBVAREXPONENT

        DOUBLE PRECISION, INTENT(IN) :: IND_ORBVAR_DERIV(NORBS,NORBS), ORB_VAR(NORBS)
        DOUBLE PRECISION, INTENT(OUT) :: PENALTY_DERIVATIVES(NORBS, NORBS)

        INTEGER :: ORIG_ORB, NEW_ORB

        PENALTY_DERIVATIVES = 0.0

        DO ORIG_ORB = 1, NORBS
            DO NEW_ORB = 1, NORBS
                PENALTY_DERIVATIVES(ORIG_ORB,NEW_ORB) = &
                    ORBVAREXPONENT * IND_ORBVAR_DERIV(ORIG_ORB,NEW_ORB) * ORB_VAR(NEW_ORB) ** (ORBVAREXPONENT-1)
            END DO
        END DO

    END SUBROUTINE GET_PENALTY_DERIVATIVES

    SUBROUTINE CALC_HESSIAN(COORDS,ROTATION,ORBVAR,ORBDIPOLES,IND_ORBVAR_DERIV,PENALTY_DERIV)

        USE COMMONS, ONLY: NROTS, NORBS
        USE MODHESS

        DOUBLE PRECISION, INTENT(IN) :: COORDS(NROTS)
        DOUBLE PRECISION, INTENT(IN) :: ROTATION(NORBS,NORBS), ORBDIPOLES(3,NORBS), IND_ORBVAR_DERIV(NORBS,NORBS)
        DOUBLE PRECISION, INTENT(IN) :: ORBVAR(NORBS), PENALTY_DERIV(NORBS,NORBS)

        DOUBLE PRECISION :: IND_ORBVAR_SECOND_DERIV(NORBS,NORBS,NORBS)
        DOUBLE PRECISION :: PENALTY_SECOND_DERIV(NORBS,NORBS,NORBS), TEMP
        DOUBLE PRECISION :: ROTMATA(NORBS,NORBS), ROTMATB(NORBS,NORBS)

        INTEGER :: TRIINDA, TRIINDB, NEW_ORB, ORIG_ORBA, ORIG_ORBB

        !IF (.NOT.ALLOCATED(HESS)) THEN
        !    PRINT *, 'Allocating Hessian with size ',NROTS,'x',NROTS
        !    ALLOCATE(HESS(NROTS,NROTS))
        !END IF
        !PRINT*,'hessian bounds:',LBOUND(HESS,DIM=1),UBOUND(HESS,DIM=1),LBOUND(HESS,DIM=2),UBOUND(HESS,DIM=2)
        !PRINT*,'hessian size:',SIZE(HESS),SHAPE(HESS)
        ! Gradually work our way up the chain rule...
        CALL CALC_ORBVAR_SECOND_DERIV(ROTATION,ORBDIPOLES,IND_ORBVAR_SECOND_DERIV)

        CALL CALC_PENALTY_SECOND_DERIV(ORBVAR,IND_ORBVAR_DERIV,IND_ORBVAR_SECOND_DERIV,PENALTY_SECOND_DERIV)


        !PRINT*,NROTS,NORBS

        !PRINT*,'IND_ORBVAR_SECOND_DERIV='
        !DO TRIINDA=1,NORBS
        !    DO TRIINDB=1,NORBS
        !        PRINT*,TRIINDA,TRIINDB,IND_ORBVAR_SECOND_DERIV(TRIINDA,TRIINDB,1:NORBS)
        !    END DO
        !END DO
        !PRINT*,'PENALTY_SECOND_DERIV='
        !DO TRIINDA=1,NORBS
        !    DO TRIINDB=1,NORBS
        !        PRINT*,TRIINDA,TRIINDB,PENALTY_SECOND_DERIV(TRIINDA,TRIINDB,1:NORBS)
        !    END DO
        !END DO

        HESS(1:NROTS,1:NROTS) = 0.0

        ! Now calculate actual hessian!
        DO TRIINDA=1,NROTS
            DO TRIINDB=1,NROTS
                ! First deal with term proportional to penalty first deriv wrt coeffs
                CALL GET_ROTATION_SECOND_DERIVATIVE(COORDS,TRIINDA,TRIINDB,ROTMATA)

                ! Use ROTMATB as scratch space.
                CALL ELEMENTWISE_MULTIPLICATION(ROTMATA(1:NORBS,1:NORBS),PENALTY_DERIV(1:NORBS,1:NORBS),ROTMATB)
                !PRINT*,'PRODMAT SECDERIV='
                !DO NEW_ORB = 1,NORBS
                !    PRINT*,ROTMATB(NEW_ORB,1:NORBS)
                !END DO

                HESS(TRIINDA,TRIINDB) = SUM(ROTMATB)
                ! Now deal with term proportional to penalty second deriv wrt coeffs.
                CALL GET_ROTATION_DERIVATIVE(COORDS,TRIINDA,ROTMATA)
                CALL GET_ROTATION_DERIVATIVE(COORDS,TRIINDB,ROTMATB)
                DO NEW_ORB = 1,NORBS
                    DO ORIG_ORBA = 1,NORBS
                        DO ORIG_ORBB = 1,NORBS
                            HESS(TRIINDA,TRIINDB) = HESS(TRIINDA,TRIINDB) + PENALTY_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB)*ROTMATA(ORIG_ORBA,NEW_ORB)*ROTMATB(ORIG_ORBB,NEW_ORB)
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE CALC_HESSIAN

    SUBROUTINE CALC_ORBVAR_SECOND_DERIV(ROTATION,ORBDIPOLES,IND_ORBVAR_SECOND_DERIV)

        USE COMMONS, ONLY: NORBS, R2INTS, DIPINTS

        DOUBLE PRECISION, INTENT(IN) :: ROTATION(NORBS,NORBS), ORBDIPOLES(3,NORBS)
        DOUBLE PRECISION, INTENT(OUT) :: IND_ORBVAR_SECOND_DERIV(NORBS,NORBS,NORBS)

        DOUBLE PRECISION :: INTERMEDIATE_SUM(NORBS)

        INTEGER :: ORIG_ORBA, ORIG_ORBB, NEW_ORB, DIRECTION

        ! First set values according to [x2] term.
        DO ORIG_ORBA = 1,NORBS
            DO ORIG_ORBB = 1,NORBS
                IND_ORBVAR_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,1:NORBS) = 2 * R2INTS(ORIG_ORBA,ORIG_ORBB)
            END DO
        END DO

        DO NEW_ORB = 1, NORBS
            DO DIRECTION = 1,3
                ! Calc \sum_i C_ip [x]_ij as for each j value (and given p).
                INTERMEDIATE_SUM(1:NORBS) = 0.0
                DO ORIG_ORBA = 1,NORBS
                    DO ORIG_ORBB = 1,NORBS
                        INTERMEDIATE_SUM(ORIG_ORBA) = INTERMEDIATE_SUM(ORIG_ORBA) + ROTATION(ORIG_ORBB,NEW_ORB) * DIPINTS(DIRECTION,ORIG_ORBA,ORIG_ORBB)
                    END DO
                END DO
                !PRINT*,'DIRECTION=',DIRECTION
                !PRINT*,'DIPINTS=',DIPINTS(DIRECTION,:,:)
                !PRINT*,'INTERMEDIATE_SUM=',INTERMEDIATE_SUM
                DO ORIG_ORBA = 1,NORBS
                    DO ORIG_ORBB = 1,NORBS
                        IND_ORBVAR_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB) = IND_ORBVAR_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB) &
                                            - 8*INTERMEDIATE_SUM(ORIG_ORBA)*INTERMEDIATE_SUM(ORIG_ORBB) &
                                            - 4*ORBDIPOLES(DIRECTION,NEW_ORB)*DIPINTS(DIRECTION,ORIG_ORBA,ORIG_ORBB)
                    END DO
                END DO

            END DO
        END DO

    END SUBROUTINE CALC_ORBVAR_SECOND_DERIV

    SUBROUTINE CALC_PENALTY_SECOND_DERIV(ORBVAR,IND_ORBVAR_DERIV,IND_ORBVAR_SECOND_DERIV,PENALTY_SECOND_DERIV)

        USE COMMONS, ONLY: ORBVAREXPONENT, NORBS

        DOUBLE PRECISION, INTENT(IN) :: ORBVAR(NORBS), IND_ORBVAR_DERIV(NORBS,NORBS)

        DOUBLE PRECISION, INTENT(IN) :: IND_ORBVAR_SECOND_DERIV(NORBS,NORBS,NORBS)
        DOUBLE PRECISION, INTENT(OUT) :: PENALTY_SECOND_DERIV(NORBS,NORBS,NORBS)

        INTEGER :: NEW_ORB, ORIG_ORBA, ORIG_ORBB

        DO NEW_ORB = 1,NORBS
            DO ORIG_ORBA = 1,NORBS
                DO ORIG_ORBB = 1,NORBS
                    PENALTY_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB) = &
                                ORBVAREXPONENT * IND_ORBVAR_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB)*(ORBVAR(NEW_ORB))**(ORBVAREXPONENT-1)
                    IF (ORBVAREXPONENT.GT.1) THEN
                        PENALTY_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB) = PENALTY_SECOND_DERIV(ORIG_ORBA,ORIG_ORBB,NEW_ORB)&
                                +ORBVAREXPONENT*(ORBVAREXPONENT-1)*IND_ORBVAR_DERIV(ORIG_ORBA,NEW_ORB)*IND_ORBVAR_DERIV(ORIG_ORBB,NEW_ORB)*(ORBVAR(NEW_ORB))**(ORBVAREXPONENT-2)
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE CALC_PENALTY_SECOND_DERIV

    SUBROUTINE ORBITALS_FINISH()

        ! Deallocate arrays in commons.

        USE COMMONS, ONLY: R2INTS, DIPINTS

        IF (ALLOCATED(R2INTS)) DEALLOCATE(R2INTS)
        IF (ALLOCATED(DIPINTS)) DEALLOCATE(DIPINTS)

    END SUBROUTINE ORBITALS_FINISH

    SUBROUTINE CHECK_FILE_EXISTS(FNAME, EX)

        ! Check if a file of the given name exists, exit if it does not and EX=.TRUE.

        LOGICAL :: FILE_EXISTS, TEX
        LOGICAL, INTENT(IN), OPTIONAL :: EX
        CHARACTER(*), INTENT(IN) :: FNAME

        IF(PRESENT(EX)) THEN
          TEX = EX
        ELSE
          TEX = .FALSE.
        ENDIF

        INQUIRE(FILE=FNAME, EXIST=FILE_EXISTS)
        IF(.NOT. FILE_EXISTS) THEN
          WRITE(*,'(A)') "The file " // TRIM(FNAME) // " does not exist."
          IF(TEX.EQV..TRUE.) THEN
            STOP
          ENDIF
        ENDIF

    END SUBROUTINE CHECK_FILE_EXISTS

END MODULE ORBITALS_MOD
