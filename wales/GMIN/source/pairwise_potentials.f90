! This module provides a library of potential functions which operate only on a pair of atoms
! It is designed to be used with the MULTIPOT module

!All subroutines in this module must have the following signature to be used with MULTIPOT:
!        SUBROUTINE POT(X1, X2, PARAMS, PG, PAIR_ENERGY, P_HESS, GTEST, STEST)
!            DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3)
!            DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)      ! Maximum number of parameters is hardcoded here
!                                                             ! Each potential will use a subset of the elements of PARAMS
!            DOUBLE PRECISION, INTENT(OUT) :: PAIR_ENERGY
!            DOUBLE PRECISION, INTENT(OUT) :: PG(6), P_HESS(9)! Gradient and energy for the subsystem composed of this pair
!            LOGICAL, INTENT(IN)           :: GTEST, STEST    ! GTEST is true when calculating the gradient as well as energy.
!                                                             ! STEST is true when calculating the Hessian as well as the other two.
!        END SUBROUTINE POT

MODULE PAIRWISE_POTENTIALS

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: WCA_CUT=1.2599210498948731647672106072782283505702514647015079 ! = 2^(1/3)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simple LJ potential, no cutoff. Atom radius sigma is given as a variable, but will often be set to 1.
SUBROUTINE PAIRWISE_LJ(X1, X2, PARAMS, PG, PAIR_ENERGY, P_HESS, GTEST, STEST)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3), PARAMS(10)
    DOUBLE PRECISION, INTENT(OUT) :: PAIR_ENERGY
    DOUBLE PRECISION, INTENT(OUT) :: PG(6), P_HESS(6,6)
    LOGICAL, INTENT(IN)           :: GTEST, STEST

    DOUBLE PRECISION :: R2, R6, R8, R14, SIG, SIG6, SIG12  ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: G, F  ! G tensor and F tensor (see PAIRWISE_GRAD and PAIRWISE_HESSIAN, below)

    SIG = PARAMS(1)
    SIG6 = SIG**6
    SIG12 = SIG6**2

    R2 = (X1(1)-X2(1))**2+(X1(2)-X2(2))**2+(X1(3)-X2(3))**2
    R2 = 1.0D0/R2
    R6 = R2**3

    PAIR_ENERGY=4.0D0*SIG6*R6*(SIG6*R6-1.0D0)

    IF (.NOT.GTEST) RETURN

    G = -24.0D0*(2.0D0*SIG6*R6-1.0D0)*R2*SIG6*R6

    CALL PAIRWISE_GRAD(X1, X2, G, PG)

    IF (.NOT.STEST) RETURN

    R8 = R2**4
    R14 = R8*R8/R2
    F = 672.0D0*SIG12*R14-192.0D0*SIG6*R8

    CALL PAIRWISE_HESSIAN(X1,X2,G,F,R2,P_HESS)

    RETURN

END SUBROUTINE PAIRWISE_LJ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WCA potential, which is simply the LJ potential but truncated at the first minimum of the function and shifted so that
! the potential is entirely repulsive.
! The energy and gradient are both continuous at the cutoff.
! Atom radius sigma is given as a variable (SIG), but will often be set to 1.
SUBROUTINE PAIRWISE_WCA(X1, X2, PARAMS, PG, PAIR_ENERGY, P_HESS, GTEST, STEST)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3), PARAMS(10) ! Maximum number of params is hardcoded here
    DOUBLE PRECISION, INTENT(OUT) :: PAIR_ENERGY
    DOUBLE PRECISION, INTENT(OUT) :: PG(6), P_HESS(6,6)
    LOGICAL, INTENT(IN)           :: GTEST, STEST  ! GTEST is true when calculating the gradient as well as energy.
                                                   ! STEST is true when calculating the Hessian as well as the other two.

    DOUBLE PRECISION :: R2, R6, R8, R14, SIG, SIG6, SIG12  ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: G, F  ! G tensor and F tensor (see PAIRWISE_GRAD and PAIRWISE_HESSIAN, below)

    SIG = PARAMS(1)
    SIG6 = SIG**6
    SIG12 = SIG6**2

    R2 = (X1(1)-X2(1))**2+(X1(2)-X2(2))**2+(X1(3)-X2(3))**2

    IF(R2.GT.WCA_CUT*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
        PAIR_ENERGY = 0.0D0
        PG(:) = 0.0D0
        P_HESS(:,:) = 0.0D0
        RETURN
    ENDIF

    R2 = 1.0D0/R2
    R6 = R2**3

    PAIR_ENERGY=4.0D0*SIG6*R6*(SIG6*R6-1.0D0) + 1

    IF (.NOT.GTEST) RETURN

    G = -24.0D0*(2.0D0*SIG6*R6-1.0D0)*R2*SIG6*R6

    CALL PAIRWISE_GRAD(X1, X2, G, PG)

    IF (.NOT.STEST) RETURN

    R8 = R2**4
    R14 = R8*R8/R2
    F = 672.0D0*SIG12*R14-192.0D0*SIG6*R8

    CALL PAIRWISE_HESSIAN(X1,X2,G,F,R2,P_HESS)

    RETURN

END SUBROUTINE PAIRWISE_WCA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hooke's Law spring potential
! The equilibrium bond length is given as PARAMS(1) (saved as R0). Energy is given in units of the spring constant.
SUBROUTINE HARMONIC_SPRINGS(X1, X2, PARAMS, PG, PAIR_ENERGY, P_HESS, GTEST, STEST)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3), PARAMS(10)
    DOUBLE PRECISION, INTENT(OUT) :: PAIR_ENERGY
    DOUBLE PRECISION, INTENT(OUT) :: PG(6), P_HESS(6,6)
    LOGICAL, INTENT(IN)           :: GTEST, STEST  ! GTEST is true when calculating the gradient as well as energy.
                                                   ! STEST is true when calculating the Hessian as well as the other two.
    DOUBLE PRECISION :: R0, R, R2, G, F

    R0 = PARAMS(1)
    R2 = (X1(1)-X2(1))**2+(X1(2)-X2(2))**2+(X1(3)-X2(3))**2
    R = SQRT(R2)

    PAIR_ENERGY = 0.5D0*(R-R0)**2

    IF(.NOT.GTEST) RETURN

    F = R0/R
    G = 1.0D0 - F
    R2 = 1.0D0/R2      ! NOTE: R2 must be the INVERSE square distance between the two particles (cf. LJ and WCA potentials)

    CALL PAIRWISE_GRAD(X1,X2,G,PG)

        IF (.NOT.STEST) RETURN

    CALL PAIRWISE_HESSIAN(X1,X2,G,F,R2,P_HESS)

    RETURN

END SUBROUTINE HARMONIC_SPRINGS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END OF PAIRWISE POTENTIALS
! Next, two functions that are general to all isotropic pairwise potentials. To be clear, by "isotropic", I mean
! "depends only on the distance between the two particles"
! Another note: I refer to the "G and F tensors", to be in line with comments in other source files. In the
! present case, where there are only 2 atoms, the tensor is represented by a 2x2 matrix in which two elements are 0 and
! the other two are equal. So I haven't bothered to represent it as a real tensor, but simply as a single scalar for the
! unique value in the matrix.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For all isotropic potentials the gradient is calculated in the same way. We simply need to know the positions of the
! two particles, and the G tensor: G = (1/r)dU/dr for an isotropic potential U(r)
SUBROUTINE PAIRWISE_GRAD(X1, X2, G, PG)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3), G
    DOUBLE PRECISION, INTENT(OUT) :: PG(6)

    ! Gradient with respect to atom 1 coordinates
    PG(1) = G*(X1(1)-X2(1))
    PG(2) = G*(X1(2)-X2(2))
    PG(3) = G*(X1(3)-X2(3))
    ! Gradient with respect to atom 2 coordinates
    PG(4) = G*(X2(1)-X1(1))
    PG(5) = G*(X2(2)-X1(2))
    PG(6) = G*(X2(3)-X1(3))

    RETURN
END SUBROUTINE PAIRWISE_GRAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For all isotropic potentials the Hessian is calculated in the same way. We simply need to know the positions of the
! two particles, the G tensor (see above) and the F tensor: F = r*d[(1/r)dU/dr]/dr
! NOTE: R2 must be the INVERSE square distance between the two particles (cf. LJ and WCA potentials)
SUBROUTINE PAIRWISE_HESSIAN(X1,X2,G,F,R2,P_HESS)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3), G, F, R2
    DOUBLE PRECISION, INTENT(OUT) :: P_HESS(6,6)
    INTEGER :: J1,J2

!       Compute the hessian. First are the entirely diagonal terms (6 of them)
!       These correspond to 2nd derivatives wrt the same coordinate twice
    DO J1=1,3
        P_HESS(J1,J1) = F*R2*(X1(J1)-X2(J1))**2 + G
        P_HESS(3+J1,3+J1) = F*R2*(X2(J1)-X1(J1))**2 + G
    ENDDO

!       Next the terms where x_i and x_j are on the same atom but different cartesian directions
!       (12 of these but we only calculate the upper RH block, then symmetrise later)
    DO J1=1,3
        DO J2=J1+1,3
            P_HESS(J1,J2) = F*R2*(X1(J1)-X2(J1))*(X1(J2)-X2(J2))
            P_HESS(3+J1,3+J2) = F*R2*(X2(J1)-X1(J1))*(X2(J2)-X1(J2))
        ENDDO
    ENDDO

!       Different atoms, same cartesian direction (6 of these but we only calculate the upper RH block, then symmetrise later)
    DO J1=1,3
        P_HESS(J1,3+J1) = -F*R2*(X1(J1)-X2(J1))**2 - G
    ENDDO

!       Finally, different atoms and different coordinates (12 of these, we only calculate 6)
    DO J1=1,3
        DO J2=J1+1,3
            P_HESS(J1,3+J2) = -F*R2*(X1(J1)-X2(J1))*(X1(J2)-X2(J2))
            P_HESS(J2,3+J1) = P_HESS(J1,3+J2)
        ENDDO
    ENDDO

!     Symmetrise Hessian
    P_HESS(2,1)=P_HESS(1,2)
    P_HESS(3,1)=P_HESS(1,3)
    P_HESS(3,2)=P_HESS(2,3)
    P_HESS(4,1)=P_HESS(1,4)
    P_HESS(4,2)=P_HESS(2,4)
    P_HESS(4,3)=P_HESS(3,4)
    P_HESS(5,1)=P_HESS(1,5)
    P_HESS(5,2)=P_HESS(2,5)
    P_HESS(5,3)=P_HESS(3,5)
    P_HESS(5,4)=P_HESS(4,5)
    P_HESS(6,1)=P_HESS(1,6)
    P_HESS(6,2)=P_HESS(2,6)
    P_HESS(6,3)=P_HESS(3,6)
    P_HESS(6,4)=P_HESS(4,6)
    P_HESS(6,5)=P_HESS(5,6)

!    DO J1=1,6
!       DO J2=J1+1,6
!          P_HESS(J2,J1)=P_HESS(J1,J2)
!       ENDDO
!    ENDDO

    RETURN
END SUBROUTINE PAIRWISE_HESSIAN

END MODULE PAIRWISE_POTENTIALS
