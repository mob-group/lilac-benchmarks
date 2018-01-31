! This module provides a library of potential functions for isotropic potentials which return a local copy of the Hessian.
! Here, "isotropic" means that the energy depends only on the distance between particles.
! It is designed to be used with the MULTIPOT module

!All subroutines in this module must have the following signature to be used with MULTIPOT:
!        SUBROUTINE POT(TMP_NATOMS, X, PARAMS, TMP_ENERGY, TMP_G, TMP_HESS, GTEST, STEST)
!            INTEGER, INTENT(IN)           :: TMP_NATOMS
!            DOUBLE PRECISION, INTENT(IN)  :: X(3*TMP_NATOMS)
!            DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)      ! Maximum number of parameters is hardcoded here
!                                                             ! Each potential will use a subset of the elements of PARAMS
!            DOUBLE PRECISION, INTENT(OUT) :: TMP_ENERGY
!            DOUBLE PRECISION, INTENT(OUT) :: TMP_G(3*TMP_NATOMS), TMP_HESS(3*TMP_NATOMS,3*TMP_NATOMS)
!            LOGICAL, INTENT(IN)           :: GTEST, STEST    ! GTEST is true when calculating the gradient as well as energy.
!                                                             ! STEST is true when calculating the Hessian as well as the other two.
!        END SUBROUTINE POT

! An exception to this rule is the EXCLUDE potentials, which need additional information passed in. They are handled differently
! (see multipot.f90)

MODULE ISOTROPIC_POTENTIALS

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: WCA_CUT=1.2599210498948731647672106072782283505702514647015079 ! = 2^(1/3)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simple LJ potential, no cutoff. Atom radius sigma is given as a variable, but will often be set to 1.
! This is almost an exact copy of ljdiff.f, but with sigma added in and the Hessian stored locally to be returned.
SUBROUTINE ISO_LJ(N, X, PARAMS, TMP_ENERGY, TMP_G, TMP_HESS, GTEST, STEST)
    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: N ! Number of atoms
    DOUBLE PRECISION, INTENT(IN)  :: X(3*N)
    DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)  ! Maximum number of parameters is hardcoded here
    DOUBLE PRECISION, INTENT(OUT) :: TMP_ENERGY
    DOUBLE PRECISION, INTENT(OUT) :: TMP_G(3*N), TMP_HESS(3*N,3*N)
    LOGICAL, INTENT(IN)           :: GTEST, STEST

    ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: DIST, R2(N,N), R6, R8(N,N), R14(N,N), SIG, SIG6, SIG12
    DOUBLE PRECISION :: G(N,N), F(N,N)  ! G tensor and F tensor (see ISOTROPIC_GRAD and ISOTROPIC_HESSIAN, below)
    INTEGER :: J1, J2, J3, J4

    SIG = PARAMS(1)
    SIG6 = SIG**6
    SIG12 = SIG6**2

    TMP_ENERGY=0.0D0

    ! The arrangement of these IF statements seems slightly odd, but it's been done to make sure we only need to set up one pair
    ! of DO loops for each call.
    IF (GTEST) THEN
        IF (STEST) THEN  ! Gradient + Hessian
            DO J1=1,N
                J3=3*J1
                R2(J1,J1)=0.0D0
                R8(J1,J1)=0.0D0
                R14(J1,J1)=0.0D0
                G(J1,J1)=0.0D0
                F(J1,J1)=0.0D0
                DO J2=J1+1,N
                    J4=3*J2

                    R2(J2,J1)=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                    R2(J2,J1)=1.0D0/R2(J2,J1)
                    R6=R2(J2,J1)**3
                    TMP_ENERGY=TMP_ENERGY+SIG6*R6*(SIG6*R6-1.0D0)

                    ! Set up storage arrays to use in the gradient- and hessian-calculating routines
                    R8(J2,J1)=R2(J2,J1)**4
                    R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
                    R2(J1,J2)=R2(J2,J1)
                    G(J2,J1)=-24.0D0*(2.0D0*SIG6*R6-1.0D0)*R2(J1,J2)*SIG6*R6
                    G(J1,J2)=G(J2,J1)
                    F(J2,J1)=672.0D0*R14(J2,J1)*SIG12-192.0D0*R8(J2,J1)*SIG6
                    F(J1,J2)=F(J2,J1)
                ENDDO
            ENDDO
        ELSE             ! Energy + Gradient only
            DO J1=1,N
                J3=3*J1
                G(J1,J1)=0.0D0
                DO J2=J1+1,N
                    J4=3*J2

                    DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                    DIST=1.0D0/DIST
                    R6=DIST**3

                    TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0)

                    G(J2,J1)=-24.0D0*(2.0D0*R6*SIG6-1.0D0)*DIST*R6*SIG6
                    G(J1,J2)=G(J2,J1)
                ENDDO
            ENDDO
        ENDIF
    ELSE                ! Energy only
        DO J1=1,N
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,N
                J4=3*J2

                DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                DIST=1.0D0/DIST
                R6=DIST**3

                TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0)
            ENDDO
        ENDDO

    ENDIF
    TMP_ENERGY=4.0D0*TMP_ENERGY

    IF (.NOT.GTEST) RETURN

    ! G should already be set
    CALL ISOTROPIC_GRAD(N, X, G, TMP_G)

    IF (.NOT.STEST) RETURN  ! It is assumed we will never need the Hessian without also needing the gradient.

    CALL ISOTROPIC_HESSIAN(N, X, G, F, R2, TMP_HESS)

    RETURN

END SUBROUTINE ISO_LJ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WCA potential, which is simply the LJ potential but truncated at the first minimum of the function and shifted so that
! the potential is entirely repulsive.
! The energy and gradient are both continuous at the cutoff.
! Atom radius sigma is given as a variable (SIG), but will often be set to 1.
SUBROUTINE ISO_WCA(N, X, PARAMS, TMP_ENERGY, TMP_G, TMP_HESS, GTEST, STEST)
    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: N ! Number of atoms
    DOUBLE PRECISION, INTENT(IN)  :: X(3*N)
    DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)  ! Maximum number of parameters is hardcoded here
    DOUBLE PRECISION, INTENT(OUT) :: TMP_ENERGY
    DOUBLE PRECISION, INTENT(OUT) :: TMP_G(3*N), TMP_HESS(3*N,3*N)
    LOGICAL, INTENT(IN)           :: GTEST, STEST

    ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: DIST, R2(N,N), R6, R8(N,N), R14(N,N), SIG, SIG6, SIG12
    DOUBLE PRECISION :: G(N,N), F(N,N)  ! G tensor and F tensor (see ISOTROPIC_GRAD and ISOTROPIC_HESSIAN, below)
    INTEGER :: J1, J2, J3, J4

    SIG = PARAMS(1)
    SIG6 = SIG**6
    SIG12 = SIG6**2

    TMP_ENERGY=0.0D0

    ! The arrangement of these IF statements seems slightly odd, but it's been done to make sure we only need to set up one pair
    ! of DO loops for each call.
    IF (GTEST) THEN
        IF (STEST) THEN  ! Gradient + Hessian
            DO J1=1,N
                J3=3*J1
                R2(J1,J1)=0.0D0
                R8(J1,J1)=0.0D0
                R14(J1,J1)=0.0D0
                G(J1,J1)=0.0D0
                F(J1,J1)=0.0D0
                DO J2=J1+1,N
                    J4=3*J2

                    R2(J2,J1)=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2

                    IF(R2(J2,J1).GT.WCA_CUT*SIG*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
                    ! We don't compute the energy for this pair of atoms. Set G and F to 0 so that the gradient and
                    ! Hessian terms will go to 0 also.
!                        write(*,*) "E+G+H. Separation>cutoff:", J1, J2, R2(J2,J1)
                        G(J1,J2) = 0.0D0
                        G(J2,J1) = 0.0D0
                        F(J1,J2) = 0.0D0
                        F(J2,J1) = 0.0D0
                        CYCLE
                    ENDIF

                    R2(J2,J1)=1.0D0/R2(J2,J1)
                    R6=R2(J2,J1)**3

                    TMP_ENERGY=TMP_ENERGY+SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
!                    write(*,*) "E+G+H. Pair, Dist, Energy:", J1, J2, 1.0D0/R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    ! Set up storage arrays to use in the gradient- and hessian-calculating routines
                    R8(J2,J1)=R2(J2,J1)**4
                    R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
                    R2(J1,J2)=R2(J2,J1)
                    G(J2,J1)=-24.0D0*(2.0D0*SIG6*R6-1.0D0)*R2(J1,J2)*SIG6*R6
                    G(J1,J2)=G(J2,J1)
                    F(J2,J1)=672.0D0*R14(J2,J1)*SIG12-192.0D0*R8(J2,J1)*SIG6
                    F(J1,J2)=F(J2,J1)
                ENDDO
            ENDDO
        ELSE             ! Energy + Gradient only
            DO J1=1,N
                J3=3*J1
                G(J1,J1)=0.0D0
                DO J2=J1+1,N
                    J4=3*J2

                    DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2

                    IF(DIST.GT.WCA_CUT*SIG*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
!                        write(*,*) "E+G. Separation>cutoff:", J1, J2, DIST
                        G(J1,J2) = 0.0D0
                        G(J2,J1) = 0.0D0
                        CYCLE
                    ENDIF

                    DIST=1.0D0/DIST
                    R6=DIST**3

                    TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0) + 0.25D0
!                    write(*,*) "E+G. Pair, Dist, Energy:", J1, J2, (X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2, SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    G(J2,J1)=-24.0D0*(2.0D0*R6*SIG6-1.0D0)*DIST*R6*SIG6
                    G(J1,J2)=G(J2,J1)
                ENDDO
            ENDDO
        ENDIF
    ELSE                ! Energy only
        DO J1=1,N
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,N
                J4=3*J2

                DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2

                IF(DIST.GT.WCA_CUT*SIG*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
!                    write(*,*) "E. Separation>cutoff:", J1, J2, R2(J2,J1)
                    CYCLE
                ENDIF

                DIST=1.0D0/DIST
                R6=DIST**3

                TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0) + 0.25D0
!                write(*,*) "E. Pair, Dist, Energy:", J1, J2, R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
            ENDDO
        ENDDO

    ENDIF
    TMP_ENERGY=4.0D0*TMP_ENERGY

    IF (.NOT.GTEST) RETURN

    ! G should already be set
    CALL ISOTROPIC_GRAD(N, X, G, TMP_G)

    IF (.NOT.STEST) RETURN  ! It is assumed we will never need the Hessian without also needing the gradient.

    CALL ISOTROPIC_HESSIAN(N, X, G, F, R2, TMP_HESS)

    RETURN

END SUBROUTINE ISO_WCA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coulomb potential (no cutoff - use for nonperiodic potentials only)
! This subroutine is called directly from multipot.f90, so that we can avoid calculating unnecessary interactions
! (see below). The same approach is used for the "EXCLUDE" potentials.
SUBROUTINE ISO_COULOMB(X, POTLIST, N_ATOM_PARTNERS, POTSCALE, PARAMS, ENERGY, GRAD, GTEST, STEST)

    ! Because we're passing in the entire coordinates array here, we need to know the total number of atoms
    USE COMMONS, ONLY: NATOMS

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: X(3*NATOMS)          ! The full atomistic coordinate array
    INTEGER, INTENT(IN)           :: POTLIST(:,:)         ! An array containing the interaction partners for each atom
    INTEGER, INTENT(IN)           :: N_ATOM_PARTNERS(:)   ! An array containing the number of interaction partners for each atom
    !DOUBLE PRECISION, INTENT(IN)  :: Q(:)                 ! An array containing the atom charges.
    DOUBLE PRECISION, INTENT(IN)  :: POTSCALE             ! The energy unit for this potential, which multiplies E, G and Hess.
    DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)           ! Maximum number of parameters is hardcoded here
    DOUBLE PRECISION, INTENT(INOUT) :: ENERGY               ! The energy of the configuration
    DOUBLE PRECISION, INTENT(INOUT) :: GRAD(3*NATOMS)       ! The energy gradient
    LOGICAL, INTENT(IN)           :: GTEST, STEST         ! Flags to specify whether the gradient should be calculated (GTEST) and
                                                          ! whether the Hessian should be calculated (STEST)

    ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: R, R2(NATOMS,NATOMS), R3, Q1, Q2
    DOUBLE PRECISION :: G(NATOMS,NATOMS), F(NATOMS,NATOMS)  ! G tensor and F tensor (see ISOTROPIC_GRAD and ISOTROPIC_HESSIAN, below)
    DOUBLE PRECISION :: TMP_ENERGY
    INTEGER :: J1, J2, J3, J4, J5, J6

    Q1 = PARAMS(1)
    Q2 = PARAMS(2)
    TMP_ENERGY=0.0D0

    ! The arrangement of these IF statements seems slightly odd, but it's been done to make sure we only need to set up one pair
    ! of DO loops for each call.
    IF (GTEST) THEN
        IF (STEST) THEN  ! Gradient + Hessian

            ! Rather than looping over all pairs of atoms, we loop over all interacting pairs as specified by POTLIST.
            ! The first entry in each element of POTLIST contains the index of a principal atom. All other entries contain
            ! the interaction partners for this principal atom.
            DO J6=1,UBOUND(POTLIST,1)
                J1=POTLIST(J6,1)  ! Get the actual atomic index for this atom
                J3=3*J1
                R2(J1,J1)=0.0D0
                G(J1,J1)=0.0D0
                F(J1,J1)=0.0D0

                ! Loop over the partners for atom J1
                DO J5=2,N_ATOM_PARTNERS(J6)+1  ! Have to add the +1 because the first element in POTLIST(J1,:)
                                               ! is the atom index J1, instead of a neighbour index.

                    J2=POTLIST(J6,J5)  ! Get the actual atomic index for this atom
                    ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                    J4=3*J2

                    R2(J2,J1)=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                    R2(J2,J1)=1.0D0/R2(J2,J1)  ! From now on, both R2 and R are inverse lengths
                    R=SQRT(R2(J2,J1))

                    TMP_ENERGY=TMP_ENERGY+Q1*Q2*R

!                    write(*,*) "E+G+H. Pair, Dist, Energy:", J1, J2, 1.0D0/R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    ! Set up storage arrays to use in the gradient- and hessian-calculating routines
                    R2(J1,J2) = R2(J2,J1)
                    R3 = R*R2(J2,J1)

                    G(J2,J1)=-Q1*Q2*R3
                    !G(J2,J1)=-Q(J1)*Q(J2)*R3
                    G(J1,J2)=G(J2,J1)

                    F(J2,J1)=-3*G(J2,J1)
                    F(J1,J2)=F(J2,J1)

                ENDDO
            ENDDO
        ELSE             ! Energy + Gradient only
            DO J6=1,UBOUND(POTLIST,1)
                J1=POTLIST(J6,1)
                J3=3*J1
                G(J1,J1)=0.0D0

                DO J5=2,N_ATOM_PARTNERS(J6)+1

                    J2=POTLIST(J6,J5)  ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                    J4=3*J2

                    R=SQRT((X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2)
                    R=1.0D0/R
                    R3=R**3

                    TMP_ENERGY=TMP_ENERGY+Q1*Q2*R
!                    write(*,*) "E+G. Pair, Dist, Energy:", J1, J2, (X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2, SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    G(J2,J1)=-Q1*Q2*R3
                    !G(J2,J1)=-Q(J1)*Q(J2)*R3
                    G(J1,J2)=G(J2,J1)
                ENDDO
            ENDDO
        ENDIF
    ELSE                ! Energy only
        DO J6=1,UBOUND(POTLIST,1)
            J1=POTLIST(J6,1)
            J3=3*J1

            DO J5=2,N_ATOM_PARTNERS(J6)+1

                J2=POTLIST(J6,J5)  ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                J4=3*J2

                R=SQRT((X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2)
                R=1.0D0/R

                TMP_ENERGY=TMP_ENERGY+Q1*Q2*R
!                write(*,*) "E. Pair, Dist, Energy:", J1, J2, R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
            ENDDO
        ENDDO

    ENDIF

    ENERGY = ENERGY + POTSCALE*TMP_ENERGY

    IF (.NOT.GTEST) RETURN

    ! G should already be set
    CALL EXCLUDE_ISOTROPIC_GRAD(POTLIST, N_ATOM_PARTNERS, X, G, GRAD, POTSCALE)

    IF (.NOT.STEST) RETURN  ! It is assumed we will never need the Hessian without also needing the gradient.

    ! This is currently a very inefficient implementation! I haven't implemented the use of POTLISTS yet.
    CALL EXCLUDE_ISOTROPIC_HESSIAN(POTLIST, N_ATOM_PARTNERS, X, G, F, R2, POTSCALE)

    RETURN

END SUBROUTINE ISO_COULOMB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LJ potential with support for excluded interaction sets. 
! Atom radius sigma is given as a variable, but will often be set to 1.
SUBROUTINE EXCLUDE_ISO_LJ(X, POTLIST, N_ATOM_PARTNERS, POTSCALE, PARAMS, ENERGY, GRAD, GTEST, STEST)

    ! Because we're passing in the entire coordinates array here, we need to know the total number of atoms
    USE COMMONS, ONLY: NATOMS
!    USE MODHESS

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: X(3*NATOMS)          ! The full atomistic coordinate array
    INTEGER, INTENT(IN)           :: POTLIST(:,:)         ! An array containing the interaction partners for each atom
    INTEGER, INTENT(IN)           :: N_ATOM_PARTNERS(:)   ! An array containing the number of interaction partners for each atom
    DOUBLE PRECISION, INTENT(IN)  :: POTSCALE             ! The energy unit for this potential, which multiplies E, G and Hess.
    DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)           ! Maximum number of parameters is hardcoded here
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY               ! The energy of the configuration
    DOUBLE PRECISION, INTENT(OUT) :: GRAD(3*NATOMS)       ! The energy gradient
    LOGICAL, INTENT(IN)           :: GTEST, STEST         ! Flags to specify whether the gradient should be calculated (GTEST) and
                                                          ! whether the Hessian should be calculated (STEST)

    ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: DIST, R2(NATOMS,NATOMS), R6, R8(NATOMS,NATOMS), R14(NATOMS,NATOMS), SIG, SIG6, SIG12
    DOUBLE PRECISION :: G(NATOMS,NATOMS), F(NATOMS,NATOMS)  ! G tensor and F tensor (see ISOTROPIC_GRAD and ISOTROPIC_HESSIAN, below)
    DOUBLE PRECISION :: TMP_ENERGY
    INTEGER :: J1, J2, J3, J4, J5, J6

    SIG = PARAMS(1)
    SIG6 = SIG**6
    SIG12 = SIG6**2

    TMP_ENERGY=0.0D0

!    F(:,:) = 0.0D0
                
    ! The arrangement of these IF statements seems slightly odd, but it's been done to make sure we only need to set up one pair
    ! of DO loops for each call.
    IF (GTEST) THEN
        IF (STEST) THEN  ! Gradient + Hessian

            ! Rather than looping over all pairs of atoms, we loop over all interacting pairs as specified by POTLIST.
            ! The first entry in each element of POTLIST contains the index of a principal atom. All other entries contain
            ! the interaction partners for this principal atom.
            DO J6=1,UBOUND(POTLIST,1)
                J1=POTLIST(J6,1)  ! Get the actual atomic index for this atom
                J3=3*J1
                R2(J1,J1)=0.0D0
                R8(J1,J1)=0.0D0
                R14(J1,J1)=0.0D0
                G(J1,J1)=0.0D0
                F(J1,J1)=0.0D0

                ! Loop over the partners for atom J1
                DO J5=2,N_ATOM_PARTNERS(J6)+1  ! Have to add the +1 because the first element in POTLIST(J1,:)
                                               ! is the atom index J1, instead of a neighbour index.

                    J2=POTLIST(J6,J5)  ! Get the actual atomic index for this atom
                    ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                    J4=3*J2
                   
                    R2(J2,J1)=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                    R2(J2,J1)=1.0D0/R2(J2,J1)
                    R6=R2(J2,J1)**3

                    TMP_ENERGY=TMP_ENERGY+SIG6*R6*(SIG6*R6-1.0D0)

!                    write(*,*) "E+G+H. Pair, Dist, Energy:", J1, J2, 1.0D0/R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    ! Set up storage arrays to use in the gradient- and hessian-calculating routines
                    R8(J2,J1)=R2(J2,J1)**4
                    R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
                    R2(J1,J2)=R2(J2,J1)
                    G(J2,J1)=-24.0D0*(2.0D0*SIG6*R6-1.0D0)*R2(J1,J2)*SIG6*R6
                    G(J1,J2)=G(J2,J1)
                    F(J2,J1)=672.0D0*R14(J2,J1)*SIG12-192.0D0*R8(J2,J1)*SIG6
                    F(J1,J2)=F(J2,J1)
                ENDDO
            ENDDO
        ELSE             ! Energy + Gradient only
            DO J6=1,UBOUND(POTLIST,1)
                J1=POTLIST(J6,1)
                J3=3*J1
                G(J1,J1)=0.0D0

                DO J5=2,N_ATOM_PARTNERS(J6)+1
                   
                    J2=POTLIST(J6,J5)  ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                    J4=3*J2

                    DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                    DIST=1.0D0/DIST
                    R6=DIST**3

                    TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0)
!                    write(*,*) "E+G. Pair, Dist, Energy:", J1, J2, (X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2, SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    G(J2,J1)=-24.0D0*(2.0D0*R6*SIG6-1.0D0)*DIST*R6*SIG6
                    G(J1,J2)=G(J2,J1)
                ENDDO
            ENDDO
        ENDIF
    ELSE                ! Energy only
        DO J6=1,UBOUND(POTLIST,1)
            J1=POTLIST(J6,1)
            J3=3*J1
!            G(J1,J1)=0.0D0

            DO J5=2,N_ATOM_PARTNERS(J6)+1

                J2=POTLIST(J6,J5)  ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                J4=3*J2

                DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
                DIST=1.0D0/DIST
                R6=DIST**3

                TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0)
!                write(*,*) "E. Pair, Dist, Energy:", J1, J2, R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
            ENDDO
        ENDDO

    ENDIF
    TMP_ENERGY=4.0D0*TMP_ENERGY
    ENERGY = ENERGY + POTSCALE*TMP_ENERGY

    IF (.NOT.GTEST) RETURN

    ! G should already be set
    CALL EXCLUDE_ISOTROPIC_GRAD(POTLIST, N_ATOM_PARTNERS, X, G, GRAD, POTSCALE)

    IF (.NOT.STEST) RETURN  ! It is assumed we will never need the Hessian without also needing the gradient.

    ! This is currently a very inefficient implementation! I haven't implemented the use of POTLISTS yet.
    CALL EXCLUDE_ISOTROPIC_HESSIAN(POTLIST, N_ATOM_PARTNERS, X, G, F, R2, POTSCALE)

    RETURN

END SUBROUTINE EXCLUDE_ISO_LJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WCA potential with support for excluded interactions.
! Atom radius sigma is given as a variable (SIG), but will often be set to 1.
SUBROUTINE EXCLUDE_ISO_WCA(X, POTLIST, N_ATOM_PARTNERS, POTSCALE, PARAMS, ENERGY, GRAD, GTEST, STEST)

    ! Because we're passing in the entire coordinates array here, we need to know the total number of atoms
    USE COMMONS, ONLY: NATOMS
!    USE MODHESS

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: X(3*NATOMS)          ! The full atomistic coordinate array
    INTEGER, INTENT(IN)           :: POTLIST(:,:)         ! An array containing the interaction partners for each atom
    INTEGER, INTENT(IN)           :: N_ATOM_PARTNERS(:)   ! An array containing the number of interaction partners for each atom
    DOUBLE PRECISION, INTENT(IN)  :: POTSCALE             ! The energy unit for this potential, which multiplies E, G and Hess.
    DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)           ! Maximum number of parameters is hardcoded here
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY               ! The energy of the configuration
    DOUBLE PRECISION, INTENT(OUT) :: GRAD(3*NATOMS)       ! The energy gradient
    LOGICAL, INTENT(IN)           :: GTEST, STEST         ! Flags to specify whether the gradient should be calculated (GTEST) and
                                                          ! whether the Hessian should be calculated (STEST)

    ! Various powers of the distance between the atoms, and the atom radius
    DOUBLE PRECISION :: DIST, R2(NATOMS,NATOMS), R6, R8(NATOMS,NATOMS), R14(NATOMS,NATOMS), SIG, SIG6, SIG12
    DOUBLE PRECISION :: G(NATOMS,NATOMS), F(NATOMS,NATOMS)  ! G tensor and F tensor (see ISOTROPIC_GRAD and ISOTROPIC_HESSIAN, below)
    DOUBLE PRECISION :: TMP_ENERGY
    INTEGER :: J1, J2, J3, J4, J5, J6

    SIG = PARAMS(1)
    SIG6 = SIG**6
    SIG12 = SIG6**2

    TMP_ENERGY=0.0D0
!    G(:,:)=0.0D0
!    F(:,:)=0.0D0

    ! The arrangement of these IF statements seems slightly odd, but it's been done to make sure we only need to set up one pair
    ! of DO loops for each call.
    IF (GTEST) THEN
        IF (STEST) THEN  ! Gradient + Hessian

            ! Rather than looping over all pairs of atoms, we loop over all interacting pairs as specified by POTLIST.
            ! The first entry in each element of POTLIST contains the index of a principal atom. All other entries contain
            ! the interaction partners for this principal atom.
            DO J6=1,UBOUND(POTLIST,1)
                J1=POTLIST(J6,1)  ! Get the actual atomic index for this atom

                J3=3*J1
                R2(J1,J1)=0.0D0
                R8(J1,J1)=0.0D0
                R14(J1,J1)=0.0D0
                G(J1,J1)=0.0D0
                F(J1,J1)=0.0D0

                ! Loop over the partners for atom J1
                DO J5=2,N_ATOM_PARTNERS(J6)+1  ! Have to add the +1 because the first element in POTLIST(J1,:)
                                               ! is the atom index J1, instead of a neighbour index.

                    J2=POTLIST(J6,J5)  ! Get the actual atomic index for this atom

                    ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                    J4=3*J2

                    R2(J2,J1)=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2

                    IF(R2(J2,J1).GT.WCA_CUT*SIG*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
                    ! We don't compute the energy for this pair of atoms. Set G and F to 0 so that the gradient and
                    ! Hessian terms will go to 0 also.
!                        write(*,*) "E+G+H. Separation>cutoff:", J1, J2, R2(J2,J1)
                        G(J1,J2) = 0.0D0
                        G(J2,J1) = 0.0D0
                        F(J1,J2) = 0.0D0
                        F(J2,J1) = 0.0D0
                        CYCLE
                    ENDIF

                    R2(J2,J1)=1.0D0/R2(J2,J1)
                    R6=R2(J2,J1)**3

                    TMP_ENERGY=TMP_ENERGY+SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
!                    write(*,*) "E+G+H. Pair, Dist, Energy:", J1, J2, 1.0D0/R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    ! Set up storage arrays to use in the gradient- and hessian-calculating routines
                    R8(J2,J1)=R2(J2,J1)**4
                    R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
                    R2(J1,J2)=R2(J2,J1)
                    G(J2,J1)=-24.0D0*(2.0D0*SIG6*R6-1.0D0)*R2(J1,J2)*SIG6*R6
                    G(J1,J2)=G(J2,J1)
                    F(J2,J1)=672.0D0*R14(J2,J1)*SIG12-192.0D0*R8(J2,J1)*SIG6
                    F(J1,J2)=F(J2,J1)
                ENDDO
            ENDDO
        ELSE             ! Energy + Gradient only
            DO J6=1,UBOUND(POTLIST,1)
                J1=POTLIST(J6,1)
                J3=3*J1
                G(J1,J1)=0.0D0

                DO J5=2,N_ATOM_PARTNERS(J6)+1
                   
                    J2=POTLIST(J6,J5)  ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                    J4=3*J2

                    DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2

                    IF(DIST.GT.WCA_CUT*SIG*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
!                        write(*,*) "E+G. Separation>cutoff:", J1, J2, DIST
                        G(J1,J2) = 0.0D0
                        G(J2,J1) = 0.0D0
                        CYCLE
                    ENDIF

                    DIST=1.0D0/DIST
                    R6=DIST**3

                    TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0) + 0.25D0
!                    write(*,*) "E+G. Pair, Dist, Energy:", J1, J2, (X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2, SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
                    G(J2,J1)=-24.0D0*(2.0D0*R6*SIG6-1.0D0)*DIST*R6*SIG6
                    G(J1,J2)=G(J2,J1)
                ENDDO
            ENDDO
        ENDIF
    ELSE                ! Energy only
        DO J6=1,UBOUND(POTLIST,1)
            J1=POTLIST(J6,1)
            J3=3*J1
            G(J1,J1)=0.0D0

            DO J5=2,N_ATOM_PARTNERS(J6)+1

                J2=POTLIST(J6,J5)  ! For historical reasons, the order of the variables J2, J4, J5 is odd.
                J4=3*J2

                DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2

                IF(DIST.GT.WCA_CUT*SIG*SIG) THEN   ! WCA_CUT is a PARAMETER - see top of file.
!                    write(*,*) "E. Separation>cutoff:", J1, J2, R2(J2,J1)
                    CYCLE
                ENDIF

                DIST=1.0D0/DIST
                R6=DIST**3

                TMP_ENERGY=TMP_ENERGY+R6*SIG6*(R6*SIG6-1.0D0) + 0.25D0
!                write(*,*) "E. Pair, Dist, Energy:", J1, J2, R2(J2,J1), SIG6*R6*(SIG6*R6-1.0D0) + 0.25D0
            ENDDO
        ENDDO

    ENDIF
    TMP_ENERGY=4.0D0*TMP_ENERGY
    ENERGY = ENERGY + POTSCALE*TMP_ENERGY

    IF (.NOT.GTEST) RETURN

    ! G should already be set
    CALL EXCLUDE_ISOTROPIC_GRAD(POTLIST, N_ATOM_PARTNERS, X, G, GRAD, POTSCALE)

    IF (.NOT.STEST) RETURN  ! It is assumed we will never need the Hessian without also needing the gradient.

    ! This is currently a very inefficient implementation! I haven't implemented the use of POTLISTS yet.
    CALL EXCLUDE_ISOTROPIC_HESSIAN(POTLIST, N_ATOM_PARTNERS, X, G, F, R2, POTSCALE)

    RETURN

END SUBROUTINE EXCLUDE_ISO_WCA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END OF ISOTROPIC POTENTIALS
! Next, two functions that are general to all isotropic potentials. To be clear, by "isotropic", I mean "depends only
! on the distance between the two particles"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For all isotropic potentials the gradient is calculated in the same way. We simply need to know the positions of the
! two particles, and the G tensor: G = (1/r)dU/dr for an isotropic potential U(r)
SUBROUTINE ISOTROPIC_GRAD(N, X, G, TMP_G)
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: N  ! Number of atoms
    DOUBLE PRECISION, INTENT(IN)  :: X(3*N), G(N,N)
    DOUBLE PRECISION, INTENT(OUT) :: TMP_G(3*N)
    INTEGER :: J1, J2, J3, J4
    DOUBLE PRECISION :: DUMMYX, DUMMYY, DUMMYZ, XMUL

    DO J1=1,N  ! The atom for which we are calculating the gradient
        J3=3*J1
        DUMMYX = 0.0D0
        DUMMYY = 0.0D0
        DUMMYZ = 0.0D0
        DO J4=1,N  ! Consider the interaction with each other atom in turn.
            ! This inner loop is the slow bit. Minimise the number of operations required here.
            J2=3*J4
            XMUL=G(J4,J1)  ! Only do the array-access once to save time

            DUMMYX=DUMMYX+XMUL*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL*(X(J3)-X(J2))
        ENDDO
        TMP_G(J3-2) = DUMMYX
        TMP_G(J3-1) = DUMMYY
        TMP_G(J3) = DUMMYZ
    ENDDO

    RETURN
END SUBROUTINE ISOTROPIC_GRAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For all isotropic potentials the Hessian is calculated in the same way. We simply need to know the positions of the
! two particles, the G tensor (see above) and the F tensor: F = r*d[(1/r)dU/dr]/dr
SUBROUTINE ISOTROPIC_HESSIAN(N, X, G, F, R2, TMP_HESS)
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: N  ! Number of atoms
    DOUBLE PRECISION, INTENT(IN)  :: X(3*N), G(N,N), F(N,N), R2(N,N)
    DOUBLE PRECISION, INTENT(OUT) :: TMP_HESS(3*N,3*N)
    INTEGER :: J1, J2, J3, J4, J5, J6
    DOUBLE PRECISION :: DUMMY

!       Compute the hessian. First are the entirely diagonal terms (3N of them)
!       These correspond to 2nd derivatives wrt the same coordinate twice
    DO J1=1,N
        DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
                DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*(X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)
            ENDDO
            TMP_HESS(J3,J3)=DUMMY
        ENDDO
    ENDDO

!  Next are the terms where x_i and x_j are on the same atom
!  but are different, e.g. y and z.

    DO J1=1,N
        DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
                DUMMY=0.0D0
                DO J5=1,N
                    DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)*(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
                ENDDO
                TMP_HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
        ENDDO
    ENDDO

!  Case III, different atoms, same cartesian coordinate.

    DO J1=1,N
        DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
                TMP_HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*(X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1)
            ENDDO
        ENDDO
    ENDDO

!  Case IV: different atoms and different cartesian coordinates.

    DO J1=1,N
        DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
                DO J5=1,J2-1
                    J6=3*(J4-1)+J5
                    TMP_HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
                    TMP_HESS(3*(J4-1)+J2,3*(J1-1)+J5)=TMP_HESS(J6,J3)
                ENDDO
            ENDDO
        ENDDO
    ENDDO

!  Symmetrise Hessian

    DO J1=1,3*N
        DO J2=J1+1,3*N
            TMP_HESS(J1,J2)=TMP_HESS(J2,J1)
        ENDDO
    ENDDO

    RETURN
END SUBROUTINE ISOTROPIC_HESSIAN


! For all isotropic potentials the gradient is calculated in the same way. We simply need to know the positions of the
! two particles, and the G tensor: G = (1/r)dU/dr for an isotropic potential U(r)
SUBROUTINE EXCLUDE_ISOTROPIC_GRAD(POTLIST, N_ATOM_PARTNERS, X, G, GRAD, POTSCALE)

    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: POTLIST(:,:)
    INTEGER, INTENT(IN)           :: N_ATOM_PARTNERS(:)  ! Number of atoms
    DOUBLE PRECISION, INTENT(IN)  :: X(:), G(:,:)
    DOUBLE PRECISION, INTENT(OUT) :: GRAD(:)  ! This is the real full gradient, so we don't want to zero it!
    DOUBLE PRECISION, INTENT(IN)  :: POTSCALE
    INTEGER :: J1, J2, J3, J4, J5, J6
    DOUBLE PRECISION :: DUMMYX, DUMMYY, DUMMYZ, XMUL

    DO J6=1,UBOUND(POTLIST,1)  ! The atom for which we are calculating the gradient
        J1=POTLIST(J6,1)
        J3=3*J1

        DUMMYX = 0.0D0
        DUMMYY = 0.0D0
        DUMMYZ = 0.0D0

        ! This inner loop is the slow bit. Minimise the number of operations required here.
        DO J5=2,N_ATOM_PARTNERS(J6)+1
            J4=POTLIST(J6,J5)
            J2=3*J4
            XMUL=G(J4,J1)  ! Only do the array-access once to save time

            DUMMYX=XMUL*(X(J3-2)-X(J2-2))*POTSCALE
            DUMMYY=XMUL*(X(J3-1)-X(J2-1))*POTSCALE
            DUMMYZ=XMUL*(X(J3)-X(J2))*POTSCALE

            GRAD(J3-2) = GRAD(J3-2) + DUMMYX
            GRAD(J3-1) = GRAD(J3-1) + DUMMYY
            GRAD(J3) = GRAD(J3) + DUMMYZ
            ! Symmetrise. The G tensor is symmetric (because U(r) depends only on r and r is symmetric wrt exchange of particles)
            ! So the only thing which changes is X1-X2 becomes X2-X1, i.e. a sign change.
            GRAD(J2-2) = GRAD(J2-2) - DUMMYX
            GRAD(J2-1) = GRAD(J2-1) - DUMMYY
            GRAD(J2) = GRAD(J2) - DUMMYZ

        ENDDO
    ENDDO

    RETURN
END SUBROUTINE EXCLUDE_ISOTROPIC_GRAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For all isotropic potentials the Hessian is calculated in the same way. We simply need to know the positions of the
! two particles, the G tensor (see above) and the F tensor: F = r*d[(1/r)dU/dr]/dr
SUBROUTINE EXCLUDE_ISOTROPIC_HESSIAN(POTLIST, N_ATOM_PARTNERS, X, G, F, R2, POTSCALE)
! NOTE! THIS ROUTINE DOES NOT WORK CURRENTLY! USE NUMERICAL HESSIANS INSTEAD!


!    USE COMMONS, ONLY: NATOMS
    USE MODHESS

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: POTLIST(:,:)
    INTEGER, INTENT(IN)           :: N_ATOM_PARTNERS(:)  ! Number of atoms
    DOUBLE PRECISION, INTENT(IN)  :: X(:), G(:,:), F(:,:), R2(:,:), POTSCALE
    INTEGER :: J1, J2, J3, J4, J5, J6, J7, J8
    DOUBLE PRECISION :: DUMMY

    write(*,*) "Hessian for exclusion potentials not currently working. Use numerical hessians instead."
    STOP

!       Compute the hessian. First are the entirely diagonal terms (3N of them)
!       These correspond to 2nd derivatives wrt the same coordinate twice

    DO J6=1,UBOUND(POTLIST,1)
        J1=POTLIST(J6,1)

        DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0

            DO J5=2,N_ATOM_PARTNERS(J6)+1
                J4=POTLIST(J6,J5)

                DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*(X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)
            ENDDO
            HESS(J3,J3)=HESS(J3,J3)+DUMMY*POTSCALE
        ENDDO
    ENDDO

!  Next are the terms where x_i and x_j are on the same atom
!  but are different, e.g. y and z.


    DO J6=1,UBOUND(POTLIST,1)
        J1=POTLIST(J6,1)

        DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
                DUMMY=0.0D0

                DO J7=2,N_ATOM_PARTNERS(J6)+1
                    J5=POTLIST(J6,J7)

                    DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)*(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
                ENDDO
                HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+DUMMY*POTSCALE
                HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)+DUMMY*POTSCALE ! Symmetrise
            ENDDO
        ENDDO
    ENDDO

!  Case III, different atoms, same cartesian coordinate.


    DO J6=1,UBOUND(POTLIST,1)
        J1=POTLIST(J6,1)

        DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0

            DO J5=2,N_ATOM_PARTNERS(J6)+1
                J4=POTLIST(J6,J5)

                DUMMY=-F(J4,J1)*R2(J4,J1)*(X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1)
                HESS(3*(J4-1)+J2,J3)=HESS(3*(J4-1)+J2,J3)+DUMMY*POTSCALE
                HESS(J3,3*(J4-1)+J2)=HESS(J3,3*(J4-1)+J2)+DUMMY*POTSCALE  ! Symmetrise
            ENDDO           
        ENDDO
    ENDDO

!  Case IV: different atoms and different cartesian coordinates.

    DO J7=1,UBOUND(POTLIST,1)
        J1=POTLIST(J7,1)

        DO J2=1,3
            J3=3*(J1-1)+J2

            DO J8=2,N_ATOM_PARTNERS(J7)+1
                J4=POTLIST(J7,J8)

                DO J5=1,J2-1
                    J6=3*(J4-1)+J5
                    DUMMY = -F(J4,J1)*R2(J4,J1)*(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
                    HESS(J6,J3)=HESS(J6,J3)+DUMMY*POTSCALE
                    HESS(J3,J6)=HESS(J3,J6)+DUMMY*POTSCALE ! Symmetrise
                    HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(3*(J4-1)+J2,3*(J1-1)+J5)+DUMMY*POTSCALE
                    HESS(3*(J1-1)+J5,3*(J4-1)+J2)=HESS(3*(J1-1)+J5,3*(J4-1)+J2)+DUMMY*POTSCALE ! Symmetrise
                ENDDO
            ENDDO
        ENDDO
    ENDDO

    RETURN
END SUBROUTINE EXCLUDE_ISOTROPIC_HESSIAN


END MODULE ISOTROPIC_POTENTIALS
