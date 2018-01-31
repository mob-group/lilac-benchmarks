      SUBROUTINE GENRIGIDTIP (R, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, RBSITE, STCHRG, DEBUG
      USE GENRIGID

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8
      INTEGER          :: SITESPERBODY
      DOUBLE PRECISION :: G(3*NATOMS), FRQN(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R6, R12, RSS2, ABSR, DUMMY
      DOUBLE PRECISION :: RI(3), RSS(3), P(3), R(3*NATOMS)
      DOUBLE PRECISION :: DVDR(NATOMS,NATOMS), D2VDR2(NATOMS,NATOMS)
      DOUBLE PRECISION :: C12, C6, CH2O
      LOGICAL          :: GTEST, STEST

      DOUBLE PRECISION :: XR(DEGFREEDOMS), GR(DEGFREEDOMS), HR(DEGFREEDOMS,DEGFREEDOMS)

      D2VDR2(:,:) = 0.D0; DVDR(:,:) = 0.D0

      IF(NSITEPERBODY(1) == 4) THEN
        SITESPERBODY = 4
      ELSE
        WRITE(*,*) "genrigidtip> Error in setup: wrong number of sites in the water molecules"
        STOP
      ENDIF

      IF(.NOT. ALLOCATED(STCHRG)) THEN  ! This should only happen the first time
        ALLOCATE(STCHRG(4))
      ENDIF

      STCHRG(:) = (/0.D0, 0.52D0, 0.52D0, -1.04D0/)
      C6        = 2552.24D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2510.4D3
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (STEST) HESS(:,:) = 0.D0

      ! Loop over pairs of different rigid bodies - exclude intramolecular interactions
      DO J1 = 1, NRIGIDBODY

         J3 = 3*(J1-1)*SITESPERBODY

         DO J2 = J1 + 1, NRIGIDBODY

            J4 = 3*(J2-1)*SITESPERBODY

!     O-O LJ CONTRIBUTION

            J5 = (J1-1)*SITESPERBODY+1   ! Atom index for the first O atom
            J6 = (J2-1)*SITESPERBODY+1
            J7 = J3 + 1                  ! Start position for the coords of the first O atom
            J8 = J4 + 1
            RSS(:) = R(J7:J7+2) - R(J8:J8+2)

            RSS2   = DOT_PRODUCT(RSS(:),RSS(:))
            ABSR   = DSQRT(RSS2)
            R2     = 1.D0/RSS2
            R6     = R2*R2*R2
            R12    = R6*R6

            ENERGY = ENERGY + C12*R12 - C6*R6

            IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R
!               DVDR(J1,J2) =-6.D0*(2.D0*C12*R12 - C6*R6)*R2
!               DVDR(J2,J1) = DVDR(J1,J2)
               DVDR(J5,J6) =-6.D0*(2.D0*C12*R12 - C6*R6)*R2
               DVDR(J6,J5) = DVDR(J5,J6)

!               G(J7:J7+2) = G(J7:J7+2) + DVDR(J1,J2)*RSS(:)
!               G(J8:J8+2) = G(J8:J8+2) - DVDR(J1,J2)*RSS(:)
                G(J7:J7+2) = G(J7:J7+2) + DVDR(J5,J6)*RSS(:)
                G(J8:J8+2) = G(J8:J8+2) - DVDR(J5,J6)*RSS(:)

!               D2VDR2(J1,J2) = (168.D0*C12*R12 - 48.D0*C6*R6)*R2*R2
               D2VDR2(J5,J6) = (168.D0*C12*R12 - 48.D0*C6*R6)*R2*R2
            ENDIF

            IF (STEST) THEN
!                D2VDR2(J2,J1) = D2VDR2(J1,J2)
                D2VDR2(J6,J5) = D2VDR2(J5,J6)
            ENDIF


            DO I = 2, SITESPERBODY

               J7 = 3*SITESPERBODY*(J1-1) + 3*(I-1) + 1
                  ! Offset for this RB   !  Offset for this atom
               J5 = (J1-1)*SITESPERBODY + I  ! Atom index for this atom

               DO J = 2, SITESPERBODY

                  J8 = 3*SITESPERBODY*(J2-1) + 3*(J-1) + 1
                  J6 = (J2-1)*SITESPERBODY + J

                  RSS(:) = R(J7:J7+2) - R(J8:J8+2)
                  RSS2   = DOT_PRODUCT(RSS(:),RSS(:))
                  R2     = 1.D0/RSS2
                  ABSR   = DSQRT(RSS2)

                  ENERGY = ENERGY + CH2O*STCHRG(I)*STCHRG(J)/ABSR

                  IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R
!                     DVDR(J1+I,J2+J) =-CH2O*STCHRG(I)*STCHRG(J)*R2/ABSR
!                     DVDR(J2+J,J1+I) = DVDR(J1+I,J2+J)
                     DVDR(J5,J6) = -CH2O*STCHRG(I)*STCHRG(J)*R2/ABSR
                     DVDR(J6,J5) = DVDR(J5,J6)

!                     G(J7:J7+2)  = G(J7:J7+2) + DVDR(J1+I,J2+J)*RSS(:)
!                     G(J8:J8+2)  = G(J8:J8+2) - DVDR(J1+I,J2+J)*RSS(:)
                     G(J7:J7+2)  = G(J7:J7+2) + DVDR(J5,J6)*RSS(:)
                     G(J8:J8+2)  = G(J8:J8+2) - DVDR(J5,J6)*RSS(:)

!                     D2VDR2(J1+I,J2+J) = 3.D0*CH2O*STCHRG(I)*STCHRG(J)*R2*R2/ABSR
                     D2VDR2(J5,J6) = 3.D0*CH2O*STCHRG(I)*STCHRG(J)*R2*R2/ABSR
                  ENDIF

                  IF (STEST) THEN

!                     D2VDR2(J2+J,J1+I) = D2VDR2(J1+I,J2+J)
                    D2VDR2(J6,J5) = D2VDR2(J5,J6)

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

      ENDDO

      IF (STEST) THEN

         DO J1 = 1, NATOMS

            J3 = 3*J1

            DO J2 = 1, NATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2

               RSS(:) = R(J3-2:J3) - R(J4-2:J4)

!     [1] THREE COMPLETELY DIAGONAL TERMS: SAME ATOM, SAME COORDINATES

!     xi,xi
               HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2(J1,J2)*RSS(1)*RSS(1) + DVDR(J1,J2)
!     yi,yi
               HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2(J1,J2)*RSS(2)*RSS(2) + DVDR(J1,J2)
!     zi,zi
               HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2(J1,J2)*RSS(3)*RSS(3) + DVDR(J1,J2)


!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME ATOM, DIFFERENT COORDINATES

!     xi,yi
               DUMMY           = D2VDR2(J1,J2)*RSS(1)*RSS(2)
               HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
               HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
               DUMMY           = D2VDR2(J1,J2)*RSS(2)*RSS(3)
               HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
               HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
               DUMMY           = D2VDR2(J1,J2)*RSS(3)*RSS(1)
               HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
               HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY


!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT ATOMS, SAME COORDINATE

!     xi,xj
               HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2(J1,J2)*RSS(1)*RSS(1) - DVDR(J1,J2)
!     yi,yj
               HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2(J1,J2)*RSS(2)*RSS(2) - DVDR(J1,J2)
!     zi,zj
               HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2(J1,J2)*RSS(3)*RSS(3) - DVDR(J1,J2)

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT ATOMS, DIFFERENT COORDINATES

!     xi,yj
               DUMMY           = - D2VDR2(J1,J2)*RSS(1)*RSS(2)
               HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
               HESS(J4-1,J3-2) = HESS(J4-1,J3-2) + DUMMY
!     yi,zj
               DUMMY           = - D2VDR2(J1,J2)*RSS(2)*RSS(3)
               HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
               HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
               DUMMY           = - D2VDR2(J1,J2)*RSS(3)*RSS(1)
               HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
               HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY

            ENDDO

         ENDDO

      ENDIF

!write(*,*) "2nd derivative matrix"
!DO J1=1,UBOUND(D2VDR2,1)
!   DO J2=1,UBOUND(D2VDR2,2)
!      write(*,*) J1, J2, D2VDR2(J1,J2)
!   ENDDO
!ENDDO
!   STOP

! sn402: remove
!    IF(DEBUG .AND. STEST) THEN
!      DO J1 = 1,3*NATOMS
!        DO J2 = 1,3*NATOMS
!            IF(ABS(HESS(J1,J2)-HESS(J2,J1)) .GT. 1.0E-7) THEN
!                write(*,*) "genrigidtip> Asymmetric Hessian, coords ", J1, J2
!            ENDIF
!        ENDDO
!      ENDDO
!      write(*,*) "Cartesian HESSIAN:"
!      DO J1=1,3*NATOMS
!          write(*,*) HESS(J1,:)
!     ENDDO
!    ENDIF

      END SUBROUTINE GENRIGIDTIP

!     ---------------------------------------------------------------------------------------------
