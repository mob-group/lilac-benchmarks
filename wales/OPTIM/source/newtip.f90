      SUBROUTINE NEWTIP (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, STCHRG
      USE KEY, ONLY: TIPID  

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), FRQN(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R6, R12, RSS2, ABSR, DUMMY
      DOUBLE PRECISION :: RI(3), RSS(3), P(3), R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DVDR(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), D2VDR2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
! hk286
      DOUBLE PRECISION :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
      DOUBLE PRECISION :: D2RMI10(3,3), D2RMI20(3,3), D2RMI30(3,3), D2RMI120(3,3), D2RMI230(3,3), D2RMI310(3,3)
! hk286
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R1(NRBSITES*NATOMS/2,3), D2R2(NRBSITES*NATOMS/2,3), D2R3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R12(NRBSITES*NATOMS/2,3), D2R23(NRBSITES*NATOMS/2,3), D2R31(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DOTI1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTI2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTI3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      DOUBLE PRECISION :: DOTJ1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTJ2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTJ3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      DOUBLE PRECISION :: C12, C6, CH2O
      LOGICAL          :: GTEST, STEST

      D2VDR2(:,:) = 0.D0; DVDR(:,:) = 0.D0
      DOTI1(:,:) = 0.D0; DOTI2(:,:) = 0.D0; DOTI3(:,:) = 0.D0
      DOTJ1(:,:) = 0.D0; DOTJ2(:,:) = 0.D0; DOTJ3(:,:) = 0.D0

      CALL DEFTIP4(C12, C6, CH2O)
 
      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (STEST) HESS(:,:) = 0.D0

      ! REALNATOMS is the number of TIP4P molecules (=NOPT/6 because 3xCoM+3xAA degrees of freedom per molecule)
      ! NOPT = 3*NATOMS (number of degrees of freedom for an atomic system), so REALNATOMS=NATOMS/2
      REALNATOMS = NATOMS/2
      ! OFFSET is the number of CoM coords.
      OFFSET     = 3*REALNATOMS
 
      ! The potential is defined by pairwise site-site interactions.
      ! In order to compute it, we need to obtain cartesian coordinates for the system.
      ! Do this here.
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3) ! CoM coords for molecule J1
         P  = X(J5-2:J5) ! AA coords for molecule J1

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, STEST)

! hk286
         IF ( RBAANORMALMODET ) THEN
            P = (/0.0D0, 0.0D0, 0.0D0/)
            CALL RMDFAS(P, RMI0, DRMI10, DRMI20, DRMI30, D2RMI10, D2RMI20, D2RMI30, D2RMI120, D2RMI230, D2RMI310, GTEST, STEST)
         ENDIF

         ! Loop over sites (i.e. atoms) in this molecule
         DO J2 = 1, NRBSITES

            ! J4 is an index for this site relative to a complete list of sites in all molecules
            J4        = NRBSITES*(J1-1) + J2
            ! R(J4,:) contains the cartesian coordinates for atom J4
            R(J4,:)   = RI(:) + MATMUL(RMI,RBSITE(J2,:))

            IF (GTEST .OR. STEST) THEN

! hk286
               IF ( RBAANORMALMODET ) THEN
                  DR1(J4,:) = MATMUL(DRMI10,MATMUL(RMI,RBSITE(J2,:)))
                  DR2(J4,:) = MATMUL(DRMI20,MATMUL(RMI,RBSITE(J2,:)))
                  DR3(J4,:) = MATMUL(DRMI30,MATMUL(RMI,RBSITE(J2,:)))
               ELSE
                  DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
                  DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
                  DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))
               ENDIF

            ENDIF

            IF (STEST) THEN

! hk286
               IF ( RBAANORMALMODET ) THEN
                  D2R1(J4,:) = MATMUL(D2RMI10,MATMUL(RMI,RBSITE(J2,:)))
                  D2R2(J4,:) = MATMUL(D2RMI20,MATMUL(RMI,RBSITE(J2,:)))
                  D2R3(J4,:) = MATMUL(D2RMI30,MATMUL(RMI,RBSITE(J2,:)))
                  
                  D2R12(J4,:) = MATMUL(D2RMI120,MATMUL(RMI,RBSITE(J2,:)))
                  D2R23(J4,:) = MATMUL(D2RMI230,MATMUL(RMI,RBSITE(J2,:)))
                  D2R31(J4,:) = MATMUL(D2RMI310,MATMUL(RMI,RBSITE(J2,:)))
               ELSE
                  D2R1(J4,:) = MATMUL(D2RMI1,RBSITE(J2,:))
                  D2R2(J4,:) = MATMUL(D2RMI2,RBSITE(J2,:))
                  D2R3(J4,:) = MATMUL(D2RMI3,RBSITE(J2,:))
                  
                  D2R12(J4,:) = MATMUL(D2RMI12,RBSITE(J2,:))
                  D2R23(J4,:) = MATMUL(D2RMI23,RBSITE(J2,:))
                  D2R31(J4,:) = MATMUL(D2RMI31,RBSITE(J2,:))
               ENDIF

            ENDIF

         ENDDO

      ENDDO

!            write(*,*) "Cartesian Coordinates of the system"
!            do j6=1,UBOUND(R,1)
!               write(*,*) R(j6,:)
!            enddo
!            stop

      ! Now compute the actual potential and derivatives (analogous to the G and F tensors in ljpshift, for example)
      ! Loop over molecules
      DO J1 = 1, REALNATOMS 

         ! Base index for the CoM coords
         J3 = 3*J1
         ! Base index for the AA coords
         J5 = OFFSET + J3

         ! Loop over other molecules - we only compute interactions between different molecules.
         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     O-O LJ CONTRIBUTION

            ! J7 and J8 are site indices for the O atoms on the two different molecules
            J7 = NRBSITES*(J1-1) + 1
            J8 = NRBSITES*(J2-1) + 1
            ! Compute the displacement between the O atoms
            RSS(:) = R(J7,:) - R(J8,:)
            ! Construct the LJ energy term
            RSS2   = DOT_PRODUCT(RSS(:),RSS(:))
            ABSR   = DSQRT(RSS2)
            R2     = 1.D0/RSS2 ! Note, from now on the R2, R6, R12 terms are all inverse lengths
            R6     = R2*R2*R2
            R12    = R6*R6
            ENERGY = ENERGY + C12*R12 - C6*R6

            IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R  (referred to as the G tensor in some other potential functions)
               DVDR(J7,J8) =-6.D0*(2.D0*C12*R12 - C6*R6)*R2
               DVDR(J8,J7) = DVDR(J7,J8)
!               write(*,*) "Atoms ", J7, J8
!               write(*,*) "Derivative ", DVDR(J7, J8)

               DOTI1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J7,:))
               DOTI2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J7,:))
               DOTI3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J7,:))

               DOTJ1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J8,:))
               DOTJ2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J8,:))
               DOTJ3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J8,:))

               ! Calculate the gradient here.
               ! CoM terms are quite simple, just add up the site terms (J7,J8)
               ! that belong to each molecule (J3,J4)
               G(J3-2:J3) = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)

               ! AA terms are more complicated to obtain from cartesian terms.
               G(J5-2) = G(J5-2) + DVDR(J7,J8)*DOTI1(J7,J8)
               G(J5-1) = G(J5-1) + DVDR(J7,J8)*DOTI2(J7,J8)
               G(J5)   = G(J5)   + DVDR(J7,J8)*DOTI3(J7,J8)

               G(J6-2) = G(J6-2) - DVDR(J7,J8)*DOTJ1(J7,J8)
               G(J6-1) = G(J6-1) - DVDR(J7,J8)*DOTJ2(J7,J8)
               G(J6)   = G(J6)   - DVDR(J7,J8)*DOTJ3(J7,J8)

               ! D2VDR2 = 1/r d/dr(1/r dV/dr) = 1/r d/dr(DVDR) is different from the 
               ! F tensor referred to in the LJ potential, for instance.
               D2VDR2(J7,J8) = (168.D0*C12*R12 - 48.D0*C6*R6)*R2*R2
            ENDIF

            IF (STEST) THEN

                D2VDR2(J8,J7) = D2VDR2(J7,J8)
                DOTI1(J8,J7)  =-DOTJ1(J7,J8)
                DOTI2(J8,J7)  =-DOTJ2(J7,J8)
                DOTI3(J8,J7)  =-DOTJ3(J7,J8)
                DOTJ1(J8,J7)  =-DOTI1(J7,J8)
                DOTJ2(J8,J7)  =-DOTI2(J7,J8)
                DOTJ3(J8,J7)  =-DOTI3(J7,J8)

            ENDIF

            ! Now compute the coulomb part of the interaction and derivative tensors.
            ! Because there are multiple coulomb sites on each molecule, we need a loop over these sites.
            DO I = 2, NRBSITES 
               ! Index for this particular interaction site (either an H atom or the dummy charge)
               J7 = NRBSITES*(J1-1) + I

               DO J = 2, NRBSITES 

                  J8     = NRBSITES*(J2-1) + J
                  ! Site-site displacement vector
                  RSS(:) = R(J7,:) - R(J8,:)
                  ! Compute the energy
                  RSS2   = DOT_PRODUCT(RSS(:),RSS(:))
                  R2     = 1.D0/RSS2  ! This is an inverse (square) length
                  ABSR   = DSQRT(RSS2) ! This is a real length
                  ENERGY = ENERGY + CH2O*STCHRG(I)*STCHRG(J)/ABSR

                  ! Calculate the gradient and second derivatives as before.
                  IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R = -CH20 * (q_i * q_j / (r_ij)^3)
                     DVDR(J7,J8) =-CH2O*STCHRG(I)*STCHRG(J)*R2/ABSR
                     DVDR(J8,J7) = DVDR(J7,J8)

                     DOTI1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                     DOTI2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                     DOTI3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J7,:))

                     DOTJ1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J8,:))
                     DOTJ2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J8,:))
                     DOTJ3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J8,:))

                     ! Calculate the gradient here.
                     ! CoM terms are quite simple, just add up the site terms (J7,J8)
                     ! that belong to each molecule (J3,J4)
                     G(J3-2:J3)  = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
                     G(J4-2:J4)  = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)

                     ! AA terms are more complicated to obtain from cartesian terms.
                     G(J5-2)     = G(J5-2) + DVDR(J7,J8)*DOTI1(J7,J8)
                     G(J5-1)     = G(J5-1) + DVDR(J7,J8)*DOTI2(J7,J8)
                     G(J5)       = G(J5)   + DVDR(J7,J8)*DOTI3(J7,J8)

                     G(J6-2)     = G(J6-2) - DVDR(J7,J8)*DOTJ1(J7,J8)
                     G(J6-1)     = G(J6-1) - DVDR(J7,J8)*DOTJ2(J7,J8)
                     G(J6)       = G(J6)   - DVDR(J7,J8)*DOTJ3(J7,J8)

                     ! D2VDR2 = 1/r d/dr(1/r dv/dr)
                     D2VDR2(J7,J8) = 3.D0*CH2O*STCHRG(I)*STCHRG(J)*R2*R2/ABSR
                  ENDIF

                  IF (STEST) THEN

                     D2VDR2(J8,J7) = D2VDR2(J7,J8)
                     DOTI1(J8,J7)  =-DOTJ1(J7,J8)
                     DOTI2(J8,J7)  =-DOTJ2(J7,J8)
                     DOTI3(J8,J7)  =-DOTJ3(J7,J8)
                     DOTJ1(J8,J7)  =-DOTI1(J7,J8)
                     DOTJ2(J8,J7)  =-DOTI2(J7,J8)
                     DOTJ3(J8,J7)  =-DOTI3(J7,J8)

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

      ENDDO

!      IF (.FALSE.) THEN
!         write(*,*) "Calculating Cartesian Hessian for debugging."
!         CALL Cart_hess(R,DVDR,D2VDR2)
!         STOP
!      ENDIF

      ! Now use the derivatives saved in the previous section to compute the Hessian.
      ! Note, the form of the expressions here is slightly different to that used for ljpshift etc., 
      ! because the 2nd derivative tensor has been defined differently.
      IF (STEST) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1  ! Index of the last CoM coord for this molecule
            J5 = OFFSET + J3  ! Index of the last AA coord for this molecule

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4

               DO I = 1, NRBSITES

                  J7 = NRBSITES*(J1 - 1) + I

                  DO J = 1, NRBSITES 

                     J8 = NRBSITES*(J2 - 1) + J

                     RSS(:) = R(J7,:) - R(J8,:) 
                 
!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES

!     xi,xi
                     HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2(J7,J8)*RSS(1)*RSS(1) + DVDR(J7,J8)
!     yi,yi             
                     HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2(J7,J8)*RSS(2)*RSS(2) + DVDR(J7,J8)
!     zi,zi
                     HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2(J7,J8)*RSS(3)*RSS(3) + DVDR(J7,J8)
!     pi1,pi1
                     HESS(J5-2,J5-2) = HESS(J5-2,J5-2) + D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTI1(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR1(J7,:),DR1(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R1(J7,:))
!     pi2,pi2
                     HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTI2(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R2(J7,:))
!     pi3,pi3
                     HESS(J5,J5)     = HESS(J5,J5) + D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTI3(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R3(J7,:))

!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME MOLECULE, DIFFERENT COORDINATES

!     xi,yi
                     DUMMY           = D2VDR2(J7,J8)*RSS(1)*RSS(2)
                     HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
                     HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
                     DUMMY           = D2VDR2(J7,J8)*RSS(2)*RSS(3)
                     HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
                     HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
                     DUMMY           = D2VDR2(J7,J8)*RSS(3)*RSS(1)
                     HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
                     HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
!     xi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(1) + DVDR(J7,J8)*DR1(J7,1)
                     HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
                     HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(2) + DVDR(J7,J8)*DR1(J7,2)
                     HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
                     HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(3) + DVDR(J7,J8)*DR1(J7,3)
                     HESS(J3,J5-2)   = HESS(J3,J5-2) + DUMMY
                     HESS(J5-2,J3)   = HESS(J5-2,J3) + DUMMY
!     xi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(1) + DVDR(J7,J8)*DR2(J7,1)
                     HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
                     HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(2) + DVDR(J7,J8)*DR2(J7,2)
                     HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
                     HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(3) + DVDR(J7,J8)*DR2(J7,3)
                     HESS(J3,J5-1)   = HESS(J3,J5-1) + DUMMY
                     HESS(J5-1,J3)   = HESS(J5-1,J3) + DUMMY
!     xi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(1) + DVDR(J7,J8)*DR3(J7,1)
                     HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
                     HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(2) + DVDR(J7,J8)*DR3(J7,2)
                     HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
                     HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(3) + DVDR(J7,J8)*DR3(J7,3)
                     HESS(J3,J5)     = HESS(J3,J5) + DUMMY
                     HESS(J5,J3)     = HESS(J5,J3) + DUMMY
!     pi1,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTI2(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR2(J7,:),DR1(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R12(J7,:))
                     HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
                     HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTI3(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR3(J7,:),DR2(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R23(J7,:))
                     HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
                     HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTI1(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR1(J7,:),DR3(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R31(J7,:))
                     HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
                     HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY
                 
!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE

!     xi,xj
                     HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2(J7,J8)*RSS(1)*RSS(1) - DVDR(J7,J8)
!     yi,yj
                     HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2(J7,J8)*RSS(2)*RSS(2) - DVDR(J7,J8)
!     zi,zj
                     HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2(J7,J8)*RSS(3)*RSS(3) - DVDR(J7,J8)
!     pi1,pj1
                     HESS(J5-2,J6-2) = HESS(J5-2,J6-2) - D2VDR2(J7,J8)*DOTJ1(J7,J8)*DOTI1(J7,J8) &
                                     - DVDR(J7,J8)*DOT_PRODUCT(DR1(J8,:),DR1(J7,:))
!     pi2,pj2
                     HESS(J5-1,J6-1) = HESS(J5-1,J6-1) - D2VDR2(J7,J8)*DOTJ2(J7,J8)*DOTI2(J7,J8) &
                                     - DVDR(J7,J8)*DOT_PRODUCT(DR2(J8,:),DR2(J7,:))
!     pi3,pj3
                     HESS(J5,J6)     = HESS(J5,J6)     - D2VDR2(J7,J8)*DOTJ3(J7,J8)*DOTI3(J7,J8) &
                                    - DVDR(J7,J8)*DOT_PRODUCT(DR3(J8,:),DR3(J7,:))

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES

!     xi,yj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(1)*RSS(2)
                     HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
                     HESS(J4-1,J3-2) = HESS(J4-1,J3-2) + DUMMY
!     yi,zj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(2)*RSS(3)
                     HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
                     HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(3)*RSS(1)
                     HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
                     HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY


!     xi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(1) - DVDR(J7,J8)*DR1(J8,1)
                     HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
                     HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(2) - DVDR(J7,J8)*DR1(J8,2)
                     HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
                     HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(3) - DVDR(J7,J8)*DR1(J8,3)
                     HESS(J3,J6-2)   = HESS(J3,J6-2) + DUMMY
                     HESS(J6-2,J3)   = HESS(J6-2,J3) + DUMMY
!     xi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(1) - DVDR(J7,J8)*DR2(J8,1)
                     HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
                     HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(2) - DVDR(J7,J8)*DR2(J8,2)
                     HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
                     HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(3) - DVDR(J7,J8)*DR2(J8,3)
                     HESS(J3,J6-1)   = HESS(J3,J6-1) + DUMMY
                     HESS(J6-1,J3)   = HESS(J6-1,J3) + DUMMY
!     xi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(1) - DVDR(J7,J8)*DR3(J8,1)
                     HESS(J3-2,J6)   = HESS(J3-2,J6) + DUMMY
                     HESS(J6,J3-2)   = HESS(J6,J3-2) + DUMMY
!     yi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(2) - DVDR(J7,J8)*DR3(J8,2)
                     HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
                     HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(3) - DVDR(J7,J8)*DR3(J8,3)
                     HESS(J3,J6)     = HESS(J3,J6) + DUMMY
                     HESS(J6,J3)     = HESS(J6,J3) + DUMMY
!     pi1,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTJ2(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR2(J8,:),DR1(J7,:))
                     HESS(J5-2,J6-1) = HESS(J5-2,J6-1) + DUMMY
                     HESS(J6-1,J5-2) = HESS(J6-1,J5-2) + DUMMY
!     pi2,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTJ3(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR3(J8,:),DR2(J7,:))
                     HESS(J5-1,J6)   = HESS(J5-1,J6) + DUMMY
                     HESS(J6,J5-1)   = HESS(J6,J5-1) + DUMMY
!     pi3,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTJ1(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR1(J8,:),DR3(J7,:))
                     HESS(J5,J6-2)   = HESS(J5,J6-2) + DUMMY
                     HESS(J6-2,J5)   = HESS(J6-2,J5) + DUMMY

                  ENDDO

               ENDDO

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

!write(*,*) "Final gradient"   ! sn402
!write(*,*) G
!write(*,*) "Derivative matrix"
!write(*,*) DVDR
!write(*,*) "Hessian"
!!write(*,*) HESS
!DO J1 = 1, 3*NATOMS, 3
!    WRITE(*,*) HESS(J1:J1+2,:)
!    write(*,*) " "
!ENDDO

      END SUBROUTINE NEWTIP 

!     ---------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP4(C12, C6, CH2O)
!     TIP4P water

      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, STCHRG

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.D0

      RBSITE(2,1) = 0.D0
      RBSITE(2,2) = SIN(0.5D0*THETA)*ROH
      RBSITE(2,3) = COS(0.5D0*THETA)*ROH

      RBSITE(3,1) = 0.D0
      RBSITE(3,2) = -SIN(0.5D0*THETA)*ROH
      RBSITE(3,3) = COS(0.5D0*THETA)*ROH

      RBSITE(4,1) = 0.D0
      RBSITE(4,2) = 0.D0
      RBSITE(4,3) = ROM

      STCHRG(:) = (/0.D0, 0.52D0, 0.52D0, -1.04D0/)
      C6        = 2552.24D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2510.4D3
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*RBSITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         RBSITE(I,:) = RBSITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP4

!     ---------------------------------------------------------------------------------------------


      SUBROUTINE INERTIANTIP(RMI, KBLOCK, TMASS)

        USE COMMONS
        USE KEY
        IMPLICIT NONE
        DOUBLE PRECISION :: TMASS, KBLOCK(3,3), RMI(3,3)       
        DOUBLE PRECISION :: MS(3), DR(3)
        INTEGER          :: I, J, J2
        
        MS(:)  = (/16.D0, 1.D0, 1.D0/)
        TMASS  = 18.D0
        CALL DEFTIP4(DC6CC, DC6CC, DC6CC)

        DO J2 = 1, NRBSITES - 1
           DR(:)  = MATMUL(RMI(:,:),RBSITE(J2,:))
           DO I = 1, 3
              KBLOCK(I,I) = KBLOCK(I,I) + MS(J2)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
              DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                 KBLOCK(I,J) = KBLOCK(I,J) - MS(J2)*DR(I)*DR(J)
              ENDDO
           ENDDO
        ENDDO
        
      END SUBROUTINE INERTIANTIP

!     ---------------------------------------------------------------------------------------------

! hk286
      SUBROUTINE VISUALISEMODESTIP4 (AP, X)

        USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
        USE MODHESS
        IMPLICIT NONE
        
        CHARACTER(LEN=20) :: MYFILENAME, MYFILENAME2, ISTR 
        INTEGER :: OFFSET, SMODE, I, J, J1, J2, J3, J5, J9
        DOUBLE PRECISION :: XTEMP(3*NATOMS), XXTEMP(3*NATOMS), X(3*NATOMS), AP(3*NATOMS,3*NATOMS)
        DOUBLE PRECISION :: P(3), KBLOCK(3,3), KBEGNV(3), RMI(3,3), RRMI(3,3), DRMI(3,3), DR(3), MS(3), RBCOORDS(3)
        DOUBLE PRECISION :: TV(NATOMS*9/2), TV2(NATOMS*9/2), TV3(NATOMS*9/2), TV4(NATOMS*9/2), TV5(NATOMS*9/2), TV6(NATOMS*9/2)
        INTEGER, PARAMETER :: LWORK = 10000 ! the dimension is set arbitrarily
        DOUBLE PRECISION :: WORK(LWORK)
        INTEGER          :: INFO, TFRAME
        
        MS(:)  = (/16.D0, 1.D0, 1.D0/)
        OFFSET = 3*NATOMS/2
        TFRAME = 20

! Here are the normal modes you want to be computed - the first six are zero
        DO SMODE = 7, 3*NATOMS
           WRITE (ISTR, '(I10)') SMODE
! Output position of the atoms
           MYFILENAME =trim(adjustl(ISTR))//"VisTIP4.xyz"
! Vectors of the normal mode oscillations
           MYFILENAME2=trim(adjustl(ISTR))//"DrawVector.tcl"
           OPEN(UNIT=26,FILE=trim(adjustl(MYFILENAME )), STATUS="unknown", form="formatted")
           OPEN(UNIT=36,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
           
! 20 frames are used for visualization
           DO J9 = 1, TFRAME
              WRITE(26,'(I6)') (NATOMS/2)*(NRBSITES-1)
              WRITE(26,*) "FRAME ", J9
              DO I = 1, 3*NATOMS
! This gives you DX where DX is taken along the direction of an eigenvector
                 XTEMP(I) = 0.20D0 * AP(I,SMODE) * COS(8.D0*DATAN(1.D0)/TFRAME * J9)
              ENDDO
              DO J1 = 1, NATOMS/2
                 J3 = 3*J1
! For the COM coordinates, this is simple addition X_minimum + DX
                 XXTEMP(J3-2:J3) = X(J3-2:J3) + XTEMP(J3-2:J3)
                 J5 = OFFSET + J3
! Now the rotational part 
                 P  = X(J5-2:J5)
                 KBLOCK(:,:) = 0.D0
                 
! Computes the rotation matrix which brings the reference geometry in stationary frame 
!to the atom positions in current minimum geometry
                 CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)
                 
! Computing inertia matrix in the moving frame, i.e. current minimum geometry
                 DO J2 = 1, NRBSITES - 1
                    DR(:)  = MATMUL(RMI(:,:),RBSITE(J2,:))
                    DO I = 1, 3
                       KBLOCK(I,I) = KBLOCK(I,I) + MS(J2)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
                       DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                          KBLOCK(I,J) = KBLOCK(I,J) - MS(J2)*DR(I)*DR(J)
                       ENDDO
                    ENDDO
                 ENDDO
! Diagonalise inertia matrix
                 CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
                 
! Going from the diagonalised rotation coordinates to per rigid body angle-axis coordinates in the moving frame
                 XXTEMP(J5-2) = KBLOCK(1,1)*XTEMP(J5-2) + KBLOCK(1,2)*XTEMP(J5-1) + KBLOCK(1,3)*XTEMP(J5  )
                 XXTEMP(J5-1) = KBLOCK(2,1)*XTEMP(J5-2) + KBLOCK(2,2)*XTEMP(J5-1) + KBLOCK(2,3)*XTEMP(J5  )
                 XXTEMP(J5  ) = KBLOCK(3,1)*XTEMP(J5-2) + KBLOCK(3,2)*XTEMP(J5-1) + KBLOCK(3,3)*XTEMP(J5  )
                 
! Computing the rotation matrix which rotates the current minimum geometry to its osciallated position, due to normal mode motion
                 P = XXTEMP(J5-2:J5)
                 CALL RMDFAS(P, RRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)
                 
                 DO J2 = 1, NRBSITES - 1
! Position of the atom = shifted position of the COM + position of atom from the centre of mass
                    RBCOORDS(1:3) = XXTEMP(J3-2:J3) + MATMUL(RRMI(:,:),MATMUL(RMI(:,:),RBSITE(J2,:)))
                    IF (J2 == 1) THEN
                       WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ELSE
                       WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ENDIF
                    
! Constructing the visualisation output
                    IF (J9 .EQ. 10) THEN
                       TV( 3*((J1-1)*(NRBSITES-1)+J2)-2:3*((J1-1)*(NRBSITES-1)+J2) ) = RBCOORDS(1:3)
                    ELSEIF (J9 .EQ. 20) THEN
                       TV( 3*((J1-1)*(NRBSITES-1)+J2)-2:3*((J1-1)*(NRBSITES-1)+J2) ) = RBCOORDS(1:3) &
                     - TV( 3*((J1-1)*(NRBSITES-1)+J2)-2:3*((J1-1)*(NRBSITES-1)+J2) )
                    ELSEIF (J9 .EQ. 5) THEN
                       TV2( 3*((J1-1)*(NRBSITES-1)+J2)-2:3*((J1-1)*(NRBSITES-1)+J2) ) = RBCOORDS(1:3)
                    END IF
                    
                 ENDDO
                 
              ENDDO
           ENDDO
           CLOSE (UNIT = 26)
      
           TV(:) = 6*TV(:)
           TV3(:) = TV2(:)-TV(:)/2
           TV4(:) = TV2(:)+TV(:)/2 
           TV5(:) = TV3(:) - TV(:)/3
           TV6(:) = TV4(:) + TV(:)/3
           DO J1 = 1, NATOMS/2
              DO J2 = 1, NRBSITES - 1
                 J9 = ((J1-1)*(NRBSITES-1)+J2)            
                 WRITE(36,'(A15,F6.3,A1,F6.3,A1,F6.3,A3,F6.3,A1,F6.3,A1,F6.3,A15)') "draw cylinder {",TV3(3*J9-2)," &
  &              ",TV3(3*J9-1)," ",TV3(3*J9),"} {",TV4(3*J9-2)," ",TV4(3*J9-1)," ",TV4(3*J9),"} radius 0.10"
                 WRITE(36,'(A15,F6.3,A1,F6.3,A1,F6.3,A3,F6.3,A1,F6.3,A1,F6.3,A15)') "draw cone {",    TV3(3*J9-2)," &
  &              ",TV3(3*J9-1)," ",TV3(3*J9),"} {",TV5(3*J9-2)," ",TV5(3*J9-1)," ",TV5(3*J9),"} radius 0.15"
                 WRITE(36,'(A15,F6.3,A1,F6.3,A1,F6.3,A3,F6.3,A1,F6.3,A1,F6.3,A15)') "draw cone {",    TV4(3*J9-2)," &
  &              ",TV4(3*J9-1)," ",TV4(3*J9),"} {",TV6(3*J9-2)," ",TV6(3*J9-1)," ",TV6(3*J9),"} radius 0.15"
              ENDDO
           ENDDO
           CLOSE (UNIT = 36)
           
        ENDDO
        
      END SUBROUTINE VISUALISEMODESTIP4


      subroutine cart_hess(R, G, F)
        ! Calculate (and print) the Hessian in cartesian coordinates using the tools developed for MULTIPOT. 
        ! This should only be used for debugging!

        USE ISOTROPIC_POTENTIALS, ONLY: ISOTROPIC_HESSIAN
        USE COMMONS, ONLY: NATOMS, NRBSITES

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: R(NRBSITES*NATOMS/2,3) ! Cartesian coordinates of the system
        ! Cartesian G and F tensors (defined in the body of the code, given above)
        DOUBLE PRECISION, INTENT(IN) :: G(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), F(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
        DOUBLE PRECISION :: X(3*NRBSITES*NATOMS/2)
        DOUBLE PRECISION :: R2DUMMY(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
        DOUBLE PRECISION :: TMP_HESS(3*NRBSITES*NATOMS/2,3*NRBSITES*NATOMS/2)
        INTEGER :: J1, J2
        write(*,*) "In cart_hess"

        DO J1=1,NRBSITES*NATOMS/2
           X(3*J1-2:3*J1)=R(J1,:)
        ENDDO
        R2DUMMY(:,:)=1.0D0

        CALL ISOTROPIC_HESSIAN(NRBSITES*NATOMS/2,X,G,F,R2DUMMY,TMP_HESS)

        DO J1=1,3*NRBSITES*NATOMS/2
           DO J2=1,3*NRBSITES*NATOMS/2
              write(*,*) J1, J2, TMP_HESS(J1,J2)
           ENDDO
        ENDDO

     end subroutine cart_hess
