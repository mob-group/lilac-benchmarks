! -----------------------------------------------------------------------------
! dj337
! Anisotropic potential for periodic benzene systems. Potential from:
! T.S.Totton, A.J.Misquitta, M.Kraft, J. Chem. Theory Comput. 6 (2010) 683–695.

! Long-range electrostatic interactions are computed using Ewald summation.

! Implemented within the GENRIGID framework for rigid bodies.

! Gradients wrt lattice parameters are computed when BOXDERIV keyword used. This
! subroutine is for triclinic systems. Use BENZGENRIGIDEWALD_ORTHO for orthorhombic
! cells.
! -----------------------------------------------------------------------------

      SUBROUTINE BENZGENRIGIDEWALD(X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NCARBON, RBSTLA, RHOCC0, RHOCC10, RHOCC20, &
     &                   RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, &
     &                   RHOCH20, ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, &
     &                   EWALDREALC, BOX_PARAMS, BOX_PARAMSGRAD

      ! adapted to the genrigid framework
      USE GENRIGID, ONLY: NRIGIDBODY, ATOMRIGIDCOORDT, TRANSFORMCTORIGID, NSITEPERBODY, &
     &                    MAXSITE, SITESRIGIDBODY, TRANSFORMRIGIDTOC, TRANSFORMGRAD, INVERSEMATRIX

      ! use Ewald summation to compute electrostatics
      USE EWALD
      ! fractional coordinates, triclinic cells, box derivatives
      USE CARTDIST
      USE BOX_DERIVATIVES

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET, FCT(6), L, M, N, IDX
      INTEGER          :: NEWALDREAL(3)
      DOUBLE PRECISION :: X(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: G(3*NATOMS)
      DOUBLE PRECISION :: XR(3*NATOMS), XC(3*NATOMS), G3C(3*NATOMS), G3(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: ENERGY
      DOUBLE PRECISION :: R2, R6, ABSRIJ, DVDR, ENERGY1, ENERGY2, ENERGY3
      DOUBLE PRECISION :: DMPFCT_SHIFT, EXPFCT_SHIFT, VSHIFT1, VSHIFT2, EWALDREALC2
      DOUBLE PRECISION :: RI(3), RR(3), RSS(3), NR(3), P(3), EI(3), EJ(3), FRIJ(3), TIJ(3), TJI(3) 
      DOUBLE PRECISION :: R(MAXSITE*NRIGIDBODY,3), E(3*MAXSITE*NRIGIDBODY,3), xdum(3*natoms)
      DOUBLE PRECISION :: DR1(MAXSITE*NRIGIDBODY,3), DR2(MAXSITE*NRIGIDBODY,3), DR3(MAXSITE*NRIGIDBODY,3)
      DOUBLE PRECISION :: DE1(3*MAXSITE*NRIGIDBODY,3), DE2(3*MAXSITE*NRIGIDBODY,3), DE3(3*MAXSITE*NRIGIDBODY,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DCADR(3), DCBDR(3)
      DOUBLE PRECISION :: RHOCC, RHOHH, RHOCH, COSTA, COSTB, DMPFCT, DDMPDR, EXPFCT 
      DOUBLE PRECISION :: DRIJDPI(3), DRIJDPJ(3), DCADPI(3), DCBDPI(3), DCADPJ(3), DCBDPJ(3)
      DOUBLE PRECISION :: H(3,3), H_grad(3,3,6), H_inverse(3,3), rrfrac(3), rssfracmin(3), rrcom(3), rcomfrac(3)
      double precision :: rrcomfrac(3), rcomfracmin(3), rssfrac(3), vol, v_fact, dv_fact(3), c(3), s(3)
      double precision :: reciplatvec(3,3), reciplatvec_grad(3,3,6)
      DOUBLE PRECISION, PARAMETER :: B = 1.6485D0
      double precision, parameter :: pi = 3.141592654d0
      integer, parameter          :: image_cutoff = 5
      LOGICAL          :: GTEST, keep_angles

      ! check if combination of triclinic cell angles is realistic
      if (boxderivt) then
         keep_angles = check_angles(box_params(4:6))
         if (.not.keep_angles) then
            ! reject the structure if unrealistic cell angles
            call reject(energy, g)
            return
         endif
      endif

      call build_H(H, H_grad, gtest) ! H matrix
      call inversematrix(H, H_inverse) ! inverse of H matrix
      call get_volume(vol) ! cell volume
      call get_reciplatvec(reciplatvec,reciplatvec_grad, .false.) ! reciprocal lattice vectors

      ! figure out how many lattice vectors to sum over
      newaldreal(1) = floor(ewaldrealc*dsqrt(sum(reciplatvec(1,:)**2))/(2.0d0*pi) + 0.5d0)
      newaldreal(2) = floor(ewaldrealc*dsqrt(sum(reciplatvec(2,:)**2))/(2.0d0*pi) + 0.5d0)
      newaldreal(3) = floor(ewaldrealc*dsqrt(sum(reciplatvec(3,:)**2))/(2.0d0*pi) + 0.5d0)

      ! reject structure if would have to sum over more than five lattice vectors
      ! NOTE: this is because cells have tendency to flatten out... it takes a long time
      ! to compute these structures and there are probably equivalent ones of a more
      ! regular shape
      if (boxderivt) then
         if (.not. all(newaldreal.le.image_cutoff)) then
            call reject(energy, g)
            return
         endif
      endif

      ! factorials
      FCT(1) = 1; FCT(2) = 2; FCT(3) = 6; FCT(4) = 24; FCT(5) = 120; FCT(6) = 720
      ! initialize energy values
      ! energy1 is due to short-range anisotropic interactions
      ! energy2 is due to damped dispersion
      ! energy3 is due to long-range electrostatics (computed using Ewald)
      ENERGY = 0.D0; ENERGY1 = 0.D0; ENERGY2 = 0.D0; ENERGY3 = 0.D0

      ! initialize gradient if GTEST true
      IF (GTEST) G(:) = 0.D0
      IF (GTEST) G3C(:) = 0.D0

      ! dj337: check if input coordinates are cartesian
      ! assumes ATOMRIGIDCOORDT is correct
      IF (ATOMRIGIDCOORDT) THEN ! if input is cartesian
         ! convert to rigidbody coordinates
         XR(:) = 0.D0
         CALL TRANSFORMCTORIGID(X, XR)
         if (boxderivt) then
            call frac2cart_rb_tri(nrigidbody, xdum, xr, H)
            x(:) = xdum(:)
         else
            x(:) = xr(:)
         endif
      ENDIF

      if (boxderivt) then
         ! compute v factor
         c(:) = dcos(box_params(4:6))
         s(:) = dsin(box_params(4:6))
         ! v_fact is factor related to volume
         v_fact = dsqrt(1.0d0 - c(1)**2-c(2)**2-c(3)**2 + 2.0d0*c(1)*c(2)*c(3))
         dv_fact(1) = s(1)*(c(1) - c(2)*c(3))/v_fact
         dv_fact(2) = s(2)*(c(2) - c(1)*c(3))/v_fact
         dv_fact(3) = s(3)*(c(3) - c(1)*c(2))/v_fact
      endif

      EWALDREALC2 = EWALDREALC**2 ! real-space cutoff

      ! OFFSET is number of CoM coords (3*NRIGIDBODY)
      OFFSET     = 3*NRIGIDBODY

      ! Computing Cartesian coordinates for the system.  
      DO J1 = 1, NRIGIDBODY

         J3 = 3*J1
         J5 = OFFSET + J3
         ! CoM coords for rigid body J1
         RI = X(J3-2:J3)
         ! AA coords for rigid body J1
         P  = X(J5-2:J5)

         ! calculates rotation matrix (RMI)
         ! also calculates derivatives if GTEST is true
         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         ! loop over sites in the rigid body
         DO J2 = 1, NSITEPERBODY(J1)

            ! J4 is index for site J2 relative to a complete list of all sites in all rigid bodies
            ! dj337: assumes that same number of sites per rigid body (i.e. NSITEPERBODY(J1) == MAXSITE)
            J4      = MAXSITE*(J1-1) + J2
            ! R(J4,:) contains Cartesian coordinates for site J4
            R(J4,:) = RI(:) + MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
            ! E(J4,:) contains Z-axis in local axis system for site J4 
            E(J4,:) = MATMUL(RMI(:,:),RBSTLA(J2,:))

            IF (GTEST) THEN

               ! calculate derivative wrt coordinates
               DR1(J4,:) = MATMUL(DRMI1(:,:),SITESRIGIDBODY(J2,:,J1))
               DR2(J4,:) = MATMUL(DRMI2(:,:),SITESRIGIDBODY(J2,:,J1))
               DR3(J4,:) = MATMUL(DRMI3(:,:),SITESRIGIDBODY(J2,:,J1))

               ! calculate derivative wrt local axis
               DE1(J4,:) = MATMUL(DRMI1(:,:),RBSTLA(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBSTLA(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBSTLA(J2,:))

            ENDIF ! gtest

         ENDDO ! sites

      ENDDO ! rigid bodies

      ! Now compute the actual potential.
      ! loop over rigid bodies (A)
      DO J1 = 1, NRIGIDBODY - 1

         J3 = 3*J1
         J5 = OFFSET + J3
         ! CoM coords for rigid body J1
         RI(:)  = X(J3-2:J3)

         ! loop over sites in the rigid body J1
         DO I = 1, NSITEPERBODY(J1)

            ! J7 is index for site I
            J7    = MAXSITE*(J1-1) + I
            ! EI is Z-axis for site I
            EI(:) = E(J7,:)

            ! loop over rigid bodies (B)   
            DO J2 = J1 + 1, NRIGIDBODY

               J4 = 3*J2
               J6 = OFFSET + J4

               ! loop over sites in the rigid body J2
               DO J = 1, NSITEPERBODY(J2)

                  ! J8 is index for site J
                  J8     = MAXSITE*(J2-1) + J
                  ! EJ is Z-axis for site J
                  EJ(:)  = E(J8,:)
                  rr(:) = r(j7,:) - r(j8,:)
                  ! convert to fractional coordinates
                  rrfrac(:) = matmul(H_inverse, rr(:))
                  ! minimum image convention
                  rssfracmin(1) = rrfrac(1) - anint(rrfrac(1))
                  rssfracmin(2) = rrfrac(2) - anint(rrfrac(2))
                  rssfracmin(3) = rrfrac(3) - anint(rrfrac(3))

                  if (gtest.and.boxderivt) then
                     ! get center of mass separation vector
                     rrcom(:) = x(j3-2:j3) - x(j4-2:j4)
                     ! convert to fractional coordinates
                     rrcomfrac(:) = matmul(H_inverse, rrcom(:))
                     ! minimum image convention
                     rcomfracmin(1) = rrcomfrac(1) - anint(rrfrac(1))
                     rcomfracmin(2) = rrcomfrac(2) - anint(rrfrac(2))
                     rcomfracmin(3) = rrcomfrac(3) - anint(rrfrac(3))
                  endif

                  ! sum over lattice vectors
                  do l = -newaldreal(1), newaldreal(1)
                  rssfrac(1) = rssfracmin(1) + l

                     do m = -newaldreal(2), newaldreal(2)
                     rssfrac(2) = rssfracmin(2) + m

                        do n = -newaldreal(3), newaldreal(3)
                        rssfrac(3) = rssfracmin(3) + n

                        ! convert to absolute coordinates
                        rss(:) = matmul(H, rssfrac(:))

                        ! get COM vector
                        if (gtest.and.boxderivt) then
                           rcomfrac(1) = rcomfracmin(1) + l
                           rcomfrac(2) = rcomfracmin(2) + m
                           rcomfrac(3) = rcomfracmin(3) + n
                        endif
                     
                        R2     = DOT_PRODUCT(RSS(:),RSS(:))
                        ! check if distance within cutoff
                        IF (R2 < EWALDREALC2) THEN
                           ! ABSRIJ is site-site separation between I and J
                           ABSRIJ = DSQRT(R2)
                           ! NR is unit site-site vector from sites I to J
                           NR(:)  = RSS(:)/ABSRIJ
                           R2     = 1.D0/R2
                           R6     = R2*R2*R2
         
      !     CALCULATE THE DISPERSION DAMPING FACTOR
         
                           ! initialize sum for the damping function and vertical shift
                           DMPFCT = 1.D0
                           DMPFCT_SHIFT = 1.D0
                           ! initialize sum for the derivative of damping function
                           DDMPDR = B
         
                           ! calculate sums
                           DO K = 1, 6
         
                              DMPFCT = DMPFCT + (B*ABSRIJ)**K/FLOAT(FCT(K))
                              DMPFCT_SHIFT = DMPFCT_SHIFT + (B*EWALDREALC)**K/FLOAT(FCT(K))
                              IF (K > 1) DDMPDR = DDMPDR + (B**K)*(ABSRIJ)**(K-1)/FLOAT(FCT(K-1))
         
                           END DO
         
                           EXPFCT = DEXP(-B*ABSRIJ)
                           EXPFCT_SHIFT = DEXP(-B*EWALDREALC)
                           ! DDMPDR is derivative of damping function with factor 1/Rab
                           DDMPDR = (B*EXPFCT*DMPFCT - EXPFCT*DDMPDR)/ABSRIJ
                           ! DMPFCT is damping function
                           DMPFCT = 1.D0 - EXPFCT*DMPFCT
                           ! DMPFCT_SHIFT is vertical shift for damping function
                           DMPFCT_SHIFT = 1.D0 - EXPFCT_SHIFT*DMPFCT_SHIFT
         
      !     NOW CALCULATE RHOAB
         
                           ! calculate cos(theta) 
                           COSTA      =-DOT_PRODUCT(NR(:),EI(:))
                           COSTB      = DOT_PRODUCT(NR(:),EJ(:))
         
                           ! calculate terms relevant to derivatives
                           IF (GTEST) THEN
         
                              ! derivative of cos(theta) wrt r_ij
                              DCADR(:)   =-EI(:)/ABSRIJ - COSTA*R2*RSS(:)
                              DCBDR(:)   = EJ(:)/ABSRIJ - COSTB*R2*RSS(:)
         
                              ! derivative of r_ij wrt pi
                              DRIJDPI(1) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                              DRIJDPI(2) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                              DRIJDPI(3) = DOT_PRODUCT(RSS(:),DR3(J7,:))
         
                              ! derivative of r_ij wrt pj
                              DRIJDPJ(1) =-DOT_PRODUCT(RSS(:),DR1(J8,:))
                              DRIJDPJ(2) =-DOT_PRODUCT(RSS(:),DR2(J8,:))
                              DRIJDPJ(3) =-DOT_PRODUCT(RSS(:),DR3(J8,:))
         
                              ! derivative of cos(theta) wrt pi
                              DCADPI(1)  =-DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE1(J7,:)) & 
                                         - COSTA*R2*DRIJDPI(1)
                              DCADPI(2)  =-DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE2(J7,:)) &
                                         - COSTA*R2*DRIJDPI(2)
                              DCADPI(3)  =-DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE3(J7,:)) &
                                         - COSTA*R2*DRIJDPI(3)
                              DCBDPI(1)  = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(1)
                              DCBDPI(2)  = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(2)
                              DCBDPI(3)  = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(3)
                         
                              ! derivative of cos(theta) wrt pj
                              DCADPJ(1)  = DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(1)
                              DCADPJ(2)  = DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(2)
                              DCADPJ(3)  = DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(3)
         
                              DCBDPJ(1)  =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE1(J8,:)) &
                                         - COSTB*R2*DRIJDPJ(1)
                              DCBDPJ(2)  =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE2(J8,:)) &
                                         - COSTB*R2*DRIJDPJ(2)
                              DCBDPJ(3)  =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE3(J8,:)) &
                                         - COSTB*R2*DRIJDPJ(3)
         
                           ENDIF
           
                           ! calculate if I and J are both carbons 
                           IF (I <= NCARBON .AND. J <= NCARBON) THEN
         
                              ! calculate rho_cc
                              RHOCC   = RHOCC0 + RHOCC10*(COSTA + COSTB) + RHOCC20*(1.5D0*COSTA*COSTA & 
                                      + 1.5D0*COSTB*COSTB - 1.D0)
                              ! ENERGY1 is energy due to short-range anisotropic interactions
                              ! calculate vertical shift for first term
                              EXPFCT  = KKJ*DEXP(-ALPHACC*(ABSRIJ - RHOCC))
                              VSHIFT1 = KKJ*DEXP(-ALPHACC*(EWALDREALC - RHOCC))
                              ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                              ! ENERGY2 is energy due to damped dispersion
                              ! calculate vertical shift for second term
                              VSHIFT2 = DC6CC*DMPFCT_SHIFT/(EWALDREALC**6)
                              ENERGY2 = ENERGY2 - DC6CC*DMPFCT*R6 + VSHIFT2
         
                              IF (GTEST) THEN
         
                                 ! DVDR is derivative of dispersion damping factor energy with factor 1/Rab
                                 DVDR    = 6.D0*DC6CC*R6*R2*DMPFCT - DC6CC*R6*DDMPDR 
                                 ! FRIJ is derivative of ENERGY1 wrt r_ij with factor 1/Rab
                                 FRIJ(:) = ALPHACC*EXPFCT*(-NR(:) + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADR(:) &
                                         + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDR(:))
                                 ! TIJ is derivative of ENERGY1 wrt pi with factor 1/Rab
                                 TIJ(:)  = ALPHACC*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADPI(:) &
                                         + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPI(:))
                                 ! TJI is derivative of ENERGY1 wrt pj with factor 1/Rab
                                 TJI(:)  = ALPHACC*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADPJ(:) &
                                         + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPJ(:)) 
         
                              ENDIF
         
                           ! calculate if I and J are both hydorgens
                           ELSEIF (I > NCARBON .AND. J > NCARBON) THEN
         
                              RHOHH  = RHOHH0 + RHOHH10*(COSTA + COSTB) + RHOHH20*(1.5D0*COSTA*COSTA      &
                                     + 1.5D0*COSTB*COSTB - 1.D0) 
                              EXPFCT  = KKJ*DEXP(-ALPHAHH*(ABSRIJ - RHOHH))
                              VSHIFT1 = KKJ*DEXP(-ALPHAHH*(EWALDREALC - RHOHH))
                              ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                              VSHIFT2 = DC6HH*DMPFCT_SHIFT/(EWALDREALC**6)
                              ENERGY2 = ENERGY2 - DC6HH*DMPFCT*R6 + VSHIFT2
         
                              IF (GTEST) THEN
         
                                 DVDR    = 6.D0*DC6HH*R6*R2*DMPFCT - DC6HH*R6*DDMPDR 
                                 FRIJ(:) = ALPHAHH*EXPFCT*(-NR(:) + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADR(:) &
                                         + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDR(:))
                                 TIJ(:)  = ALPHAHH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADPI(:) &
                                         + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPI(:))
                                 TJI(:)  = ALPHAHH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADPJ(:) &
                                         + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPJ(:))
         
                              ENDIF
         
                           ! calculate if I is carbon and J is hydrogen
                           ELSE IF (I <= NCARBON .AND. J > NCARBON) THEN 
         
                              RHOCH  = RHOCH0 + RHOC10H*COSTA + RHOCH10*COSTB + RHOC20H*(1.5D0*COSTA*COSTA &
                                     - 0.5D0) + RHOCH20*(1.5D0*COSTB*COSTB - 0.5D0)
                              EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                              VSHIFT1 = KKJ*DEXP(-ALPHACH*(EWALDREALC - RHOCH))
                              ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                              VSHIFT2 = DC6CH*DMPFCT_SHIFT/(EWALDREALC**6)
                              ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6 + VSHIFT2
         
                              IF (GTEST) THEN
                           
                                 DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                                 FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADR(:) &
                                         + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDR(:))
                                 TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADPI(:) &
                                         + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPI(:))
                                 TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADPJ(:) &
                                         + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPJ(:))
         
                              ENDIF
         
                           ELSE !IF(I > NCARBON .AND. J <= NCARBON) THEN
         
                              RHOCH  = RHOCH0 + RHOCH10*COSTA + RHOC10H*COSTB + RHOCH20*(1.5D0*COSTA*COSTA &
                                     - 0.5D0) + RHOC20H*(1.5D0*COSTB*COSTB - 0.5D0)
                              EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                              VSHIFT1 = KKJ*DEXP(-ALPHACH*(EWALDREALC - RHOCH))
                              ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                              VSHIFT2 = DC6CH*DMPFCT_SHIFT/(EWALDREALC**6)
                              ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6 + VSHIFT2
         
                              IF (GTEST) THEN
         
                                 DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                                 FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADR(:) &
                                         + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDR(:))
                                 TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADPI(:) &
                                         + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPI(:))
                                 TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADPJ(:) &
                                         + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPJ(:))
         
                              ENDIF
         
                           ENDIF
         
                           IF (GTEST) THEN
         
                              ! total gradient wrt CoM coords for rigid body J1
                              G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:) + FRIJ(:)
                              ! total gradient wrt CoM coords for rigid body J2
                              G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:) - FRIJ(:)

                              ! total gradient wrt AA coords for rigid body J1
                              G(J5-2:J5) = G(J5-2:J5) + DVDR*DRIJDPI(:) + TIJ(:)
                              ! total gradient wrt AA coords for rigid body J2
                              G(J6-2:J6) = G(J6-2:J6) + DVDR*DRIJDPJ(:) + TJI(:)

                              ! gradients wrt cell parameters
                              if (boxderivt) then
                                 do idx = 1, 6
                                    box_paramsgrad(idx) = box_paramsgrad(idx) + dot_product((dvdr*rss(1:3) + frij(1:3)), matmul(H_grad(:,:,idx), rcomfrac(:)))
                                 enddo
                              endif ! box derivatives

                           ENDIF ! gtest
      
                        ENDIF ! within cutoff

                     enddo ! n
                  enddo ! m
               enddo ! l

               ENDDO ! sites j

            ENDDO ! rigid bodies J
 
         ENDDO ! sites i

      ENDDO ! rigid bodies I

! INCLUDE CONTRIBUTION OF RIGID BODY WITH PERIODIC IMAGE OF ITSELF

      ! loop over rigidbodies
      do j1 = 1, nrigidbody
         j3 = 3*j1
         j5 = offset + j3
         ri(:) = x(j3-2:j3)

         ! loop over sites i
         do i = 1, nsiteperbody(j1)
            j7 = maxsite*(j1-1) + i
            ei(:) = e(j7,:)

            ! loop over sites j
            do j = 1, nsiteperbody(j1)
               j8 = maxsite*(j1-1) + j
               ej(:) = e(j8,:)

               ! site-site separation vector
               rr(:) = r(j7,:) - r(j8,:)
               ! convert to fractional
               rrfrac(:) = matmul(H_inverse, rr(:))

               ! sum over lattice vectors
               do l = -newaldreal(1), newaldreal(1)
                  do m = -newaldreal(2), newaldreal(2)
                     do n = -newaldreal(3), newaldreal(3)

                     ! if not on same rigid body
                     if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then

                        rssfrac(1) = rrfrac(1) + l
                        rssfrac(2) = rrfrac(2) + m
                        rssfrac(3) = rrfrac(3) + n
                        ! convert back to absolute
                        rss(:) = matmul(H, rssfrac(:))

                        ! get vector between COM
                        if (gtest.and.boxderivt) then
                           rcomfrac(1) = l
                           rcomfrac(2) = m
                           rcomfrac(3) = n
                        endif

                        r2 = dot_product(rss(:), rss(:))
                        if (r2 < ewaldrealc2) then

                        ! absolute site-site distance
                        absrij = dsqrt(r2)
                        nr(:) = rss(:)/absrij
                        r2 = 1.d0/r2
                        r6 = r2*r2*r2

                        ! CALCULATE DISPERSION DAMPING FACTOR

                        ! initialize sum for the damping function and vertical shift
                        DMPFCT = 1.D0
                        DMPFCT_SHIFT = 1.D0
                        ! initialize sum for the derivative of damping function
                        DDMPDR = B

                        ! calculate sums
                        DO K = 1, 6

                           DMPFCT = DMPFCT + (B*ABSRIJ)**K/FLOAT(FCT(K))
                           DMPFCT_SHIFT = DMPFCT_SHIFT + (B*EWALDREALC)**K/FLOAT(FCT(K))
                           IF (K > 1) DDMPDR = DDMPDR + (B**K)*(ABSRIJ)**(K-1)/FLOAT(FCT(K-1))

                        END DO

                        EXPFCT = DEXP(-B*ABSRIJ)
                        EXPFCT_SHIFT = DEXP(-B*EWALDREALC)
                        ! DDMPDR is derivative of damping function with factor 1/Rab
                        DDMPDR = (B*EXPFCT*DMPFCT - EXPFCT*DDMPDR)/ABSRIJ
                        ! DMPFCT is damping function
                        DMPFCT = 1.D0 - EXPFCT*DMPFCT
                        ! DMPFCT_SHIFT is vertical shift for damping function
                        DMPFCT_SHIFT = 1.D0 - EXPFCT_SHIFT*DMPFCT_SHIFT

                        ! CALCULATE RHOAB
                        ! calculate cos(theta) 
                        COSTA      =-DOT_PRODUCT(NR(:),EI(:))
                        COSTB      = DOT_PRODUCT(NR(:),EJ(:))

                        ! calculate terms relevant to derivatives
                        IF (GTEST) THEN

                           ! derivative of cos(theta) wrt r_ij
                           DCADR(:)   =-EI(:)/ABSRIJ - COSTA*R2*RSS(:)
                           DCBDR(:)   = EJ(:)/ABSRIJ - COSTB*R2*RSS(:)

                           ! derivative of r_ij wrt pi
                           DRIJDPI(1) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                           DRIJDPI(2) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                           DRIJDPI(3) = DOT_PRODUCT(RSS(:),DR3(J7,:))

                           ! derivative of r_ij wrt pj
                           DRIJDPJ(1) =-DOT_PRODUCT(RSS(:),DR1(J8,:))
                           DRIJDPJ(2) =-DOT_PRODUCT(RSS(:),DR2(J8,:))
                           DRIJDPJ(3) =-DOT_PRODUCT(RSS(:),DR3(J8,:))

                           ! derivative of cos(theta) wrt pi
                           DCADPI(1)  =-DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE1(J7,:)) & 
                                      - COSTA*R2*DRIJDPI(1)
                           DCADPI(2)  =-DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE2(J7,:)) &
                                      - COSTA*R2*DRIJDPI(2)
                           DCADPI(3)  =-DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE3(J7,:)) &
                                      - COSTA*R2*DRIJDPI(3)
                           DCBDPI(1)  = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(1)
                           DCBDPI(2)  = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(2)
                           DCBDPI(3)  = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(3)

                           ! derivative of cos(theta) wrt pj
                           DCADPJ(1)  = DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(1)
                           DCADPJ(2)  = DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(2)
                           DCADPJ(3)  = DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(3)

                           DCBDPJ(1)  =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE1(J8,:)) &
                                      - COSTB*R2*DRIJDPJ(1)
                           DCBDPJ(2)  =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE2(J8,:)) &
                                      - COSTB*R2*DRIJDPJ(2)
                           DCBDPJ(3)  =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE3(J8,:)) &
                                      - COSTB*R2*DRIJDPJ(3)

                        ENDIF

                        ! calculate if I and J are both carbons 
                        IF (I <= NCARBON .AND. J <= NCARBON) THEN

                           ! calculate rho_cc
                           RHOCC   = RHOCC0 + RHOCC10*(COSTA + COSTB) + RHOCC20*(1.5D0*COSTA*COSTA & 
                                   + 1.5D0*COSTB*COSTB - 1.D0)
                           ! ENERGY1 is energy due to short-range anisotropic interactions
                           ! calculate vertical shift for first term
                           EXPFCT  = KKJ*DEXP(-ALPHACC*(ABSRIJ - RHOCC))
                           VSHIFT1 = KKJ*DEXP(-ALPHACC*(EWALDREALC - RHOCC))
                           ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                           ! ENERGY2 is energy due to damped dispersion
                           ! calculate vertical shift for second term
                           VSHIFT2 = DC6CC*DMPFCT_SHIFT/(EWALDREALC**6)
                           ENERGY2 = ENERGY2 - DC6CC*DMPFCT*R6 + VSHIFT2

                           IF (GTEST) THEN

                              ! DVDR is derivative of dispersion damping factor energy with factor 1/Rab
                              DVDR    = 6.D0*DC6CC*R6*R2*DMPFCT - DC6CC*R6*DDMPDR 
                              ! FRIJ is derivative of ENERGY1 wrt r_ij with factor 1/Rab
                              FRIJ(:) = ALPHACC*EXPFCT*(-NR(:) + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADR(:) &
                                      + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDR(:))
                              ! TIJ is derivative of ENERGY1 wrt pi with factor 1/Rab
                              TIJ(:)  = ALPHACC*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADPI(:) &
                                      + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPI(:))
                              ! TJI is derivative of ENERGY1 wrt pj with factor 1/Rab
                              TJI(:)  = ALPHACC*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCC10 + 3.D0*RHOCC20*COSTA)*DCADPJ(:) &
                                      + (RHOCC10 + 3.D0*RHOCC20*COSTB)*DCBDPJ(:)) 

                           ENDIF

                        ! calculate if I and J are both hydorgens
                        ELSEIF (I > NCARBON .AND. J > NCARBON) THEN

                           RHOHH  = RHOHH0 + RHOHH10*(COSTA + COSTB) + RHOHH20*(1.5D0*COSTA*COSTA      &
                                  + 1.5D0*COSTB*COSTB - 1.D0)
                           EXPFCT  = KKJ*DEXP(-ALPHAHH*(ABSRIJ - RHOHH))
                           VSHIFT1 = KKJ*DEXP(-ALPHAHH*(EWALDREALC - RHOHH))
                           ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                           VSHIFT2 = DC6HH*DMPFCT_SHIFT/(EWALDREALC**6)
                           ENERGY2 = ENERGY2 - DC6HH*DMPFCT*R6 + VSHIFT2

                           IF (GTEST) THEN

                              DVDR    = 6.D0*DC6HH*R6*R2*DMPFCT - DC6HH*R6*DDMPDR 
                              FRIJ(:) = ALPHAHH*EXPFCT*(-NR(:) + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADR(:) &
                                      + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHAHH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADPI(:) &
                                      + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHAHH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOHH10 + 3.D0*RHOHH20*COSTA)*DCADPJ(:) &
                                      + (RHOHH10 + 3.D0*RHOHH20*COSTB)*DCBDPJ(:))

                           ENDIF

                        ! calculate if I is carbon and J is hydrogen
                        ELSE IF (I <= NCARBON .AND. J > NCARBON) THEN 

                           RHOCH  = RHOCH0 + RHOC10H*COSTA + RHOCH10*COSTB + RHOC20H*(1.5D0*COSTA*COSTA &
                                  - 0.5D0) + RHOCH20*(1.5D0*COSTB*COSTB - 0.5D0)
                           EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                           VSHIFT1 = KKJ*DEXP(-ALPHACH*(EWALDREALC - RHOCH))
                           ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                           VSHIFT2 = DC6CH*DMPFCT_SHIFT/(EWALDREALC**6)
                           ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6 + VSHIFT2

                           IF (GTEST) THEN

                              DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                              FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADR(:) &
                                      + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADPI(:) &
                                      + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOC10H + 3.D0*RHOC20H*COSTA)*DCADPJ(:) &
                                      + (RHOCH10 + 3.D0*RHOCH20*COSTB)*DCBDPJ(:))

                           ENDIF

                        ELSE !IF(I > NCARBON .AND. J <= NCARBON) THEN

                           RHOCH  = RHOCH0 + RHOCH10*COSTA + RHOC10H*COSTB + RHOCH20*(1.5D0*COSTA*COSTA &
                                  - 0.5D0) + RHOC20H*(1.5D0*COSTB*COSTB - 0.5D0)
                           EXPFCT  = KKJ*DEXP(-ALPHACH*(ABSRIJ - RHOCH))
                           VSHIFT1 = KKJ*DEXP(-ALPHACH*(EWALDREALC - RHOCH))
                           ENERGY1 = ENERGY1 + EXPFCT - VSHIFT1
                           VSHIFT2 = DC6CH*DMPFCT_SHIFT/(EWALDREALC**6)
                           ENERGY2 = ENERGY2 - DC6CH*DMPFCT*R6 + VSHIFT2

                           IF (GTEST) THEN

                              DVDR    = 6.D0*DC6CH*R6*R2*DMPFCT - DC6CH*R6*DDMPDR 
                              FRIJ(:) = ALPHACH*EXPFCT*(-NR(:) + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADR(:) &
                                      + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDR(:))
                              TIJ(:)  = ALPHACH*EXPFCT*(-DRIJDPI(:)/ABSRIJ + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADPI(:) &
                                      + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPI(:))
                              TJI(:)  = ALPHACH*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + (RHOCH10 + 3.D0*RHOCH20*COSTA)*DCADPJ(:) &
                                      + (RHOC10H + 3.D0*RHOC20H*COSTB)*DCBDPJ(:))

                           ENDIF

                        ENDIF


                        IF (GTEST) THEN

                           ! total gradient wrt AA coords for rigid body J1
                           G(J5-2:J5) = G(J5-2:J5) + DVDR*DRIJDPI(:) + TIJ(:)
                           ! total gradient wrt AA coords for rigid body J2
                           G(J5-2:J5) = G(J5-2:J5) + DVDR*DRIJDPJ(:) + TJI(:)

                           ! gradients wrt cell parameters
                           if (boxderivt) then
                              do idx = 1, 6
                                    box_paramsgrad(idx) = box_paramsgrad(idx) + dot_product((dvdr*rss(1:3) + frij(1:3)), matmul(H_grad(:,:,idx), rcomfrac(:)))
                              enddo
                           endif ! box derivatives

                        ENDIF ! gtest
                        endif ! central box
                    endif ! within cutoff
                  enddo ! n
               enddo ! m
            enddo ! l
            enddo ! sites j
         enddo ! sites i
      enddo ! rigid bodies

      ! convert to cartesian coordinates
      XC(:) = 0.D0
      if (boxderivt) then
         xdum(:) = x(:)
         call cart2frac_rb_tri(nrigidbody, xdum, x, H_inverse)
      endif
      CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XC, X)
      ! restore cartesian rigid body coordinates
      if (boxderivt) x(:) = xdum(:)

      ! ENERGY3 and G3 are energy and gradient due to electrostatics
      ! computed using Ewald summation
      CALL EWALDSUM(1, XC, G3C, ENERGY3, GTEST)

      ! convert Ewald contribution of gradient to rigidbody coordinates
      IF (GTEST) G3(:) = 0.D0
      CALL TRANSFORMGRAD(G3C, X, G3)

      ! dj337: if input was cartesian, convert back to cartesian
      ! assumes ATOMRIGIDCOORDT is correct
      IF (ATOMRIGIDCOORDT) THEN

         ! convert to cartesian coordinates
         if (boxderivt) then
            xdum(:) = x(:)
            call cart2frac_rb_tri(nrigidbody, xdum, x, H_inverse)
         endif
         CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XR, X)
         X(:) = XR(:)
      ENDIF

      ! add WCA-style repulsion to keep cell volume away from zero
      if (boxderivt) call constrain_volume(v_fact, dv_fact, energy1, box_paramsgrad(4:6), gtest)

      ! sum energies / gradients and convert to kJ/mol
      ENERGY = (ENERGY1 + ENERGY2 + ENERGY3)*2625.499D0
      IF (GTEST) G(:) = (G(:) + G3(:))*2625.499D0
      if (gtest) box_paramsgrad(1:6) = box_paramsgrad(1:6)*2625.499D0

      END SUBROUTINE BENZGENRIGIDEWALD

!     ----------------------------------------------------------------------------------------------
!
!      SUBROUTINE DEFPAHARIGID()
!
!      USE COMMONS, ONLY: RHOCC0, RHOCC10, RHOCC20,  RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, RHOCH20, &
!                         ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, CCKJ
!
!      IMPLICIT NONE
! 
!      ALPHACC = 1.861500D0
!      ALPHAHH = 1.431200D0
!      ALPHACH = 1.775600D0
!
!      DC6CC    = 30.469D0
!      DC6HH    = 5.359D0
!      DC6CH    = 12.840D0
!
!      RHOCC0  = 5.814700D0
!      RHOCC10 = 0.021700D0
!      RHOCC20 =-0.220800D0
!
!      RHOHH0  = 4.486200D0
!      RHOHH10 =-0.271800D0
!      RHOHH20 = 0.0D0
!
!      RHOCH0  = 5.150500D0
!      RHOC10H = 0.021700D0
!      RHOCH10 =-0.271800D0
!      RHOC20H =-0.220800D0
!      RHOCH20 = 0.0D0
!
!      KKJ     = 1.D-03
!      CCKJ    = 1.D0   !1389.354848D0
!
!      END SUBROUTINE DEFPAHARIGID
!
!!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBENZENERIGIDEWALD()

      USE COMMONS, ONLY: NRBSITES, SITE, RBSTLA, STCHRG 
      USE GENRIGID, ONLY: NRIGIDBODY

      IMPLICIT NONE
 
      INTEGER :: J1

!!     C6H6

!!     D6h reference geometry: C-C: 1.397 angstrom; C-H: 1.087 angstrom
!      ! dj337: note that all done in atomic units

      SITE(1,:)  = (/ 2.63923430843701,   0.00000000000000,   0.00000000000000/)
      SITE(2,:)  = (/ 1.31961715421850,  -2.28564395764590,   0.00000000000000/)
      SITE(3,:)  = (/-1.31961715421850,  -2.28564395764590,   0.00000000000000/)
      SITE(4,:)  = (/-2.63923430843701,   0.00000000000000,   0.00000000000000/)
      SITE(5,:)  = (/-1.31961715421850,   2.28564395764590,   0.00000000000000/)
      SITE(6,:)  = (/ 1.31961715421850,   2.28564395764590,   0.00000000000000/)
      SITE(7,:)  = (/ 4.69338981379532,   0.00000000000000,   0.00000000000000/)
      SITE(8,:)  = (/ 2.34669490689766,  -4.06459480860986,   0.00000000000000/)
      SITE(9,:)  = (/-2.34669490689766,  -4.06459480860986,   0.00000000000000/)
      SITE(10,:) = (/-4.69338981379532,   0.00000000000000,   0.00000000000000/)
      SITE(11,:) = (/-2.34669490689766,   4.06459480860986,   0.00000000000000/)
      SITE(12,:) = (/ 2.34669490689766,   4.06459480860986,   0.00000000000000/)

      RBSTLA(1,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO H1
      RBSTLA(2,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO H2
      RBSTLA(3,:)  = SITE(9,:)  - SITE(3,:)                 ! Z FROM C3 TO H3
      RBSTLA(4,:)  = SITE(10,:) - SITE(4,:)                 ! Z FROM C4 TO H4
      RBSTLA(5,:)  = SITE(11,:) - SITE(5,:)                 ! Z FROM C5 TO H5
      RBSTLA(6,:)  = SITE(12,:) - SITE(6,:)                 ! Z FROM C6 TO H6
      RBSTLA(7,:)  = SITE(7,:)  - SITE(1,:)                 ! Z FROM C1 TO H1!
      RBSTLA(8,:)  = SITE(8,:)  - SITE(2,:)                 ! Z FROM C2 TO H2!
      RBSTLA(9,:)  = SITE(9,:) -  SITE(3,:)                 ! Z FROM C3 TO H3!
      RBSTLA(10,:) = SITE(10,:) - SITE(4,:)                 ! Z FROM C4 TO H4!
      RBSTLA(11,:) = SITE(11,:) - SITE(5,:)                 ! Z FROM C5 TO H5!
      RBSTLA(12,:) = SITE(12,:) - SITE(6,:)                 ! Z FROM C6 TO H6!

      DO J1 = 1, NRBSITES
 
         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      DO J1 = 1, NRIGIDBODY
         STCHRG(12*(J1-1)+1:12*(J1-1)+6) = -0.11114D0
         STCHRG(12*(J1-1)+7:12*(J1-1)+12) = 0.11114D0
      ENDDO

      END SUBROUTINE DEFBENZENERIGIDEWALD
