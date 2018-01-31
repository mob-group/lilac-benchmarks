!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!----------------------------------------------------------------------------------------------!
!                                                                                              !
! PAP                                                                                          !
!                                                                                              !
! Calculates Energy and gradients for the patch-antipatch potential, using the rigid body      !
! angle axis framework described in:                                                           !
! 'Simulations of rigid bodies in an angle axis framework', Dwaipayan Chakrabarti and          !
! David Wales, Phys. Chem. Chem. Phys., 2009, 11, 1970-1976                                    !
! X: the positions and orientations of the bodies                                              !
! G: the gradients of the potential energy surface for changes in positions and orientations   !
! ENERGY: the potential energy for the configuration stored in X                               !
! GTEST: logical, true if gradients need to be calculated                                      !
! STEST: logical, true if second derivatives are required                                      !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE PAP(X, G, ENERGY, GTEST, STEST)

! NATOMS: twice the number of bodies, one for position and one for orientation
! NRBSITES: the number of patches and antipatches per body (plus one)
! RBSTLA: the directions of the patch and antipatch vectors
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA

! PAPALP, PAPEPS, PAPS, PAPCD: PAP potential parameters specified in the data file
      USE KEY, ONLY: PAPALP, PAPEPS, PAPS, PAPCD

! HESS: the hessian matrix
      USE MODHESS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) ::  G(3*NATOMS), ENERGY
      LOGICAL, INTENT(IN) :: GTEST, STEST

! I, J: iterators over bodies
! IT, IR: indices of translational and rotational coordinates of I
! JT, JR: indices of translational and rotational coordinates of J
! S, T: indices of patches in body I and J
! S2, T2: iterators over patches in I and J
! TLOW, TUP: lower and upper bounds for T2
! REALATOMS: the number of bodies
! OFFSET: the indicial offset to the start of orientations in X
      INTEGER :: I, IT, IR, J, JT, JR, S, S2, T, T2, TLOW, TUP, REALATOMS, OFFSET

! E: the orientations of patches and antipatches in the lab frame
! P: the rotation vector for a body
! RI, RJ: the position vecotr of a body
! ES, ET: the orientation vector of a patch
      DOUBLE PRECISION :: E(NATOMS*(NRBSITES-1)/2,3), P(3), RI(3), RJ(3), ES(3), ET(3)

! RMI: the rotation matrix for body I
! DRMIk: the derivative RMI with respect to the kth component of P
! D2RMIk: the second derivate of RMI with respect to the kth component of P
! D2RMIkl: the mixed second derivative of RMI with respect to the kth and lth components of P 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)

! DEk: matrix product of DRMIk with E, as in equation (12) of the paper
! D2Ek: matrix product of D2RMIk with E
! D2Ekl: matrix product of D2RMIkl with E
      DOUBLE PRECISION :: DE1(NATOMS*(NRBSITES-1)/2,3), DE2(NATOMS*(NRBSITES-1)/2,3), DE3(NATOMS*(NRBSITES-1)/2,3)
      DOUBLE PRECISION :: D2E1(NATOMS*(NRBSITES-1)/2,3), D2E2(NATOMS*(NRBSITES-1)/2,3), D2E3(NATOMS*(NRBSITES-1)/2,3)
      DOUBLE PRECISION :: D2E12(NATOMS*(NRBSITES-1)/2,3), D2E23(NATOMS*(NRBSITES-1)/2,3), D2E31(NATOMS*(NRBSITES-1)/2,3) 

! RIJ: vector displacement between bodies I and J
! NR: normalised RIJ
! RIJSQ: squared distance between I and J
! ABSRIJ: distance between I and J
! R2: ABSRIJ**-2
! R4: ABSRIJ**-4
      DOUBLE PRECISION :: RIJ(3), NR(3), RIJSQ, ABSRIJ, R2, R4 

! LJSIGMASQ: squared Lennard-Jones sigma parameter
! LJN: exponent in Lennard-Jones potential
! RLJN: (sigma over rij) to the power of LJN
! R2LJN: RLJN squared
! DUMMY: very short term store for a term
      DOUBLE PRECISION :: LJSIGMASQ, LJN, RLJN, R2LJN, DUMMY

! ANGFAC: factor in the angular potential phi, PI/(1-PAPCD)
! INVS: the reciprocal of the pap s parameter
! LAMBDA: longest distance with full PAP attraction
! PI: the mathematical constant
! DELR: ABSRIJ - LAMBDA
      DOUBLE PRECISION :: ANGFAC, INVS, LAMBDA, PI, DELR

! DVDR: derivative of LJ potential wrt rij, divided by rij
! D2VDR2: derivative of DVDR wrt rij, divided by rij
! WP: distance dependent part of PAP potential
! DWPDR: derivative of WP wrt rij divided by rij
! D2WPDR2: derivative of DWPDR wrt rij divided by rij
! COSWP, SINWP: terms in WP and its derivatives
! ALP, BET: dot product of ES, ET with NR
! ARGS, ARGT: argument of the Cos in angular part of PAP potential
! PHIS, PHIT: angluar parts of PAP potential
      DOUBLE PRECISION :: DVDR, D2VDR2, WP, DWPDR, D2WPDR2, SINWP, COSWP, ALP, BET, ARGS, ARGT, PHIS, PHIT

! DADR, DBDR: derivative of A, B wrt RI
! DPHISDR, DPHITDR: derivative of PHIS, PHIT wrt RI
! DPHISDR, DPHITDR: derivatives of DPHISDR, DPHITDR wrt RI
      DOUBLE PRECISION :: DADR(3), DBDR(3), DPHISDR, DPHITDR, D2PHISDR2, D2PHITDR2

! D2ADRk2, D2BDRk2: second derivative of ALP, BET wrt ri(k) twice
! D2ADRkl, D2BDRkl: mixed second derivative of ALP, BET wrt ri(k) and ri(l)
! DADPI, DBDPJ: derivative of ALP, BET wrt PI, PJ (cross terms are zero)
      DOUBLE PRECISION :: D2ADX2, D2ADY2, D2ADZ2, D2BDX2, D2BDY2, D2BDZ2
      DOUBLE PRECISION :: D2ADXY, D2ADYZ, D2ADZX, D2BDXY, D2BDYZ, D2BDZX
      DOUBLE PRECISION :: DADPI(3), DBDPJ(3)

! DPHIxDPy: derivative of PHI for patch x on body y wrt P for y (the cross terms are zero)
      DOUBLE PRECISION :: DPHISDPI(3)!, DPHITDPJ(3)

! Initialise potential parameters
      LJSIGMASQ = 1.D0 + 1.D0/PAPALP**(1.D0/3.D0)
      LJN       = 23
      PI        = 4.D0*DATAN(1.D0)
      ANGFAC    = PI/(1.D0 - PAPCD)
      INVS      = 1.D0/PAPS
      LAMBDA    = SQRT(1.D0 + (2.D0/PAPALP)**(1.D0/3.D0))

! Initialise energy, gradients and hessian
      ENERGY = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (STEST) HESS(:,:) = 0.D0

! Find REALATOMS as the actual number of bodies, and OFFSET as the start of rotational coordinates
      REALATOMS = NATOMS/2
      OFFSET    = 3*REALATOMS

! I Loop over all bodies  
      DO I = 1, REALATOMS

! Set IT as the index of the third translational coordinate of the Ith body, and IR as the index of the third rotational
! coordinate of the Ith body
! Set EI to be the rotation vector of body I
        IT = 3*I
        IR = OFFSET + IT
        P(:)  = X(IR-2:IR)

! Calculate the rotation matrix and derivatives thereof
        CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, STEST)

! S2 Loop over all patches and antipatches for body I
        DO S2 = 1, (NRBSITES-1)

! Set S as the index of the S2th patch of body I
          S = (NRBSITES-1)*(I-1) + S2
! Calculate E from the rotation matrix acting on the patch orientation in the body frame
          E(S,:) = MATMUL(RMI(:,:),RBSTLA(S2,:))

          IF (GTEST .OR. STEST) THEN
! If derivatives are required, calculate the derivatives of the rotation matrix for body I acting on the orientation of
! the Sth patch in the body frame
            DE1(S,:) = MATMUL(DRMI1(:,:),RBSTLA(S2,:))
            DE2(S,:) = MATMUL(DRMI2(:,:),RBSTLA(S2,:))
            DE3(S,:) = MATMUL(DRMI3(:,:),RBSTLA(S2,:))
          ENDIF ! End IF GTEST or STEST

          IF (STEST) THEN
! If second derivatives are required, calculate the second derivatives of the rotation matrix acting on the orientation
! of the Sth patch in the body frame
            D2E1(S,:) = MATMUL(D2RMI1(:,:),RBSTLA(S2,:))
            D2E2(S,:) = MATMUL(D2RMI2(:,:),RBSTLA(S2,:))
            D2E3(S,:) = MATMUL(D2RMI3(:,:),RBSTLA(S2,:))

            D2E12(S,:) = MATMUL(D2RMI12(:,:),RBSTLA(S2,:))
            D2E23(S,:) = MATMUL(D2RMI23(:,:),RBSTLA(S2,:))
            D2E31(S,:) = MATMUL(D2RMI31(:,:),RBSTLA(S2,:))
          ENDIF ! End IF STEST
          
        ENDDO ! End loop over patches
      ENDDO ! End loop over bodies

! I Loop over all bodies
      DO I = 1, REALATOMS  

! Set IT as the index of the third translational coordinate of I, IR as the index of the third rotational coordinate of I
! and RI as the position vector of I
        IT = 3*I
        IR = IT + OFFSET
        RI(:) = X(IT-2:IT)

! J Loop over all bodies
        DO J = 1, REALATOMS

! Bodies don't interact with themselves
          IF (I == J) CYCLE

! Set JT as the index of the third translational coordinate of J, JR as the index of the third rotational coordinate of J
! and RJ as the position vector of J
          JT = 3*J
          JR = JT + OFFSET
          RJ(:) = X(JT-2:JT)

! LJ contribution to energy and gradients
! Work out various useful powers of the distance from I to J
          RIJ(:) = RI(:) - RJ(:)
          RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
          ABSRIJ = SQRT(RIJSQ)
          NR(:) = RIJ(:)/ABSRIJ
          DELR   = ABSRIJ - LAMBDA
          R2 = 1.D0/RIJSQ
          R4 = R2*R2
! Construct the LJ terms
          RLJN = (LJSIGMASQ*R2)**(LJN/2.D0)
          R2LJN = RLJN*RLJN
! Energy contribution, only added half the time, to avoid double counting
          IF (J .GT. I) THEN
            ENERGY = ENERGY + 4.D0*(R2LJN-RLJN)
          ENDIF

! Gradient contribution
          IF (GTEST .OR. STEST) THEN
             DVDR = 4.D0*LJN*(RLJN - 2.D0*R2LJN)*R2
             G(IT-2:IT)  = G(IT-2:IT) + DVDR*RIJ(:)
!             G(JT-2:JT)  = G(JT-2:JT) - DVDR*RIJ(:) jwrm2> removed to avoid double counting
          ENDIF
! Hessian contributions
          IF (STEST) THEN
             D2VDR2 = 4.D0*LJN*((4.D0*LJN+4.D0)*R2LJN - (LJN+2.D0)*RLJN)*R4
! There is only a contribution when both derivatives are positional
! 3 completely diagonal terms, same molecule same coordinate
! xi, xi
             HESS(IT-2,IT-2) = HESS(IT-2,IT-2) + D2VDR2*RIJ(1)*RIJ(1) + DVDR
! yi, yi
             HESS(IT-1,IT-1) = HESS(IT-1,IT-1) + D2VDR2*RIJ(2)*RIJ(2) + DVDR
! zi, zi
             HESS(IT,IT)     = HESS(IT,IT)     + D2VDR2*RIJ(3)*RIJ(3) + DVDR
! 6 off diagonal terms on the diagonal block, same molecule different coordinate
! xi, yi
             DUMMY = D2VDR2*RIJ(1)*RIJ(2)
             HESS(IT-2,IT-1) = HESS(IT-2,IT-1) + DUMMY
             HESS(IT-1,IT-2) = HESS(IT-1,IT-2) + DUMMY
! yi, zi
             DUMMY = D2VDR2*RIJ(2)*RIJ(3)
             HESS(IT-1,IT) = HESS(IT-1,IT) + DUMMY
             HESS(IT,IT-1) = HESS(IT,IT-1) + DUMMY
! zi, xi
             DUMMY = D2VDR2*RIJ(3)*RIJ(1)
             HESS(IT,IT-2) = HESS(IT,IT-2) + DUMMY
             HESS(IT-2,IT) = HESS(IT-2,IT) + DUMMY
! 3 diagonal terms on off diagonal blocks, different molecule, same coordinate
! xi, xj
             HESS(IT-2,JT-2) = -D2VDR2*RIJ(1)*RIJ(1) - DVDR
! yi, yj
             HESS(IT-1,JT-1) = -D2VDR2*RIJ(2)*RIJ(2) - DVDR
! zi, zj
             HESS(IT,JT)     = -D2VDR2*RIJ(3)*RIJ(3) - DVDR
! 6 off diagonal terms on off diagonal blocks, different molecule different coordinate
! xi, yj
             DUMMY = -D2VDR2*RIJ(1)*RIJ(2)
             HESS(IT-2,JT-1) = HESS(IT-2,JT-1) + DUMMY
             HESS(JT-1,IT-2) = HESS(JT-1,IT-2) + DUMMY
! yi, zj
             DUMMY = -D2VDR2*RIJ(2)*RIJ(3)
             HESS(IT-1,JT) = HESS(IT-1,JT) + DUMMY
             HESS(JT,IT-1) = HESS(JT,IT-1) + DUMMY
! zi, xj
             DUMMY = -D2VDR2*RIJ(3)*RIJ(1)
             HESS(IT,JT-2) = HESS(IT,JT-2) + DUMMY
             HESS(JT-2,IT) = HESS(JT-2,IT) + DUMMY
          ENDIF

! Calculate the distance dependence of the PAP, same for all patches for a pair of bodies
          IF (DELR > INVS) THEN
! Distance is greater than LAMBDA+INVS, so no attraction
! No point proceeding with the calculation for this pair
            WP = 0.D0
            CYCLE
          ELSE IF (DELR < 0.D0) THEN
! Distance less than LAMBDA, so full attraction
            WP =-1.D0 
            IF (GTEST .OR. STEST) DWPDR = 0.D0
            IF (STEST) D2WPDR2 = 0.D0         
          ELSE
! Distance is between LAMBDA and LAMBDA+INVS, so potential varies from -1 to 0
            COSWP = COS(PI*DELR*PAPS)
            SINWP = SIN(PI*DELR*PAPS)
            WP =-0.5D0*(1.D0 + COSWP)
            IF (GTEST .OR. STEST) DWPDR = 0.5D0*PI*PAPS*SINWP/ABSRIJ
            IF (STEST) D2WPDR2 = 0.5D0*PI*PAPS*R2*(PI*PAPS*COSWP - SINWP/ABSRIJ) 
          ENDIF


! Loop over patches in I
          DO S2 = 1, (NRBSITES-1)

! Set bounds on T for patch-antipatch
            IF (S2 .LE. (NRBSITES-1)/2) THEN
              TLOW = (NRBSITES-1)/2 + 1
              TUP  = NRBSITES-1
            ELSE
              TLOW = 1
              TUP = (NRBSITES-1)/2
            ENDIF

! Loop over patches in J
            DO T2 = TLOW, TUP

! Find the indices of the third coordinate of patches S2 and T2
              S = (I-1)*(NRBSITES-1) + S2
              T = (J-1)*(NRBSITES-1) + T2
! Set ES and ET as the direction vectors for the patches
              ES(:) = E(S,:)
              ET(:) = E(T,:)
! Set ALP as ES dot nr and BET as ET dot nr
              ALP = -DOT_PRODUCT(NR(:), ES(:))
              BET =  DOT_PRODUCT(NR(:), ET(:))

! Calculate the values of PHI for the angular potential
              IF (ALP < PAPCD) THEN
! Angle greater than patch width, so no attraction
                PHIS =-1.D0
              ELSE
! Angle less than patch width, so some attraction
                ARGS = ANGFAC*(ALP-PAPCD)
                PHIS = -COS(ARGS)
              ENDIF 
              IF (BET < PAPCD) THEN
! Angle greater than patch width, so no attraction
                PHIT =-1.D0
              ELSE
! Angle less than patch width, so some attraction
                ARGT = ANGFAC*(BET-PAPCD)
                PHIT = -COS(ARGT)
              ENDIF 

! Add the energy contribution, only included half the time, to avoid double counting
              IF (J .GT. I) THEN
                ENERGY = ENERGY + 0.25D0*PAPEPS*(1.D0 + PHIS)*(1.D0 + PHIT)*WP
              ENDIF

! Now for the gradients
              IF (GTEST .OR. STEST) THEN

! Calculate the derivates of ALP and BET wrt RI
! For signs, remember that ALP already includes a minus
                DADR(:) = -ALP*RIJ(:)*R2 - ES(:)/ABSRIJ
                DBDR(:) = -BET*RIJ(:)*R2 + ET(:)/ABSRIJ

! Translational is only non zero if both PHIs are not -1
                IF ((ALP >= PAPCD) .AND. (BET >= PAPCD)) THEN
! Find the derivatives of the PHIs wrt a change in RI
                  DPHISDR = ANGFAC*SIN(ARGS)
                  DPHITDR = ANGFAC*SIN(ARGT)

! Add the contribution to the gradient of a change in RI, but not vice versa to avoid double counting
                  G(IT-2:IT) = G(IT-2:IT) + 0.25D0*PAPEPS*((1.D0 + PHIT)*WP*DPHISDR*DADR(:)        &
     &                   + (1.D0 + PHIS)*WP*DPHITDR*DBDR(:) + (1.D0+PHIS)*(1.D0+PHIT)*DWPDR*RIJ(:))
! jwrm2> removed to avoid double counting
!                  G(JT-2:JT) = G(JT-2:JT) - 0.25D0*PAPEPS*((1.D0 + PHIT)*WP*DPHISDR*DADR(:)       &
!     &                   + (1.D0 + PHIS)*WP*DPHITDR*DBDR(:) + (1.D0+PHIS)*(1.D0+PHIT)*DWPDR*RIJ(:))
                ELSE
                  DPHISDR = 0.D0
                  DPHITDR = 0.D0
                ENDIF ! End IF translational is non zero

! Rotational for I non zero if PHIS not -1
                IF (ALP >= PAPCD) THEN
                  DPHISDR = ANGFAC*SIN(ARGS)
! Find the derivatives of PHIS wrt a change in each of the rotational coordinates of body I
                  DPHISDPI(1) = -ANGFAC*SIN(ARGS)*DOT_PRODUCT(NR(:),DE1(S,:))
                  DPHISDPI(2) = -ANGFAC*SIN(ARGS)*DOT_PRODUCT(NR(:),DE2(S,:))
                  DPHISDPI(3) = -ANGFAC*SIN(ARGS)*DOT_PRODUCT(NR(:),DE3(S,:))
                  G(IR-2:IR) = G(IR-2:IR) + 0.25D0*PAPEPS*(1.D0 + PHIT)*WP*DPHISDPI(:)
                ELSE
                  DPHISDPI(:) = 0.D0
                ENDIF ! End IF rotational for I non zero
! Rotational for J non zero if PHIT not -1
                IF (BET >= PAPCD) THEN
                  DPHITDR = ANGFAC*SIN(ARGT)
! Find the derivatives of PHIT wrt a change in each of the rotational coordinates of body J
! jwrm2> removed to avoid double counting
!                  DPHITDPJ(1) = ANGFAC*SIN(ARGT)*DOT_PRODUCT(NR(:),DE1(T,:))
!                  DPHITDPJ(2) = ANGFAC*SIN(ARGT)*DOT_PRODUCT(NR(:),DE2(T,:))
!                  DPHITDPJ(3) = ANGFAC*SIN(ARGT)*DOT_PRODUCT(NR(:),DE3(T,:))
!                  G(JR-2:JR) = G(JR-2:JR) + 0.25D0*PAPEPS*(1.D0 + PHIS)*WP*DPHITDPJ(:)
                ELSE
!                  DPHITDPJ(:) = 0.D0
                ENDIF ! End IF rotational for J non zero

              ENDIF ! END IF GTEST or STEST

! And finally the hessian
              IF (STEST) THEN

! Factors in the second derivative of phi wrt ri
                IF (ALP >= PAPCD) THEN
                  D2PHISDR2 = ANGFAC*ANGFAC*COS(ARGS)
                ELSE
                  D2PHISDR2 = 0.D0
                ENDIF
                IF (BET >= PAPCD) THEN
                  D2PHITDR2 = ANGFAC*ANGFAC*COS(ARGT)
                ELSE
                  D2PHITDR2 = 0.D0
                ENDIF

! Second derivatives of ALP wrt ri
                D2ADX2 = 2.D0*ES(1)*RIJ(1)*R2/ABSRIJ + 3*ALP*RIJ(1)*RIJ(1)*R4 - ALP*R2
                D2ADY2 = 2.D0*ES(2)*RIJ(2)*R2/ABSRIJ + 3*ALP*RIJ(2)*RIJ(2)*R4 - ALP*R2
                D2ADZ2 = 2.D0*ES(3)*RIJ(3)*R2/ABSRIJ + 3*ALP*RIJ(3)*RIJ(3)*R4 - ALP*R2
                D2ADXY = R2*(ES(1)*RIJ(2)+ES(2)*RIJ(1))/ABSRIJ + 3*ALP*RIJ(1)*RIJ(2)*R4
                D2ADYZ = R2*(ES(2)*RIJ(3)+ES(3)*RIJ(2))/ABSRIJ + 3*ALP*RIJ(2)*RIJ(3)*R4
                D2ADZX = R2*(ES(3)*RIJ(1)+ES(1)*RIJ(3))/ABSRIJ + 3*ALP*RIJ(3)*RIJ(1)*R4
! Second derivatives of BET wrt ri
                D2BDX2 = - 2.D0*ET(1)*RIJ(1)*R2/ABSRIJ + 3*BET*RIJ(1)*RIJ(1)*R4 - BET*R2
                D2BDY2 = - 2.D0*ET(2)*RIJ(2)*R2/ABSRIJ + 3*BET*RIJ(2)*RIJ(2)*R4 - BET*R2
                D2BDZ2 = - 2.D0*ET(3)*RIJ(3)*R2/ABSRIJ + 3*BET*RIJ(3)*RIJ(3)*R4 - BET*R2
                D2BDXY = - R2*(ET(1)*RIJ(2)+ET(2)*RIJ(1))/ABSRIJ + 3*BET*RIJ(1)*RIJ(2)*R4
                D2BDYZ = - R2*(ET(2)*RIJ(3)+ET(3)*RIJ(2))/ABSRIJ + 3*BET*RIJ(2)*RIJ(3)*R4
                D2BDZX = - R2*(ET(3)*RIJ(1)+ET(1)*RIJ(3))/ABSRIJ + 3*BET*RIJ(3)*RIJ(1)*R4
! Derivatives of ALP wrt pi
                DADPI(1) = - DOT_PRODUCT(NR(:), DE1(S,:))
                DADPI(2) = - DOT_PRODUCT(NR(:), DE2(S,:))
                DADPI(3) = - DOT_PRODUCT(NR(:), DE3(S,:))
! Derivatives of BET wrt pj
                DBDPJ(1) = DOT_PRODUCT(NR(:), DE1(T,:))
                DBDPJ(2) = DOT_PRODUCT(NR(:), DE2(T,:))
                DBDPJ(3) = DOT_PRODUCT(NR(:), DE3(T,:))

! 6 completely diagonal terms, same molecule same coordinate
! xi, xi
                HESS(IT-2,IT-2) = HESS(IT-2,IT-2) + 0.25D0*PAPEPS*(                           &
                                  ((D2PHISDR2*DADR(1)*DADR(1) + DPHISDR*D2ADX2)*(1.D0+PHIT)   &
                                 + (D2PHITDR2*DBDR(1)*DBDR(1) + DPHITDR*D2BDX2)*(1.D0+PHIS)   &
                                 + 2.D0*DPHISDR*DADR(1)*DPHITDR*DBDR(1))*WP                   &
                                + (DPHISDR*DADR(1)*(1.D0+PHIT) + DPHITDR*DBDR(1)*(1.D0+PHIS)) &
                                  *2.D0*DWPDR*RIJ(1)                                          &
                                + (1.D0+PHIS)*(1.D0+PHIT)*(D2WPDR2*RIJ(1)*RIJ(1) + DWPDR))
! yi, yi
                HESS(IT-1,IT-1) = HESS(IT-1,IT-1) + 0.25D0*PAPEPS*(                           &
                                  ((D2PHISDR2*DADR(2)*DADR(2) + DPHISDR*D2ADY2)*(1.D0+PHIT)   &
                                 + (D2PHITDR2*DBDR(2)*DBDR(2) + DPHITDR*D2BDY2)*(1.D0+PHIS)   &
                                 + 2.D0*DPHISDR*DADR(2)*DPHITDR*DBDR(2))*WP                   &
                                + (DPHISDR*DADR(2)*(1.D0+PHIT) + DPHITDR*DBDR(2)*(1.D0+PHIS)) &
                                  *2.D0*DWPDR*RIJ(2)                                          &
                                + (1.D0+PHIS)*(1.D0+PHIT)*(D2WPDR2*RIJ(2)*RIJ(2) + DWPDR))

! zi, zi
                HESS(IT,IT)     = HESS(IT,IT)     + 0.25D0*PAPEPS*(                           &
                                  ((D2PHISDR2*DADR(3)*DADR(3) + DPHISDR*D2ADZ2)*(1.D0+PHIT)   &
                                 + (D2PHITDR2*DBDR(3)*DBDR(3) + DPHITDR*D2BDZ2)*(1.D0+PHIS)   &
                                 + 2.D0*DPHISDR*DADR(3)*DPHITDR*DBDR(3))*WP                   &
                                + (DPHISDR*DADR(3)*(1.D0+PHIT) + DPHITDR*DBDR(3)*(1.D0+PHIS)) &
                                  *2.D0*DWPDR*RIJ(3)                                          &
                                + (1.D0+PHIS)*(1.D0+PHIT)*(D2WPDR2*RIJ(3)*RIJ(3) + DWPDR))
! pi1, pi1
                HESS(IR-2,IR-2) = HESS(IR-2,IR-2) + 0.25D0*PAPEPS*WP*(1.D0+PHIT)*             &
                                  (D2PHISDR2*DADPI(1)*DADPI(1)                                &
                                   - DPHISDR*DOT_PRODUCT(NR(:),D2E1(S,:)))
! pi2, pi2
                HESS(IR-1,IR-1) = HESS(IR-1,IR-1) + 0.25D0*PAPEPS*WP*(1.D0+PHIT)*             &
                                  (D2PHISDR2*DADPI(2)*DADPI(2)                                &
                                   - DPHISDR*DOT_PRODUCT(NR(:),D2E2(S,:)))
! pi3, pi3
                HESS(IR,IR)     = HESS(IR,IR)     + 0.25D0*PAPEPS*WP*(1.D0+PHIT)*             &
                                  (D2PHISDR2*DADPI(3)*DADPI(3)                                &
                                   - DPHISDR*DOT_PRODUCT(NR(:),D2E3(S,:)))

! 30 off diagonal terms on the diagonal block, same molecule different coordinate
! xi, yi
                DUMMY = 0.25D0*PAPEPS*(                                                          &
                        ((D2PHISDR2*DADR(1)*DADR(2) + DPHISDR*D2ADXY)*(1.D0+PHIT)                &
                       + (D2PHITDR2*DBDR(1)*DBDR(2) + DPHITDR*D2BDXY)*(1.D0+PHIS)                &
                       + DPHISDR*DPHITDR*(DADR(1)*DBDR(2) + DADR(2)*DBDR(1)))*WP                 &
                      + (DPHISDR*DADR(1)*(1.D0+PHIT) + DPHITDR*DBDR(1)*(1.D0+PHIS))*DWPDR*RIJ(2) &
                      + (DPHISDR*DADR(2)*(1.D0+PHIT) + DPHITDR*DBDR(2)*(1.D0+PHIS))*DWPDR*RIJ(1) &
                      + (1.D0+PHIS)*(1.D0+PHIT)*D2WPDR2*RIJ(1)*RIJ(2))
                HESS(IT-2,IT-1) = HESS(IT-2,IT-1) + DUMMY
                HESS(IT-1,IT-2) = HESS(IT-1,IT-2) + DUMMY
! yi, zi
                DUMMY = 0.25D0*PAPEPS*(                                                          &
                        ((D2PHISDR2*DADR(2)*DADR(3) + DPHISDR*D2ADYZ)*(1.D0+PHIT)                &
                       + (D2PHITDR2*DBDR(2)*DBDR(3) + DPHITDR*D2BDYZ)*(1.D0+PHIS)                &
                       + DPHISDR*DPHITDR*(DADR(2)*DBDR(3) + DADR(3)*DBDR(2)))*WP                 &
                      + (DPHISDR*DADR(2)*(1.D0+PHIT) + DPHITDR*DBDR(2)*(1.D0+PHIS))*DWPDR*RIJ(3) &
                      + (DPHISDR*DADR(3)*(1.D0+PHIT) + DPHITDR*DBDR(3)*(1.D0+PHIS))*DWPDR*RIJ(2) &
                      + (1.D0+PHIS)*(1.D0+PHIT)*D2WPDR2*RIJ(2)*RIJ(3))
                HESS(IT-1,IT) = HESS(IT-1,IT) + DUMMY
                HESS(IT,IT-1) = HESS(IT,IT-1) + DUMMY
! zi, xi
                DUMMY = 0.25D0*PAPEPS*(                                                          &
                        ((D2PHISDR2*DADR(3)*DADR(1) + DPHISDR*D2ADZX)*(1.D0+PHIT)                &
                       + (D2PHITDR2*DBDR(3)*DBDR(1) + DPHITDR*D2BDZX)*(1.D0+PHIS)                &
                       + DPHISDR*DPHITDR*(DADR(3)*DBDR(1) + DADR(1)*DBDR(3)))*WP                 &
                      + (DPHISDR*DADR(3)*(1.D0+PHIT) + DPHITDR*DBDR(3)*(1.D0+PHIS))*DWPDR*RIJ(1) &
                      + (DPHISDR*DADR(1)*(1.D0+PHIT) + DPHITDR*DBDR(1)*(1.D0+PHIS))*DWPDR*RIJ(3) &
                      + (1.D0+PHIS)*(1.D0+PHIT)*D2WPDR2*RIJ(3)*RIJ(1))
                HESS(IT,IT-2) = HESS(IT,IT-2) + DUMMY
                HESS(IT-2,IT) = HESS(IT-2,IT) + DUMMY
! xi, pi1
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE1(S,1)/ABSRIJ - DADPI(1)*R2*RIJ(1))      &
                                           + D2PHISDR2*DADPI(1)*DADR(1))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(1)*DPHISDR*DADPI(1))                   &
                                       + DPHISDR*DADPI(1)*(1.D0 + PHIT)*DWPDR*RIJ(1))
                HESS(IT-2,IR-2) = HESS(IT-2,IR-2) + DUMMY
                HESS(IR-2,IT-2) = HESS(IR-2,IT-2) + DUMMY
! xi, pi2
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE2(S,1)/ABSRIJ - DADPI(2)*R2*RIJ(1))      &
                                           + D2PHISDR2*DADPI(2)*DADR(1))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(1)*DPHISDR*DADPI(2))                   &
                                       + DPHISDR*DADPI(2)*(1.D0 + PHIT)*DWPDR*RIJ(1))
                HESS(IT-2,IR-1) = HESS(IT-2,IR-1) + DUMMY
                HESS(IR-1,IT-2) = HESS(IR-1,IT-2) + DUMMY
! xi, pi3
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE3(S,1)/ABSRIJ - DADPI(3)*R2*RIJ(1))      &
                                           + D2PHISDR2*DADPI(3)*DADR(1))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(1)*DPHISDR*DADPI(3))                   &
                                       + DPHISDR*DADPI(3)*(1.D0 + PHIT)*DWPDR*RIJ(1))
                HESS(IT-2,IR) = HESS(IT-2,IR) + DUMMY
                HESS(IR,IT-2) = HESS(IR,IT-2) + DUMMY
! yi, pi1
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE1(S,2)/ABSRIJ - DADPI(1)*R2*RIJ(2))      &
                                           + D2PHISDR2*DADPI(1)*DADR(2))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(2)*DPHISDR*DADPI(1))                   &
                                       + DPHISDR*DADPI(1)*(1.D0 + PHIT)*DWPDR*RIJ(2))
                HESS(IT-1,IR-2) = HESS(IT-1,IR-2) + DUMMY
                HESS(IR-2,IT-1) = HESS(IR-2,IT-1) + DUMMY
! yi, pi2
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE2(S,2)/ABSRIJ - DADPI(2)*R2*RIJ(2))      &
                                           + D2PHISDR2*DADPI(2)*DADR(2))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(2)*DPHISDR*DADPI(2))                   &
                                       + DPHISDR*DADPI(2)*(1.D0 + PHIT)*DWPDR*RIJ(2))
                HESS(IT-1,IR-1) = HESS(IT-1,IR-1) + DUMMY
                HESS(IR-1,IT-1) = HESS(IR-1,IT-1) + DUMMY
! yi, pi3
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE3(S,2)/ABSRIJ - DADPI(3)*R2*RIJ(2))      &
                                           + D2PHISDR2*DADPI(3)*DADR(2))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(2)*DPHISDR*DADPI(3))                   &
                                       + DPHISDR*DADPI(3)*(1.D0 + PHIT)*DWPDR*RIJ(2))
                HESS(IT-1,IR) = HESS(IT-1,IR) + DUMMY
                HESS(IR,IT-1) = HESS(IR,IT-1) + DUMMY
! zi, pi1
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE1(S,3)/ABSRIJ - DADPI(1)*R2*RIJ(3))      &
                                           + D2PHISDR2*DADPI(1)*DADR(3))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(3)*DPHISDR*DADPI(1))                   &
                                       + DPHISDR*DADPI(1)*(1.D0 + PHIT)*DWPDR*RIJ(3))
                HESS(IT,IR-2) = HESS(IT,IR-2) + DUMMY
                HESS(IR-2,IT) = HESS(IR-2,IT) + DUMMY
! zi, pi2
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE2(S,3)/ABSRIJ - DADPI(2)*R2*RIJ(3))      &
                                           + D2PHISDR2*DADPI(2)*DADR(3))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(3)*DPHISDR*DADPI(2))                   &
                                       + DPHISDR*DADPI(2)*(1.D0 + PHIT)*DWPDR*RIJ(3))
                HESS(IT,IR-1) = HESS(IT,IR-1) + DUMMY
                HESS(IR-1,IT) = HESS(IR-1,IT) + DUMMY
! zi, pi3
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHISDR*(-DE3(S,3)/ABSRIJ - DADPI(3)*R2*RIJ(3))      &
                                           + D2PHISDR2*DADPI(3)*DADR(3))*(1.D0 + PHIT)           &
                                           + DPHITDR*DBDR(3)*DPHISDR*DADPI(3))                   &
                                       + DPHISDR*DADPI(3)*(1.D0 + PHIT)*DWPDR*RIJ(3))
                HESS(IT,IR) = HESS(IT,IR) + DUMMY
                HESS(IR,IT) = HESS(IR,IT) + DUMMY
! pi1, pi2
                DUMMY = 0.25D0*PAPEPS*WP*(1.D0+PHIT)*(D2PHISDR2*DADPI(1)*DADPI(2)                &
                         - DPHISDR*DOT_PRODUCT(NR(:),D2E12(S,:)))
                HESS(IR-2,IR-1) = HESS(IR-2,IR-1) + DUMMY
                HESS(IR-1,IR-2) = HESS(IR-1,IR-2) + DUMMY
! pi2, pi3
                DUMMY = 0.25D0*PAPEPS*WP*(1.D0+PHIT)*(D2PHISDR2*DADPI(2)*DADPI(3)                &
                         - DPHISDR*DOT_PRODUCT(NR(:),D2E23(S,:)))
                HESS(IR-1,IR) = HESS(IR-1,IR) + DUMMY
                HESS(IR,IR-1) = HESS(IR,IR-1) + DUMMY
! pi3, pi1
                DUMMY = 0.25D0*PAPEPS*WP*(1.D0+PHIT)*(D2PHISDR2*DADPI(3)*DADPI(1)                &
                         - DPHISDR*DOT_PRODUCT(NR(:),D2E31(S,:)))
                HESS(IR,IR-2) = HESS(IR,IR-2) + DUMMY
                HESS(IR-2,IR) = HESS(IR-2,IR) + DUMMY

! 6 diagonal terms on the off diagonal block, different molecule same coordinate
! xi, xj
                HESS(IT-2,JT-2) = HESS(IT-2,JT-2) - 0.25D0*PAPEPS*(                           &
                                  ((D2PHISDR2*DADR(1)*DADR(1) + DPHISDR*D2ADX2)*(1.D0+PHIT)   &
                                 + (D2PHITDR2*DBDR(1)*DBDR(1) + DPHITDR*D2BDX2)*(1.D0+PHIS)   &
                                 + 2.D0*DPHISDR*DADR(1)*DPHITDR*DBDR(1))*WP                   &
                                + (DPHISDR*DADR(1)*(1.D0+PHIT) + DPHITDR*DBDR(1)*(1.D0+PHIS)) &
                                  *2.D0*DWPDR*RIJ(1)                                          &
                                + (1.D0+PHIS)*(1.D0+PHIT)*(D2WPDR2*RIJ(1)*RIJ(1) + DWPDR))

! yi, yj
                HESS(IT-1,JT-1) = HESS(IT-1,JT-1) - 0.25D0*PAPEPS*(                           &
                                  ((D2PHISDR2*DADR(2)*DADR(2) + DPHISDR*D2ADY2)*(1.D0+PHIT)   &
                                 + (D2PHITDR2*DBDR(2)*DBDR(2) + DPHITDR*D2BDY2)*(1.D0+PHIS)   &
                                 + 2.D0*DPHISDR*DADR(2)*DPHITDR*DBDR(2))*WP                   &
                                + (DPHISDR*DADR(2)*(1.D0+PHIT) + DPHITDR*DBDR(2)*(1.D0+PHIS)) &
                                  *2.D0*DWPDR*RIJ(2)                                          &
                                + (1.D0+PHIS)*(1.D0+PHIT)*(D2WPDR2*RIJ(2)*RIJ(2) + DWPDR))
! zi, zj
                HESS(IT,JT)     = HESS(IT,JT) - 0.25D0*PAPEPS*(                               &
                                  ((D2PHISDR2*DADR(3)*DADR(3) + DPHISDR*D2ADZ2)*(1.D0+PHIT)   &
                                 + (D2PHITDR2*DBDR(3)*DBDR(3) + DPHITDR*D2BDZ2)*(1.D0+PHIS)   &
                                 + 2.D0*DPHISDR*DADR(3)*DPHITDR*DBDR(3))*WP                   &
                                + (DPHISDR*DADR(3)*(1.D0+PHIT) + DPHITDR*DBDR(3)*(1.D0+PHIS)) &
                                  *2.D0*DWPDR*RIJ(3)                                          &
                                + (1.D0+PHIS)*(1.D0+PHIT)*(D2WPDR2*RIJ(3)*RIJ(3) + DWPDR))
! pi1, pj1
                HESS(IR-2,JR-2) = HESS(IR-2,JR-2) +                                           &
                                   0.25D0*PAPEPS*WP*DPHISDR*DADPI(1)*DPHITDR*DBDPJ(1)
! pi2, pj2
                HESS(IR-1,JR-1) = HESS(IR-1,JR-1) +                                           &
                                   0.25D0*PAPEPS*WP*DPHISDR*DADPI(2)*DPHITDR*DBDPJ(2)
! pi3, pj3
                HESS(IR,JR)     = HESS(IR,JR) +                                               &
                                   0.25D0*PAPEPS*WP*DPHISDR*DADPI(3)*DPHITDR*DBDPJ(3)

! 30 completely off diagonal terms, different molecule different coordinate
! xi, yj
                DUMMY = - 0.25D0*PAPEPS*(                                                        &
                        ((D2PHISDR2*DADR(1)*DADR(2) + DPHISDR*D2ADXY)*(1.D0+PHIT)                &
                       + (D2PHITDR2*DBDR(1)*DBDR(2) + DPHITDR*D2BDXY)*(1.D0+PHIS)                &
                       + DPHISDR*DPHITDR*(DADR(1)*DBDR(2) + DADR(2)*DBDR(1)))*WP                 &
                      + (DPHISDR*DADR(1)*(1.D0+PHIT) + DPHITDR*DBDR(1)*(1.D0+PHIS))*DWPDR*RIJ(2) &
                      + (DPHISDR*DADR(2)*(1.D0+PHIT) + DPHITDR*DBDR(2)*(1.D0+PHIS))*DWPDR*RIJ(1) &
                      + (1.D0+PHIS)*(1.D0+PHIT)*D2WPDR2*RIJ(1)*RIJ(2))
                HESS(IT-2,JT-1) = HESS(IT-2,JT-1) + DUMMY
                HESS(JT-1,IT-2) = HESS(JT-1,IT-2) + DUMMY
! yi, zj
                DUMMY = - 0.25D0*PAPEPS*(                                                        &
                        ((D2PHISDR2*DADR(2)*DADR(3) + DPHISDR*D2ADYZ)*(1.D0+PHIT)                &
                       + (D2PHITDR2*DBDR(2)*DBDR(3) + DPHITDR*D2BDYZ)*(1.D0+PHIS)                &
                       + DPHISDR*DPHITDR*(DADR(2)*DBDR(3) + DADR(3)*DBDR(2)))*WP                 &
                      + (DPHISDR*DADR(2)*(1.D0+PHIT) + DPHITDR*DBDR(2)*(1.D0+PHIS))*DWPDR*RIJ(3) &
                      + (DPHISDR*DADR(3)*(1.D0+PHIT) + DPHITDR*DBDR(3)*(1.D0+PHIS))*DWPDR*RIJ(2) &
                      + (1.D0+PHIS)*(1.D0+PHIT)*D2WPDR2*RIJ(2)*RIJ(3))
                HESS(IT-1,JT) = HESS(IT-1,JT) + DUMMY
                HESS(JT,IT-1) = HESS(JT,IT-1) + DUMMY
! zi, xj
                DUMMY = - 0.25D0*PAPEPS*(                                                        &
                        ((D2PHISDR2*DADR(3)*DADR(1) + DPHISDR*D2ADZX)*(1.D0+PHIT)                &
                       + (D2PHITDR2*DBDR(3)*DBDR(1) + DPHITDR*D2BDZX)*(1.D0+PHIS)                &
                       + DPHISDR*DPHITDR*(DADR(3)*DBDR(1) + DADR(1)*DBDR(3)))*WP                 &
                      + (DPHISDR*DADR(3)*(1.D0+PHIT) + DPHITDR*DBDR(3)*(1.D0+PHIS))*DWPDR*RIJ(1) &
                      + (DPHISDR*DADR(1)*(1.D0+PHIT) + DPHITDR*DBDR(1)*(1.D0+PHIS))*DWPDR*RIJ(3) &
                      + (1.D0+PHIS)*(1.D0+PHIT)*D2WPDR2*RIJ(3)*RIJ(1))
                HESS(IT,JT-2) = HESS(IT,JT-2) + DUMMY
                HESS(JT-2,IT) = HESS(JT-2,IT) + DUMMY
! xi, pj1
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE1(T,1)/ABSRIJ - DBDPJ(1)*R2*RIJ(1))       &
                                           + D2PHITDR2*DBDPJ(1)*DBDR(1))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(1)*DPHITDR*DBDPJ(1))                   &
                                       + DPHITDR*DBDPJ(1)*(1.D0 + PHIS)*DWPDR*RIJ(1))
                HESS(IT-2,JR-2) = HESS(IT-2,JR-2) + DUMMY
                HESS(JR-2,IT-2) = HESS(JR-2,IT-2) + DUMMY
! xi, pj2
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE2(T,1)/ABSRIJ - DBDPJ(2)*R2*RIJ(1))       &
                                           + D2PHITDR2*DBDPJ(2)*DBDR(1))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(1)*DPHITDR*DBDPJ(2))                   &
                                       + DPHITDR*DBDPJ(2)*(1.D0 + PHIS)*DWPDR*RIJ(1))
                HESS(IT-2,JR-1) = HESS(IT-2,JR-1) + DUMMY
                HESS(JR-1,IT-2) = HESS(JR-1,IT-2) + DUMMY
! xi, pj3
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE3(T,1)/ABSRIJ - DBDPJ(3)*R2*RIJ(1))       &
                                           + D2PHITDR2*DBDPJ(3)*DBDR(1))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(1)*DPHITDR*DBDPJ(3))                   &
                                       + DPHITDR*DBDPJ(3)*(1.D0 + PHIS)*DWPDR*RIJ(1))
                HESS(IT-2,JR) = HESS(IT-2,JR) + DUMMY
                HESS(JR,IT-2) = HESS(JR,IT-2) + DUMMY
! yi, pj1
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE1(T,2)/ABSRIJ - DBDPJ(1)*R2*RIJ(2))       &
                                           + D2PHITDR2*DBDPJ(1)*DBDR(2))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(2)*DPHITDR*DBDPJ(1))                   &
                                       + DPHITDR*DBDPJ(1)*(1.D0 + PHIS)*DWPDR*RIJ(2))
                HESS(IT-1,JR-2) = HESS(IT-1,JR-2) + DUMMY
                HESS(JR-2,IT-1) = HESS(JR-2,IT-1) + DUMMY
! yi, pj2
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE2(T,2)/ABSRIJ - DBDPJ(2)*R2*RIJ(2))       &
                                           + D2PHITDR2*DBDPJ(2)*DBDR(2))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(2)*DPHITDR*DBDPJ(2))                   &
                                       + DPHITDR*DBDPJ(2)*(1.D0 + PHIS)*DWPDR*RIJ(2))
                HESS(IT-1,JR-1) = HESS(IT-1,JR-1) + DUMMY
                HESS(JR-1,IT-1) = HESS(JR-1,IT-1) + DUMMY
! yi, pj3
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE3(T,2)/ABSRIJ - DBDPJ(3)*R2*RIJ(2))       &
                                           + D2PHITDR2*DBDPJ(3)*DBDR(2))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(2)*DPHITDR*DBDPJ(3))                   &
                                       + DPHITDR*DBDPJ(3)*(1.D0 + PHIS)*DWPDR*RIJ(2))
                HESS(IT-1,JR) = HESS(IT-1,JR) + DUMMY
                HESS(JR,IT-1) = HESS(JR,IT-1) + DUMMY
! zi, pj1
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE1(T,3)/ABSRIJ - DBDPJ(1)*R2*RIJ(3))       &
                                           + D2PHITDR2*DBDPJ(1)*DBDR(3))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(3)*DPHITDR*DBDPJ(1))                   &
                                       + DPHITDR*DBDPJ(1)*(1.D0 + PHIS)*DWPDR*RIJ(3))
                HESS(IT,JR-2) = HESS(IT,JR-2) + DUMMY
                HESS(JR-2,IT) = HESS(JR-2,IT) + DUMMY
! zi, pj2
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE2(T,3)/ABSRIJ - DBDPJ(2)*R2*RIJ(3))       &
                                           + D2PHITDR2*DBDPJ(2)*DBDR(3))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(3)*DPHITDR*DBDPJ(2))                   &
                                       + DPHITDR*DBDPJ(2)*(1.D0 + PHIS)*DWPDR*RIJ(3))
                HESS(IT,JR-1) = HESS(IT,JR-1) + DUMMY
                HESS(JR-1,IT) = HESS(JR-1,IT) + DUMMY
! zi, pj3
                DUMMY = 0.25D0*PAPEPS*(WP*((DPHITDR*(DE3(T,3)/ABSRIJ - DBDPJ(3)*R2*RIJ(3))       &
                                           + D2PHITDR2*DBDPJ(3)*DBDR(3))*(1.D0 + PHIS)           &
                                           + DPHISDR*DADR(3)*DPHITDR*DBDPJ(3))                   &
                                       + DPHITDR*DBDPJ(3)*(1.D0 + PHIS)*DWPDR*RIJ(3))
                HESS(IT,JR) = HESS(IT,JR) + DUMMY
                HESS(JR,IT) = HESS(JR,IT) + DUMMY
! pi1, pj2
                DUMMY = 0.25D0*PAPEPS*WP*DPHISDR*DADPI(1)*DPHITDR*DBDPJ(2)
                HESS(IR-2,JR-1) = HESS(IR-2,JR-1) + DUMMY
                HESS(JR-1,IR-2) = HESS(JR-1,IR-2) + DUMMY
! pi2, pj3
                DUMMY = 0.25D0*PAPEPS*WP*DPHISDR*DADPI(2)*DPHITDR*DBDPJ(3)
                HESS(IR-1,JR) = HESS(IR-1,JR) + DUMMY
                HESS(JR,IR-1) = HESS(JR,IR-1) + DUMMY
! pi3, pj1
                DUMMY = 0.25D0*PAPEPS*WP*DPHISDR*DADPI(3)*DPHITDR*DBDPJ(1)
                HESS(IR,JR-2) = HESS(IR,JR-2) + DUMMY
                HESS(JR-2,IR) = HESS(JR-2,IR) + DUMMY
                                
              ENDIF ! End IF STEST
            ENDDO ! End loop over T2 patches
          ENDDO ! End loop over S2 patches
        ENDDO ! End loop over J bodies
      ENDDO ! End loop over I bodies

      END SUBROUTINE PAP
 
!----------------------------------------------------------------------------------------------!
!                                                                                              !
! DEFPAP                                                                                       !
!                                                                                              !
! Defines pap potential by specifying the positions of patches and antipatches                 !
! The first half of the vectors are patches, the second half anitpatches                       !
!                                                                                              !
!----------------------------------------------------------------------------------------------!
      SUBROUTINE DEFPAP()

! NRBSITES: the total number of patches and antipatches per body (plus one)
! RBSITE: positions of sites, used for distance metric
! RBSTLA: vectors specifying the directions of patches and antipatches
!         these vecotrs MUST be NORMALISED
! NATOMS: the number of linesof coordinates in the odata file
! PAPALP: the PAP distance parameter alpha
! PAPID: parameter from odata file which determines the patch distribution used
! NTSITES: the total number of sites on all particles
      USE COMMONS, ONLY: NRBSITES, RBSITE, RBSTLA, NATOMS
      USE KEY, ONLY: PAPALP, PAPID, NTSITES

      IMPLICIT NONE

! I: loop counter
      INTEGER :: I

! LFCT: the length factor used for setting the distance of the patches from the CoM
      DOUBLE PRECISION :: LFCT

      IF (PAPID == 4) THEN
! Create a tetrahedral arrangement of patches and antipatches
        NRBSITES = 5
        ALLOCATE(RBSTLA(NRBSITES-1,3))
        RBSTLA(1,:)= 1.D0/SQRT(3.D0)*(/  1.D0,  1.D0,  1.D0/)
        RBSTLA(2,:)= 1.D0/SQRT(3.D0)*(/ -1.D0, -1.D0,  1.D0/)
        RBSTLA(3,:)= 1.D0/SQRT(3.D0)*(/ -1.D0,  1.D0, -1.D0/)
        RBSTLA(4,:)= 1.D0/SQRT(3.D0)*(/  1.D0, -1.D0, -1.D0/)

      ELSE IF (PAPID == 3) THEN
        NRBSITES = 3
        ALLOCATE(RBSTLA(NRBSITES-1,3))
! Create a patch along the positive z axis and an antipatch along the negative z axis
!        RBSTLA(1,:) = (/ 0.D0, 0.D0, 1.D0/)
!        RBSTLA(2,:) = (/ 0.D0, 0.D0,-1.D0/)
! Bernal patch 1 only
!        RBSTLA(1,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
!        RBSTLA(2,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
! Bernal patch 2 only
!        RBSTLA(1,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
!        RBSTLA(2,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
! Bernal patch 3 only
        RBSTLA(1,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
        RBSTLA(2,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)

      ELSE IF (PAPID == 2) THEN
        NRBSITES = 5
        ALLOCATE(RBSTLA(NRBSITES-1,3))
! Bernal patches 1 and 2 only
        RBSTLA(1,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
        RBSTLA(2,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
        RBSTLA(3,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
        RBSTLA(4,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
! Bernal patches 2 and 3 only
!        RBSTLA(1,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
!        RBSTLA(2,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
!        RBSTLA(3,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
!        RBSTLA(4,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)
! Bernal patches 1 and 3 only
!        RBSTLA(1,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
!        RBSTLA(2,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
!        RBSTLA(3,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
!        RBSTLA(4,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)

      ELSE IF (PAPID == 1) THEN
        NRBSITES = 7
        ALLOCATE(RBSTLA(NRBSITES-1,3))
! Bernal spiral patches
        RBSTLA(1,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
        RBSTLA(2,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
        RBSTLA(3,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
        RBSTLA(4,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
        RBSTLA(5,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
        RBSTLA(6,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)

      ELSE
        PRINT *, 'PAPID not implemented. Implemented PAPIDs are 1 (Bernal 6 patch), ', &
     &           '2 (Bernal 4 patch), 3 (Bernal 2 patch) and 4 (tetrahedral).'  

      ENDIF

      ALLOCATE(RBSITE(NRBSITES,3))
      NTSITES = NATOMS*NRBSITES/2
! Fill RBSITE (used for distance metric) with CoM and patch sites
      RBSITE(1,:) = (/0.D0, 0.D0, 0.D0/)
! Calculate a third of the equilibrium body-body distance
      LFCT = SQRT(1.D0 + (2.D0/PAPALP)**(1.D0/3.D0))/3.D0
      DO I = 2, NRBSITES
        RBSITE(I,:) = RBSTLA(I-1,:)
! Renormalise so that the distance of a patch from the CoM is one third of the equilibrium 
! body-body distance
        RBSITE(I,:) = RBSITE(I,:)*LFCT
      ENDDO 

      END SUBROUTINE DEFPAP

!----------------------------------------------------------------------------------------------!
!                                                                                              !
! INERTIAPAP                                                                                   !
!                                                                                              !
! Returns moments of inertia and total mass for the purpose of frequency calculations. Since   !
! the patches are directions only, not additional particles, they have no mass. The particle   !
! is a sphere of unit mass (in reduced units) with a radius taken to be half the equilibrium   !
! pair distance.                                                                               !
!                                                                                              !
! RMI: rotation matrix for the body; not really relevant for PAP                               !
! KBLOCK: the moment of inertia tensor                                                         !
! TMASS: the total mass of the body                                                            !
!----------------------------------------------------------------------------------------------!

      SUBROUTINE INERTIAPAP(RMI, KBLOCK, TMASS)

! PAPALP: pap alpha parameter for calculating equilibrium pair distance
        USE KEY, ONLY: PAPALP

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN)  :: RMI(3,3)
        DOUBLE PRECISION, INTENT(OUT) :: TMASS, KBLOCK(3,3)
! LFCT: half the equilibrium pair distance squared
! MOMENT: scalar moment of inertia of a sphere about any axis
        DOUBLE PRECISION :: LFCT, MOMENT

! Set the total mass as unity
        TMASS  = 1.D0
! Set the value of LFCT to half the equilibrium pair distance squared
        LFCT = 0.25D0*(1.D0 + (2.D0/PAPALP)**(1.D0/3.D0))
! Set the value of MOMENT using the formula I= 0.2*M*R*R
        MOMENT = 0.2D0*TMASS*LFCT

! Set the moment of inertia tensor to zero.
! Off diagonal elements will remain zero
        KBLOCK(:,:) = 0.D0

! Set the diagonal elements to MOMENT
        KBLOCK(1,1) = MOMENT
        KBLOCK(2,2) = MOMENT
        KBLOCK(3,3) = MOMENT

      END SUBROUTINE INERTIAPAP
