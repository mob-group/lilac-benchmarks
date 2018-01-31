!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!  ===============================================================
!  Multicomponent Lennard-Jones without a cutoff.
!  Assumed reduced units with SIGMA_AA = EPSILON_AA = 1.
!  Per-atom energy is stored in array VT(NATOMS).
!
!  ds656> 6/12/2014
!  ===============================================================
!
SUBROUTINE MLJ(X,GRAD,POT,GRADT,HESST,STRESST)
  !
  USE COMMONS, ONLY : NATOMS,VT,LABELS,NSPECIES,STRESS,MYNODE,GBHT
  USE MODHESS, ONLY : HESS
  USE POT_PARAMS, ONLY : MLJ_EPS, MLJ_SIG
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (IN) :: GRADT,HESST,STRESST
  DOUBLE PRECISION, INTENT (IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT (INOUT) :: POT, GRAD(3*NATOMS)
  !
  INTEGER :: I,J,J1,J13,J2,J23,T1,T2,IP
  DOUBLE PRECISION :: DIST2, DX(3), IR6, IR12, DUMMY, &
       SIG2(NSPECIES(0),NSPECIES(0)), &
       EPS4(NSPECIES(0),NSPECIES(0)), &
       EPS24(NSPECIES(0),NSPECIES(0)), IGRAD(3), DUDR, D2UDR2
  !
  ! Zero the potential and the gradient
  VT(1:NATOMS) = 0.0D0
  POT = 0.0D0
  IF(GRADT) GRAD(:) = 0.0D0
  IF(HESST) HESS(:,:) = 0.0D0
  IF(STRESST) STRESS(:,:,:) = 0.0D0
  !
  ! This can be precomputed elsewhere...
  SIG2(:,:) = MLJ_SIG(:,:)*MLJ_SIG(:,:)
  EPS4(:,:) = 4.0D0*MLJ_EPS(:,:)
  EPS24(:,:) = 24.0D0*MLJ_EPS(:,:)
  !
  IF(GBHT) THEN
     IP=1
  ELSE
     IP=MYNODE+1
  ENDIF
  !
  ! Double loop over atom types
  DO J1=1,NATOMS-1
     J13=3*(J1-1)
     T1=LABELS(J1,IP)
     DO J2=J1+1,NATOMS
        J23 = 3*(J2-1)
        T2=LABELS(J2,IP)
        !
        !
        DIST2 = 0.0D0
        DO I = 1,3
           DX(I) = X(J13 + I) - X(J23 + I)
           DIST2 = DIST2 + DX(I)*DX(I)
        ENDDO
        !
        ! ---- Potential-specific bit -------------
        IR6 = (SIG2(T1,T2)/DIST2)**3
        IR12 = IR6*IR6
        DUMMY = EPS4(T1,T2)*(IR12 - IR6)
        !WRITE(*,*) "DUMMY FOR POT", DUMMY
        ! -----------------------------------------
        !
        VT(J1) = VT(J1) + DUMMY
        VT(J2) = VT(J2) + DUMMY
        POT = POT + DUMMY
        !
        IF(GRADT) THEN ! Calculate gradient
           ! Compute 1/R*dUdR
           DUDR = EPS24(T1,T2)*(IR6 - 2.0D0*IR12)/DIST2
           !WRITE(*,*) "DUMMY FOR GRAD", DUMMY
           DO I = 1,3
              IGRAD(I) = DUDR*DX(I)
              GRAD(J13+I) = GRAD(J13+I) + IGRAD(I)
              GRAD(J23+I) = GRAD(J23+I) - IGRAD(I)
           ENDDO
           !
           IF(STRESST) THEN
              DO I=1,3
                 DO J=I,3
                    DUMMY = IGRAD(I)*DX(J)
                    STRESS(J1,I,J) = STRESS(J1,I,J) + &
                         DUMMY
                    STRESS(J2,I,J) = STRESS(J2,I,J) + &
                         DUMMY
                    STRESS(0,I,J) = STRESS(0,I,J) + &
                         DUMMY
                 ENDDO
              ENDDO
           ENDIF
           !
           IF(HESST) THEN
              ! Compute d2UdR2
              D2UDR2 = EPS24(T1,T2)*(26.0D0*IR12-7.0D0*IR6)/DIST2
              ! What we need is (D2UDR - DUDR/DIST) / DIST
              D2UDR2 = (D2UDR2 - DUDR) / DIST2
              !
              DO I=1,3
                 DO J=I,3
                    !
                    DUMMY = DX(I)*DX(J) ! Reset DUMMY
                    DUMMY = DUMMY*D2UDR2
                    IF(I==J) DUMMY = DUMMY + DUDR
                    !
                    ! Accumulate diagonal blocks first
                    HESS(J13+I, J13+J) = &
                         HESS(J13+I, J13+J) + DUMMY
                    HESS(J23+I, J23+J) = &
                         HESS(J23+I, J23+J) + DUMMY
                    ! NOTE: J < I will be filled later 
                    ! by symmetry of Hessian
                    !                                
                    ! Off-diagonal blocks acquire 
                    ! contribution from only one
                    ! pair (J1,J2) in pairwise-additive 
                    ! potentials, so no  
                    ! accumulation is needed 
                    ! (just value assignment).
                    !
                    DUMMY = -DUMMY
                    HESS(J13+I, J23+J) = DUMMY
                    HESS(J23+J, J13+I) = DUMMY 
                    ! imposed symmetry 
                    !
                    IF(I .NE. J) THEN ! Use symmetry of DUMMY
                       HESS(J13+J, J23+I) = DUMMY
                       HESS(J23+I, J13+J) = DUMMY
                    ENDIF
                    !
                 ENDDO
              ENDDO
              !
           ENDIF ! End of HESSIAN loop
           !
        ENDIF
        !
     ENDDO  ! Double loop over atoms
  ENDDO
  !
  IF(GRADT) THEN
     !
     ! Test per-atom energies and gradients:
     ! DO T1=1,NSPECIES(0)
     !    DO G1=1,2
     !       WRITE(*,'(A,3(1X,I5))') "BLJ_CLUST> T,G,NAT(T,G)=", &
     !            T1, G1, ATOMLISTS(T1,G1,0)
     !       DO J2=1,ATOMLISTS(T1,G1,0)
     !          J1=ATOMLISTS(T1,G1,J2)
     !          J13=3*(J1-1)
     !          WRITE(*,'(A,1X, I4, 4(1X,F12.6))') "BLJ_CLUST>", &
     !               J1, VT(J1), GRAD(J13+1), GRAD(J13+2), GRAD(J13+3)
     !       ENDDO
     !    ENDDO
     ! ENDDO
     !
     IF(STRESST) THEN
        DO I=1,3
           DO J=I,3
              STRESS(0,J,I) = STRESS(0,I,J) ! Impose symmetry
           ENDDO
        ENDDO
     ENDIF
     !
     IF(HESST) THEN
        ! Use symmetry to fill skipped entries in diagonal blocks
        DO J1=1,NATOMS
           DO I = 1,3 ! coordinates of atom 1
              DO J = I+1,3 ! coordinates of atom 2 
                 HESS(J1+J, J1+I) = HESS(J1+I, J1+J)
              ENDDO
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE MLJ
