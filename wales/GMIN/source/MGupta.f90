!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!   Loop structure recoded by J.A. Elliott 2009
!
!   GMIN is free software; you can redistribute it and/or modIFy
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
!   along with this program; IF not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!******************************************************************
!   Multi-component Gupta potential.
!   -ds656 (Dec 2014)
!******************************************************************
!
SUBROUTINE MGUPTA(X,GRAD,POT,GRADT,HESST,STRESST)
  !
  USE COMMONS, ONLY : NATOMS, VT, NSPECIES, ATOMLISTS,&
       INVATOMLISTS, STRESS, MYNODE, LABELS, GBHT
  USE MODHESS, ONLY : HESS
  USE POT_PARAMS, ONLY : MGUPTA_A, MGUPTA_XI, MGUPTA_P, &
       MGUPTA_Q, MGUPTA_R0, MGUPTA_M2Q, MGUPTA_XI_SQ, &
       MGUPTA_MP_DIVBY_R0, MGUPTA_M2Q_DIVBY_R0
  !
  IMPLICIT NONE
  !
  ! Passed variables -----------------------------------------
  LOGICAL, INTENT(IN) :: GRADT, HESST, STRESST
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(OUT) :: GRAD(3*NATOMS), POT
  !
  ! Local variables -----------------------------------------
  INTEGER :: I,J,J1,J2,J3,J13,J23,J33,T1,T2,T3,IP
  DOUBLE PRECISION :: DIST, DX(3), RHO(NATOMS), DUMMY, &
       RIJ(NATOMS,NATOMS), DUDR, D2UDR2, IGRAD(3), &
       REP_IJ(NATOMS,NATOMS), ATT_IJ(NATOMS,NATOMS), &
       REP_DUM1, REP_DUM2, REP_DUM3, ATT_DUM1, ATT_DUM2, &
       ATT_DUM3, RHO_DUM, RHO4HESS(NATOMS)
  !
  ! Do some initialisations
  RIJ(1:NATOMS,1:NATOMS)=0.0D0
  REP_IJ(1:NATOMS,1:NATOMS)=0.0D0
  ATT_IJ(1:NATOMS,1:NATOMS)=0.0D0
  RHO(1:NATOMS)=0.0D0
  RHO4HESS(1:NATOMS)=0.0D0
  VT(1:NATOMS) = 0.0D0
  !
  ! Atoms labels from global LABELS(1:NATOMS,IP)
  IF(GBHT) THEN
     IP=1
  ELSE
     IP=MYNODE+1
  ENDIF
  !
  DO J1=1,NATOMS-1
     J13 = 3*(J1-1)
     T1=LABELS(J1,IP)
     DO J2=J1+1,NATOMS
        J23 = 3*(J2-1)
        T2=LABELS(J2,IP)
        !
        DIST=0.0D0
        DO I=1,3
           DIST = DIST + (X(J23+I)-X(J13+I))**2
        ENDDO
        DIST = DSQRT(DIST)
        !
        RIJ(J2,J1)=DIST ! store distance
        RIJ(J1,J2)=DIST ! impose symmetry
        !
        DUMMY = MGUPTA_A(T1,T2)*&
             DEXP(MGUPTA_MP_DIVBY_R0(T1,T2)*DIST + &
             MGUPTA_P(T1,T2))
        REP_IJ(J1,J2) = DUMMY
        REP_IJ(J2,J1) = DUMMY
        !
        VT(J1) = VT(J1) + DUMMY
        VT(J2) = VT(J2) + DUMMY
        !
        DUMMY = MGUPTA_XI_SQ(T1,T2)*&
             DEXP(MGUPTA_M2Q_DIVBY_R0(T1,T2)*DIST - &
             MGUPTA_M2Q(T1,T2))
        ATT_IJ(J1,J2) = DUMMY
        ATT_IJ(J2,J1) = DUMMY
        !
        ! Accumulate many-body potential term:
        RHO(J1) = RHO(J1) + DUMMY
        RHO(J2) = RHO(J2) + DUMMY
        !
        IF (GRADT.OR.HESST.OR.STRESST) THEN
           !
           ! Precompute derivatives of rep and att here
           REP_IJ(J1,J2) = MGUPTA_MP_DIVBY_R0(T1,T2)*&
                REP_IJ(J1,J2)*2.0D0  ! need the 2 here!
           REP_IJ(J2,J1) = REP_IJ(J1,J2)
           ATT_IJ(J1,J2) = MGUPTA_M2Q_DIVBY_R0(T1,T2)*&
                ATT_IJ(J1,J2)
           ATT_IJ(J2,J1) = ATT_IJ(J1,J2)
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  ! Now store per-atom energies in array VT and sum over all atoms
  !
  POT=0.0D0
  DO J1=1,NATOMS
     RHO(J1)=DSQRT(RHO(J1)) ! square root density
     VT(J1) = VT(J1) - RHO(J1) ! accummulate per-atom energy
     POT = POT + VT(J1)       ! accumulate potential energy
     IF(GRADT.OR.GRADT.OR.STRESST) THEN
        RHO(J1)=1.0D0/RHO(J1) 
        IF(HESST) RHO4HESS(J1) = 0.25D0*RHO(J1)**3
        RHO(J1)=0.5D0*RHO(J1)
     ENDIF
  ENDDO
  !
  ! Calculate gradient terms, if required
  !
  IF (GRADT.OR.HESST.OR.STRESST) THEN
     ! 
     GRAD(:) = 0.0D0
     IF(HESST) HESS(:,:) = 0.0D0
     IF(STRESST) STRESS(:,:,:) = 0.0D0
     !
     DO J1=1,NATOMS-1
        J13 = 3*(J1-1)
        T1=LABELS(J1,IP)
        DO J2=J1+1,NATOMS
           J23 = 3*(J2-1)
           T2=LABELS(J2,IP)
           !
           RHO_DUM = RHO(J1) + RHO(J2) ! factor of 1/2 is included
           DUDR = REP_IJ(J1,J2) - RHO_DUM*ATT_IJ(J1,J2)
           ! Want 1/R*dU/dR
           DIST = RIJ(J1,J2)
           DUDR = DUDR/DIST
           !
           ! Cartesian components of distance
           DO I=1,3
              DX(I)=X(J13+I)-X(J23+I)
              IGRAD(I)=DUDR*DX(I) 
              GRAD(J13+I) = GRAD(J13+I) + IGRAD(I)
              GRAD(J23+I) = GRAD(J23+I) - IGRAD(I)
           ENDDO
           !
           IF(STRESST) THEN ! Accumulate local stresses
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
              !
              DUMMY=1.0D0/DIST**2 ! Reset DUMMY
              REP_DUM1=REP_IJ(J1,J2)/DIST
              REP_DUM2=MGUPTA_MP_DIVBY_R0(T1,T2)*&
                   REP_IJ(J1,J2)
              REP_DUM2=(REP_DUM2-REP_DUM1)*DUMMY
              ATT_DUM1=ATT_IJ(J1,J2)/DIST
              ATT_DUM2=MGUPTA_M2Q_DIVBY_R0(T1,T2)*&
                   ATT_IJ(J1,J2)
              ATT_DUM2=(ATT_DUM2-ATT_DUM1)*DUMMY
              !
              DO I = 1,3 ! coordinates of atom J1
                 DO J= 1,3 ! coordinate of atom J2 (/=J1)
                    !                                
                    DUMMY = DX(I)*DX(J) ! Reset DUMMY
                    !
                    REP_DUM3 = REP_DUM2*DUMMY
                    ATT_DUM3 = ATT_DUM2*DUMMY
                    !
                    IF (I == J) THEN
                       REP_DUM3 = REP_DUM3 + REP_DUM1
                       ATT_DUM3 = ATT_DUM3 + ATT_DUM1
                    ENDIF
                    !
                    D2UDR2 = REP_DUM3 - RHO_DUM*ATT_DUM3
                    DUDR = D2UDR2
                    !
                    ! Now go through triplets... FECK!
                    DO J3=1,NATOMS                    
                       J33 = 3*(J3-1)
                       T3 = LABELS(J3,IP)
                       !
                       IF ((J3.NE.J2).AND.(J3.NE.J1)) THEN
                          DUMMY = RHO4HESS(J3)*ATT_IJ(J1,J3)* &
                               ATT_IJ(J2,J3)/(RIJ(J1,J3)*RIJ(J2,J3))
                          D2UDR2 = D2UDR2 - &
                               DUMMY*(X(J13+I)-X(J33+I))* &
                               (X(J23+J)-X(J33+J))
                          DUDR = DUDR - &
                               DUMMY*(X(J23+I)-X(J33+I))* &
                               (X(J13+J)-X(J33+J))
                       ENDIF
                       !
                       IF (J3.NE.J2) THEN
                          DUMMY = RHO4HESS(J2)* &
                               ATT_IJ(J1,J2)*ATT_IJ(J3,J2)/ &
                               (RIJ(J1,J2)*RIJ(J3,J2))
                          D2UDR2 = D2UDR2 + &
                               DUMMY*(X(J13+I)-X(J23+I))* &
                               (X(J33+J)-X(J23+J))
                          DUDR = DUDR + &
                               DUMMY*(X(J33+I)-X(J23+I))* &
                               (X(J13+J)-X(J23+J))
                       ENDIF
                       !
                       IF (J3.NE.J1) THEN
                          DUMMY = RHO4HESS(J1)* &
                               ATT_IJ(J3,J1)*ATT_IJ(J2,J1)/ &
                               (RIJ(J3,J1)*RIJ(J2,J1))
                          D2UDR2 = D2UDR2 + &
                               DUMMY*(X(J33+I)-X(J13+I))* &
                               (X(J23+J)-X(J13+J))
                          DUDR = DUDR + &
                               DUMMY*(X(J23+I)-X(J13+I))* &
                               (X(J33+J)-X(J13+J))
                       ENDIF
                       !
                    ENDDO
                    !
                    ! Accumulate diagonal blocks first
                    HESS(J13+I, J13+J) = &
                         HESS(J13+I, J13+J) + D2UDR2
                    HESS(J23+I, J23+J) = &
                         HESS(J23+I, J23+J) + DUDR ! D2UDR2 /= DUDR
                    !
                    ! Compute (not accumulate!) off-diagonal blocks,
                    ! imposing symmetry over J1 and J2.
                    D2UDR2 = -D2UDR2
                    HESS(J13+I, J23+J) = D2UDR2
                    HESS(J23+J, J13+I) = D2UDR2
                    !
                 ENDDO
              ENDDO
              !
           ENDIF ! End of Hessian(J1,J2) block
           !
        ENDDO
     ENDDO
     !
     IF(STRESST) THEN
        DO I=1,3
           DO J=I,3
              STRESS(0,J,I) = STRESS(0,I,J) ! Impose symmetry
           ENDDO
        ENDDO
     ENDIF
     !
     ! IF(HESST) THEN
     !    write(*,*) 'Hessian 6x6 block:'
     !    DO J1=1,6
     !       WRITE(*,'(6(1X,F10.5))') (HESS(J1,J2), J2=1,6)
     !    ENDDO              
     ! ENDIF
     !
  ENDIF 
  !
  RETURN
  !
END SUBROUTINE MGUPTA
