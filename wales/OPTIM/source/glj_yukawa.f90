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
!-----------------------------------------------------------------------
! 
! Energy, gradient and Hessian for generalised Lennard-Jones + Yukawa.
! [See Mossa et al. Langmuir 20, 10756 (2004).]
!
! 31/5/2013 ds656
!
SUBROUTINE GLJYPOT(X,GRAD,POT,GTEST,HTEST)
  !
  USE COMMONS, ONLY : NATOMS, PARAM1, PARAM2, PARAM3, DEBUG
  USE MODHESS, ONLY : HESS
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT (IN) :: GTEST, HTEST
  DOUBLE PRECISION, INTENT (IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT (OUT) :: POT, GRAD(3*NATOMS)
  !
  INTEGER :: J1, J2, J13, J23, I, J
  DOUBLE PRECISION :: EPS, SIG, GLJ_EXP, YUK_A, YUK_XI
  DOUBLE PRECISION :: DIST, DIST2, IDIST, IDIST2, DUMMY
  DOUBLE PRECISION :: DX(3), GLJ_ATT, GLJ_REP, YUK_EXP, DVDR, D2VDR2
  DOUBLE PRECISION :: POTIJ, POTIJ_YUK, EPS4, FACT, FACT2, FACT3
  !
  EPS = 1.0D0
  SIG = 1.0D0
  GLJ_EXP = PARAM1 ! Exponent for GLJ
  YUK_A = PARAM2 ! Prefactor for Yukawa (in units of EPS)
  YUK_XI = PARAM3 ! Length scale for Yukawa (in units of SIG)
  !
  IF(DEBUG) THEN
     WRITE(*,'(3(A,1X,F10.3))') 'gljy> ALPHA=', GLJ_EXP, ' A=', YUK_A, &
          ' XI=', YUK_XI
  ENDIF
  !
  EPS4 = 4.0D0*EPS
  !
  POT=0.0D0
  !
  IF (GTEST .AND. .NOT. HTEST) THEN ! With gradient only (no Hessian)
     !
     IF(DEBUG) WRITE(*,*) 'gljy> Call with gradient and no Hessian'
     !
     GRAD(:) = 0.0D0
     FACT = EPS4*GLJ_EXP
     !
     DO J1=1,NATOMS-1
        J13=3*(J1-1)
        DO J2=J1+1,NATOMS
           J23=3*(J2-1)
           !
           DIST2 = 0.0D0
           DO I = 1,3
              DX(I) = X(J13 + I) - X(J23 + I)
              DIST2 = DIST2 + DX(I)*DX(I)
           ENDDO
           DIST = DSQRT(DIST2)
           !
           GLJ_ATT = (SIG/DIST)**GLJ_EXP
           GLJ_REP = GLJ_ATT*GLJ_ATT
           DUMMY = DIST /  YUK_XI
           YUK_EXP = YUK_A*DEXP(-DUMMY)
           !
           POTIJ = EPS4*(GLJ_REP - GLJ_ATT) + YUK_EXP/DUMMY
           !
           POT = POT + POTIJ
           !
           DVDR = - ( FACT*(2.0D0*GLJ_REP - GLJ_ATT) + &
                YUK_EXP*(1.0D0 + 1.0D0/DUMMY) ) / DIST 
           !
           !IF(DEBUG) WRITE(*,*) 'gljy> DVDR = ', DVDR
           DO I = 1,3
              DX(I) = DX(I) / DIST
              GRAD(J13+I) = GRAD(J13+I) + DVDR*DX(I)
              GRAD(J23+I) = GRAD(J23+I) - DVDR*DX(I)
           ENDDO
           !
        ENDDO
     ENDDO
     !
  ELSEIF (GTEST .AND. HTEST) THEN ! With gradient and Hessian
     !
     IF(DEBUG) WRITE(*,*) 'gljy> Call with gradient and Hessian'
     !
     GRAD(:) = 0.0D0
     HESS(:,:) = 0.0D0
     FACT = EPS4*GLJ_EXP
     FACT3 = 1.0D0 + GLJ_EXP
     FACT2 = FACT*FACT3
     FACT3 = 2.0D0*(1.0D0 + GLJ_EXP/FACT3)
     !
     DO J1=1,NATOMS-1
        J13=3*(J1-1)
        DO J2=J1+1,NATOMS
           J23=3*(J2-1)
           !
           DIST2 = 0.0D0
           DO I = 1,3
              DX(I) = X(J13 + I) - X(J23 + I)
              DIST2 = DIST2 + DX(I)*DX(I)
           ENDDO
           DIST = DSQRT(DIST2)
           !
           GLJ_ATT = (SIG/DIST)**GLJ_EXP
           GLJ_REP = GLJ_ATT*GLJ_ATT
           DUMMY = DIST /  YUK_XI ! Set DUMMY
           YUK_EXP = YUK_A*DEXP(-DUMMY)
           POTIJ_YUK = YUK_EXP / DUMMY
           !
           POTIJ = EPS4*(GLJ_REP - GLJ_ATT) + POTIJ_YUK
           !
           POT = POT + POTIJ
           !
           ! --- Gradient ---
           !
           DUMMY = DUMMY + 1 ! Update DUMMY
           DVDR = - ( FACT*(2.0D0*GLJ_REP - GLJ_ATT) + &
                POTIJ_YUK*DUMMY ) / DIST2
           ! NOTE: DVDR here is actually DVDR/DIST
           !
           DO I = 1,3
              GRAD(J13+I) = GRAD(J13+I) + DVDR*DX(I)
              GRAD(J23+I) = GRAD(J23+I) - DVDR*DX(I)
           ENDDO
           !
           ! --- Hessian ---
           !
           DUMMY = DUMMY*DUMMY + 1 ! Update DUMMY
           D2VDR2 = ( POTIJ_YUK*DUMMY + & ! Yukawa
                FACT2*(FACT3*GLJ_REP - GLJ_ATT) ) / DIST2 ! GLJ 
           !
           ! What we actually need is (D2VDR - DVDR/DIST) / DIST2
           D2VDR2 = (D2VDR2 - DVDR) / DIST2
           !
           DO I = 1,3 ! coordinates of atom 1
              DO J = I,3 ! coordinates of atom 2 (need only J >= I)
                 !
                 DUMMY = DX(I)*DX(J) ! Reset DUMMY
                 DUMMY = DUMMY*D2VDR2
                 IF (I == J) DUMMY = DUMMY + DVDR
                 !
                 ! Accumulate diagonal blocks first
                 HESS(J13+I, J13+J) = HESS(J13+I, J13+J) + DUMMY
                 HESS(J23+I, J23+J) = HESS(J23+I, J23+J) + DUMMY
                 ! NOTE: J < I will be filled later by symmetry of Hessian
                 !
                 ! Off-diagonal blocks acquire contribution from only one
                 ! pair (J1,J2) in pairwise-additive potentials, so no
                 ! accumulation is needed (just value assignment).
                 !
                 DUMMY = -DUMMY
                 HESS(J13+I, J23+J) = DUMMY 
                 HESS(J23+J, J13+I) = DUMMY ! impose symmetry
                 !
                 IF(I .NE. J) THEN ! Use symmetry of DUMMY w.r.t. I,J.
                    HESS(J13+J, J23+I) = DUMMY
                    HESS(J23+I, J13+J) = DUMMY
                 ENDIF
                 !
              ENDDO
           ENDDO           
           !
        ENDDO
     ENDDO
     !
     ! Use symmetry to fill the skipped entries in diagonal blocks
     DO J1=1,NATOMS
        DO I = 1,3 ! coordinates of atom 1
           DO J = I+1,3 ! coordinates of atom 2
              HESS(J1+J, J1+I) = HESS(J1+I, J1+J)
           ENDDO
        ENDDO
     ENDDO
     !
  ELSE ! Without gradient or Hessian
     !
     IF(DEBUG) WRITE(*,*) 'gljy> Call without gradient or Hessian.'
     !
     DO J1=1,NATOMS-1
        J13=3*(J1-1)
        DO J2=J1+1,NATOMS
           J23=3*(J2-1)
           !
           DIST2 = 0.0D0
           DO I = 1,3
              DX(I) = X(J13 + I) - X(J23 + I)
              DIST2 = DIST2 + DX(I)*DX(I)
           ENDDO
           DIST = DSQRT(DIST2)
           !
           GLJ_ATT = (SIG/DIST)**GLJ_EXP
           GLJ_REP = GLJ_ATT*GLJ_ATT
           DUMMY = DIST /  YUK_XI
           YUK_EXP = YUK_A*DEXP(-DUMMY)
           !
           POTIJ = EPS4*(GLJ_REP - GLJ_ATT) + YUK_EXP/DUMMY
           !
           POT = POT + POTIJ
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE GLJYPOT
